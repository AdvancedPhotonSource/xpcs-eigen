/**

Copyright (c) 2016, UChicago Argonne, LLC. All rights reserved.

Copyright 2016. UChicago Argonne, LLC. This software was produced 
under U.S. Government contract DE-AC02-06CH11357 for Argonne National 
Laboratory (ANL), which is operated by UChicago Argonne, LLC for the 
U.S. Department of Energy. The U.S. Government has rights to use, 
reproduce, and distribute this software.  NEITHER THE GOVERNMENT NOR 
UChicago Argonne, LLC MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR a
ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.  If software is 
modified to produce derivative works, such modified software should 
be clearly marked, so as not to confuse it with the version available 
from ANL.

Additionally, redistribution and use in source and binary forms, with 
or without modification, are permitted provided that the following 
conditions are met:

    * Redistributions of source code must retain the above copyright 
      notice, this list of conditions and the following disclaimer. 

    * Redistributions in binary form must reproduce the above copyright 
      notice, this list of conditions and the following disclaimer in 
      the documentation and/or other materials provided with the 
      distribution. 

    * Neither the name of UChicago Argonne, LLC, Argonne National 
      Laboratory, ANL, the U.S. Government, nor the names of its 
      contributors may be used to endorse or promote products derived 
      from this software without specific prior written permission. 

THIS SOFTWARE IS PROVIDED BY UChicago Argonne, LLC AND CONTRIBUTORS 
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT 
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS 
FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL UChicago 
Argonne, LLC OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, 
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, 
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; 
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT 
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN 
ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
POSSIBILITY OF SUCH DAMAGE.

**/

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <stdint.h>
#include <memory>
#include <map>
#include <vector>
#include <iostream>

#include <sys/stat.h>

#include "hdf5.h"
#include "gflags/gflags.h"
#include "spdlog/spdlog.h"

#include "corr.h"
#include "xpcs/configuration.h"
#include "h5_result.h"
#include "benchmark.h"
#include "xpcs/io/reader.h"
#include "xpcs/io/imm.h"
#include "xpcs/io/ufxc.h"
#include "xpcs/io/rigaku.h"
#include "xpcs/io/hdf5.h"
#include "xpcs/filter/filter.h"
#include "xpcs/filter/sparse_filter.h"
#include "xpcs/filter/dense_filter.h"
#include "xpcs/filter/stride.h"
#include "xpcs/data_structure/dark_image.h"
#include "xpcs/data_structure/sparse_data.h"
#include "xpcs/data_structure/row.h"

using namespace std;
namespace spd = spdlog; 

DEFINE_bool(g2out, false, "Write intermediate output from G2 computation");
DEFINE_bool(darkout, false, "Write dark average and std-data");
DEFINE_bool(ufxc, false, "IF the file format is from ufxc photon counting detector.");
DEFINE_bool(rigaku, false, "IF the file format is from rigaku photon counting detector.");
DEFINE_bool(hdf5, false, "IF the file format is HDF5 file format.");
DEFINE_bool(frame_threading, false, "Run twotime with frame threading.");
DEFINE_bool(transposed, false, "HDF5 flag when data is stored in col-major order");
DEFINE_int32(frameout, false, "Number of post-processed frames to write out for debuggin.");
DEFINE_string(imm, "", "The path to IMM file. By default the file specified in HDF5 metadata is used");
DEFINE_string(inpath, "", "The path prefix to replace");
DEFINE_string(outpath, "", "The path prefix to replace with");
DEFINE_string(exchange, "", "The output result path");
DEFINE_string(entry, "", "The metadata path in HDF5 file");

int main(int argc, char** argv)
{

  if (argc < 2) {
      fprintf(stderr, "Please specify a HDF5 metadata file\n");
      return 1;
  }

  xpcs::Benchmark total("Total");
 
  gflags::ParseCommandLineFlags(&argc, &argv, true);

  auto console = spd::stdout_color_mt("console");

  std::string entry = "/xpcs";

  if (!FLAGS_entry.empty())
      entry = FLAGS_entry;

  console->info("H5 metadata path {}", entry.c_str());

  xpcs::Configuration *conf = xpcs::Configuration::instance();
  conf->init(argv[1], entry);
  
  if (!FLAGS_imm.empty())
      conf->setIMMFilePath(FLAGS_imm);

  if (!FLAGS_inpath.empty() and !FLAGS_outpath.empty())
  {
      std::string file = conf->getIMMFilePath();
      std::string::size_type pos = file.find(FLAGS_inpath);

      if (pos != std::string::npos)
      {
          file.replace(file.begin()+pos,
                       file.end()-(strlen(file.c_str()) - strlen(FLAGS_inpath.c_str())),
                       FLAGS_outpath.begin(), FLAGS_outpath.end());
      }

      conf->setIMMFilePath(file);
  }

  if (!FLAGS_exchange.empty())
  {
    conf->OutputPath(FLAGS_exchange);
  }

  console->info("Processing IMM file at path {}..", conf->getIMMFilePath().c_str());
  struct stat st;
  if(stat(conf->getIMMFilePath().c_str(), &st) == 0) {
    char prefix[] = {' ', 'K', 'M', 'G', 'T'};
    unsigned long size = st.st_size;
    int suffix = 0;
    while (size >= 1024) {
       size = size / 1024;
       suffix++;
    }

    console->info("File size {0} {1}bytes", suffix > 0 ? (float)st.st_size/ pow(1024.0, suffix) : st.st_size, prefix[suffix]);
  }

  if (conf->IsTwoTime() && conf->SmoothingMethod() == xpcs::SMOOTHING_UNKNONW) {
    cerr<<"Smoothing method is not valid"<<endl;
    exit(1);
  }

  int* dqmap = conf->getDQMap();
  int *sqmap = conf->getSQMap();

  int frames = conf->getFrameTodoCount();
  int frameFrom = conf->getFrameStartTodo();
  int frameTo = frameFrom + frames; //conf->getFrameEndTodo();
  int swindow = conf->getStaticWindowSize();
  int stride_factor = conf->FrameStride();
  int average_factor = conf->FrameAverage();

  console->info("Data frames={0} stride={1} average={2}", frames, stride_factor, average_factor);
  console->debug("Frames count={0}, from={1}, todo={2}", frames, frameFrom, frameTo);

  int pixels = conf->getFrameWidth() * conf->getFrameHeight();
  int maxLevel = xpcs::Corr::calculateLevelMax(frames, conf->DelaysPerLevel());
  vector<std::tuple<int,int> > delays_per_level = xpcs::Corr::delaysPerLevel(frames, conf->DelaysPerLevel(), maxLevel);
  
  /* Print block to debug delays_per_level aka tau values. 
  int maxLevel_tmp = xpcs::Corr::calculateLevelMax(512, conf->DelaysPerLevel());
  vector<std::tuple<int,int> > delays_per_level_tmp = xpcs::Corr::delaysPerLevel(512, conf->DelaysPerLevel(), maxLevel_tmp);
  for (int i = 0; i < delays_per_level_tmp.size(); i++) {
    std::tuple<int, int> d = delays_per_level_tmp[i];
    printf("%d - %d\n", std::get<0>(d), std::get<1>(d)); 
    
  }
  return 0;*/

  float* g2s = new float[pixels * delays_per_level.size()];
  float* ips = new float[pixels * delays_per_level.size()];
  float* ifs = new float[pixels * delays_per_level.size()];

  for (int i = 0; i < (pixels * delays_per_level.size()); i++) {
    g2s[i] = 0.0f;
    ips[i] = 0.0f;
    ifs[i] = 0.0f;
  }

  xpcs::io::Reader *reader = NULL; 
  xpcs::filter::Filter *filter = NULL;

  if (FLAGS_ufxc) {
    reader = new xpcs::io::Ufxc(conf->getIMMFilePath().c_str());
  } else if (FLAGS_rigaku) {
     filter = new xpcs::filter::SparseFilter();
    reader = new xpcs::io::Rigaku(conf->getIMMFilePath().c_str(), filter);
  } else if (FLAGS_hdf5) {
    if (FLAGS_transposed)
      reader = new xpcs::io::Hdf5(conf->getIMMFilePath().c_str(), true);
    else {
      reader = new xpcs::io::Hdf5(conf->getIMMFilePath().c_str(), false);
    }
  } else {
    reader = new xpcs::io::Imm(conf->getIMMFilePath().c_str());
  }
      
  xpcs::data_structure::DarkImage *dark_image = NULL;
  {
    xpcs::Benchmark benchmark("Loading data");

    int r = 0;

    if (!reader->compression()) 
    {
      int dark_s = conf->getDarkFrameStart();
      int dark_e = conf->getDarkFrameEnd();
      int darks = conf->getDarkFrames();

      if (dark_s != dark_e) 
      {
        struct xpcs::io::ImmBlock *data = reader->NextFrames(darks);
        dark_image = new xpcs::data_structure::DarkImage(data->value, darks, pixels, conf->getFlatField());
        r += darks;
      }
    }
    
    if (frameFrom > 0 && r < frameFrom) 
    {
      reader->SkipFrames(frameFrom - r);
      r += (frameFrom - r);
    }

    if (!FLAGS_rigaku)
    {
      if (reader->compression()) 
      {
        filter = new xpcs::filter::SparseFilter();
      }
      else 
      {
        filter = new xpcs::filter::DenseFilter(dark_image);
      }

      int read_in_count = stride_factor > 1 ? stride_factor : average_factor;
      if (stride_factor > 1 && average_factor > 1)
        read_in_count = stride_factor * average_factor;

      // The last frame outside the stride will be ignored. 
      int f = 0;
      while (f < frames) {
        struct xpcs::io::ImmBlock* data = reader->NextFrames(read_in_count);
        filter->Apply(data);
        f++;
      }
    }

    if (FLAGS_frameout > 0 && FLAGS_frameout < frames) 
    {
      xpcs::data_structure::SparseData *data = filter->Data();
      int fcount = FLAGS_frameout;
      int f = 0;

      float* data_out = new float[pixels * fcount];

      for (int i = 0; i < (pixels*fcount); i++)
        data_out[i] = 0.0f;

      for (int j = 0; j < pixels; j++) 
      {
        if (!data->Exists(j)) continue;

        xpcs::data_structure::Row *row = data->Pixel(j);
        for (int x = 0; x < row->indxPtr.size(); x++) 
        {
          int f = row->indxPtr[x];
          float v = row->valPtr[x];

          if (f >= fcount) break;

          data_out[f*pixels+j] = v;
        }
      }

      xpcs::H5Result::write3DData(conf->getFilename(), 
                        conf->OutputPath(), 
                        "frames_out", 
                        conf->getFrameHeight(), 
                        conf->getFrameWidth(),
                        fcount, 
                        data_out);
    }
  }

  float* frames_sum = filter->FramesSum();
  if (conf->IsNormalizedByFramesum()) 
  {
    xpcs::Benchmark benchmark("Normalize by Frame-sum took ");
    float sum_of_framesums = 0.0f;
    float framesums_mean = 0.0f;
    for (int i = 0 ; i < frames; i++) 
    {
      sum_of_framesums += frames_sum[i+frames];
    }
    framesums_mean = sum_of_framesums / frames;

    xpcs::data_structure::SparseData *data = filter->Data();
    for (int j = 0; j < pixels; j++) 
    {
      if (!data->Exists(j)) continue;

      xpcs::data_structure::Row *row = data->Pixel(j);
      for (int x = 0; x < row->indxPtr.size(); x++) 
      {
        int f = row->indxPtr[x];
        row->valPtr[x] = row->valPtr[x] / (frames_sum[f+frames] / framesums_mean);
      }
    }
  }
  
  float* pixels_sum = filter->PixelsSum();
  for (int i = 0 ; i < pixels; i++) 
  {
    pixels_sum[i] /= frames;
  }

  xpcs::H5Result::write2DData(conf->getFilename(), 
                        conf->OutputPath(), 
                        "pixelSum", 
                        conf->getFrameHeight(), 
                        conf->getFrameWidth(), 
                        pixels_sum);

  
  xpcs::H5Result::write2DData(conf->getFilename(), 
                              conf->OutputPath(), 
                              "frameSum", 
                              2, 
                              frames, 
                              frames_sum);

  float *partitions_mean = filter->PartitionsMean();
  float *partial_part_mean = filter->PartialPartitionsMean();
  int total_static_partns = conf->getTotalStaticPartitions();
  int partitions = (int) floor((double) frames/ swindow);
  int *pixels_per_sbin = conf->PixelsPerStaticBin();
  float norm_factor = conf->getNormFactor();

  float denominator = 1.0f;
  for (int i = 0; i < total_static_partns; i++) {
    for (int j = 0; j < partitions; j++) {
      denominator = (float)pixels_per_sbin[i] * swindow;
      partial_part_mean[j * total_static_partns + i] /= denominator;
    }
  }

  for (int i = 0; i < total_static_partns; i++) {
    denominator = (float)pixels_per_sbin[i] * frames;
    partitions_mean[i] /= denominator;
  }

  xpcs::H5Result::write1DData(conf->getFilename(), 
                        conf->OutputPath(), 
                        "partition-mean-total", 
                        total_static_partns,
                        partitions_mean);

  xpcs::H5Result::write2DData(conf->getFilename(), 
                        conf->OutputPath(), 
                        "partition-mean-partial", 
                        partitions, 
                        total_static_partns,
                        partial_part_mean);
                        
  xpcs::H5Result::write1DData(conf->getFilename(), 
                        conf->OutputPath(), 
                        "partition_norm_factor", 
                        1, 
                        &norm_factor);

  xpcs::H5Result::write2DData(conf->getFilename(), 
                              conf->OutputPath(), 
                              "timestamp_clock", 
                              2, 
                              conf->getRealFrameTodoCount(), 
                              filter->TimestampClock());

  xpcs::H5Result::write2DData(conf->getFilename(), 
                              conf->OutputPath(), 
                              "timestamp_tick", 
                              2, 
                              conf->getRealFrameTodoCount(), 
                              filter->TimestampTicks());

  float *tau = new float[delays_per_level.size()];
  for (int x = 0 ; x < delays_per_level.size(); x++)
  {   
    std::tuple<int, int> value = delays_per_level[x];
    tau[x] = std::get<1>(value);
  }

  if (!conf->IsTwoTime()) {
    xpcs::H5Result::write1DData(conf->getFilename(), 
                              conf->OutputPath(), 
                              "tau", 
                              (int)delays_per_level.size(), 
                              tau);
  }
  
  {
    if (conf->IsTwoTime()) {
      xpcs::Benchmark benchmark("Computing G2 TwoTimes");
      xpcs::Corr::twotime(filter->Data(), FLAGS_frame_threading);
    } else {

      {
        xpcs::Benchmark benchmark("Computing G2 MultiTau");
        xpcs::Corr::multiTau2(filter->Data(), g2s, ips, ifs);
      }
      
      {
        xpcs::Benchmark benchmark("Normalizing Data");
        Eigen::MatrixXf G2s = Eigen::Map<Eigen::MatrixXf>(g2s, pixels, delays_per_level.size());
        Eigen::MatrixXf IPs = Eigen::Map<Eigen::MatrixXf>(ips, pixels, delays_per_level.size());
        Eigen::MatrixXf IFs = Eigen::Map<Eigen::MatrixXf>(ifs, pixels, delays_per_level.size());

        xpcs::Corr::normalizeG2s(G2s, IPs, IFs);

        if (FLAGS_g2out) {
          xpcs::Benchmark b("Writing G2s, IPs and IFs");
          xpcs::H5Result::write2DData(conf->getFilename(), conf->OutputPath(), "G2", G2s);
          xpcs::H5Result::write2DData(conf->getFilename(), conf->OutputPath(), "IP", IPs);
          xpcs::H5Result::write2DData(conf->getFilename(), conf->OutputPath(), "IF", IFs);
        }

      }
    }
    
  }
  
  if (!reader->compression() && FLAGS_darkout) {
    xpcs::Benchmark b("Writing Dark average and std image");
    if (dark_image != NULL) {
      double* dark_avg = dark_image->dark_avg();
      double* dark_std = dark_image->dark_std();
      xpcs::H5Result::write2DData(conf->getFilename(), 
                                  conf->OutputPath(), 
                                  "DarkAvg", 
                                  conf->getFrameHeight(),
                                  conf->getFrameWidth(), 
                                  dark_avg);
      
      xpcs::H5Result::write2DData(conf->getFilename(), 
                                  conf->OutputPath(), 
                                  "DarkStd", 
                                  conf->getFrameHeight(),
                                  conf->getFrameWidth(),
                                  dark_std);
    }
  }
}

