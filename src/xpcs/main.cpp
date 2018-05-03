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
#include "xpcs/io/imm_reader.h"
#include "xpcs/filter/filter.h"
#include "xpcs/filter/sparse_filter.h"
#include "xpcs/filter/dense_filter.h"
#include "xpcs/filter/stride.h"
#include "xpcs/filter/average.h"
#include "xpcs/filter/dense_average.h"
#include "xpcs/data_structure/dark_image.h"


using namespace std;
namespace spd = spdlog; 

DEFINE_bool(g2out, false, "Write intermediate output from G2 computation");
DEFINE_string(imm, "", "The path to IMM file. By default the file specified in HDF5 metadata is used");
DEFINE_string(inpath, "", "The path prefix to replace");
DEFINE_string(outpath, "", "The path prefix to replace with");
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

  int* dqmap = conf->getDQMap();
  int *sqmap = conf->getSQMap();

  int frames = conf->getFrameTodoCount();
  int frameFrom = conf->getFrameStartTodo();
  int frameTo = conf->getFrameEndTodo();
  int swindow = conf->getStaticWindowSize();
  int stride_factor = conf->FrameStride();
  int average_factor = conf->FrameAverage();

  console->info("Data frames={0} stride={1} average={2}", frames, stride_factor, average_factor);
  console->debug("Frames count={0}, from={1}, todo={2}", frames, frameFrom, frameTo);

  int pixels = conf->getFrameWidth() * conf->getFrameHeight();
  int maxLevel = xpcs::Corr::calculateLevelMax(frames, conf->DelaysPerLevel());
  vector<std::tuple<int,int> > delays_per_level = xpcs::Corr::delaysPerLevel(frames, conf->DelaysPerLevel(), maxLevel);

  float* g2s = new float[pixels * delays_per_level.size()];
  float* ips = new float[pixels * delays_per_level.size()];
  float* ifs = new float[pixels * delays_per_level.size()];

  for (int i = 0; i < (pixels * delays_per_level.size()); i++) {
    g2s[i] = 0.0f;
    ips[i] = 0.0f;
    ifs[i] = 0.0f;
  }

  // printf("Effective frames = %d\n", frames);
  float *timestamp_clock = new float[2 * frames];
  float *timestamp_tick = new float[2 * frames];

  xpcs::io::ImmReader reader(conf->getIMMFilePath().c_str());
  xpcs::filter::Filter *filter = NULL;
  
  xpcs::data_structure::DarkImage *dark_image = NULL;
  {
    xpcs::Benchmark benchmark("Loading data");

    int r = 0;

    if (!reader.compression()) {
      int dark_s = conf->getDarkFrameStart();
      int dark_e = conf->getDarkFrameEnd();
      int darks = conf->getDarkFrames();

      if (dark_s != dark_e) {
        struct xpcs::io::ImmBlock *data = reader.NextFrames(darks);
        dark_image = new xpcs::data_structure::DarkImage(data->value, darks, pixels, conf->getFlatField());
        r += darks;
      }
    }
    
    if (frameFrom > 0 && r < frameFrom) {
      reader.SkipFrames(frameFrom - r);
      r += (frameFrom - r);
    }

    if (reader.compression())
      filter = new xpcs::filter::SparseFilter();
    else
      filter = new xpcs::filter::DenseFilter(dark_image);

    xpcs::filter::Stride stride;
    xpcs::filter::Average average;
    xpcs::filter::DenseAverage dense_average;

    int f = 0;
    int read_in_count = stride_factor > 1 ? stride_factor : average_factor;
    if (stride_factor > 1 && average_factor > 1)
      read_in_count = stride_factor * average_factor;

    // The last frame outside the stride will be ignored. 
    while (r <= ((frameTo+1) - read_in_count)) {
      struct xpcs::io::ImmBlock* data = reader.NextFrames(read_in_count);

      filter->Apply(data);
      timestamp_clock[f + frames] = data->clock[0];
      timestamp_tick[f + frames] = data->ticks[0]; 
      f++;
      r += read_in_count;
    }

  }
  
  float* pixels_sum = filter->PixelsSum();
  for (int i = 0 ; i < pixels; i++) {
    pixels_sum[i] /= frames;
  }

  float* frames_sum = filter->FramesSum();

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
      denominator = pixels_per_sbin[i] * swindow * norm_factor;
      partial_part_mean[j * total_static_partns + i] /= denominator;
    }
  }

  for (int i = 0; i < total_static_partns; i++) {
    denominator = pixels_per_sbin[i] * frames * norm_factor;
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

  xpcs::H5Result::write1DData(conf->getFilename(), 
                              conf->OutputPath(), 
                              "tau", 
                              (int)delays_per_level.size(), 
                              tau);

  {
    xpcs::Benchmark benchmark("Computing G2");
    xpcs::Corr::multiTau2(filter->Data(), g2s, ips, ifs);
  }

  
  xpcs::Benchmark benchmark("Normalizing Data");
  Eigen::MatrixXf G2s = Eigen::Map<Eigen::MatrixXf>(g2s, pixels, delays_per_level.size());
  Eigen::MatrixXf IPs = Eigen::Map<Eigen::MatrixXf>(ips, pixels, delays_per_level.size());
  Eigen::MatrixXf IFs = Eigen::Map<Eigen::MatrixXf>(ifs, pixels, delays_per_level.size());

  xpcs::Corr::normalizeG2s(G2s, IPs, IFs);

  if (FLAGS_g2out) {
    xpcs::Benchmark b("Writing G2s, IPs and IFs");
    xpcs::H5Result::write2DData(conf->getFilename(), conf->OutputPath(), "G2", pixels, delays_per_level.size(), g2s);
    xpcs::H5Result::write2DData(conf->getFilename(), conf->OutputPath(), "IP", pixels, delays_per_level.size(), ips);
    xpcs::H5Result::write2DData(conf->getFilename(), conf->OutputPath(), "IF", pixels, delays_per_level.size(), ifs);
  }
}

