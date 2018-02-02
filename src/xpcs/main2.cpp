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

#include "hdf5.h"
#include "gflags/gflags.h"
#include "spdlog/spdlog.h"

#include "corr.h"
#include "xpcs/configuration.h"
#include "h5_result.h"
#include "benchmark.h"
#include "xpcs/io/imm_reader.h"
#include "xpcs/filter/sparse_filter.h"


using namespace std;
namespace spd = spdlog; 

DEFINE_bool(g2out, false, "Write intermediate output from G2 computation");
DEFINE_string(imm, "", "The path to IMM file. By default the file specified in HDF5 metadata is used");
DEFINE_string(inpath, "", "The path prefix to replace");
DEFINE_string(outpath, "", "The path prefix to replace with");
DEFINE_string(entry, "", "The metadata path in HDF5 file");

int main(int argc, char** argv)
{

  xpcs::Benchmark total("Total");
  gflags::ParseCommandLineFlags(&argc, &argv, true);

  if (argc < 2) {
      fprintf(stderr, "Please specify a HDF5 metadata file");
      return 1;
  }

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

  int* dqmap = conf->getDQMap();
  int *sqmap = conf->getSQMap();

  int frames = conf->getFrameTodoCount();
  int frameFrom = conf->getFrameStartTodo();
  int frameTo = conf->getFrameEndTodo();
  int swindow = conf->getStaticWindowSize();

  console->info("Data frames {0}..", frames);
  console->debug("Frames count={0}, from={1}, todo={2}", frames, frameFrom, frameTo);

  int pixels = conf->getFrameWidth() * conf->getFrameHeight();
  int maxLevel = xpcs::Corr::calculateLevelMax(frames, 4);
  vector<std::tuple<int,int> > delays_per_level = xpcs::Corr::delaysPerLevel(frames, 4, maxLevel);

  float* g2s = new float[pixels * delays_per_level.size()];
  float* ips = new float[pixels * delays_per_level.size()];
  float* ifs = new float[pixels * delays_per_level.size()];

  float *timestamp_clock = new float[2 * frames];
  float *timestamp_tick = new float[2 * frames];

  xpcs::filter::SparseFilter filter;
  xpcs::io::ImmReader reader(conf->getIMMFilePath().c_str());

  {
    xpcs::Benchmark benchmark("Loading data");
    int r = 0;
    if (frameFrom > 0) {
      reader.SkipFrames(frameFrom);
      r += frameFrom;
    }

    int f = 0;
    while (r <= frameTo) {
      struct xpcs::io::ImmBlock* data = reader.NextFrames();
      filter.Apply(data);
      timestamp_clock[f] = f + 1;
      timestamp_clock[f + frames] = data->clock[0];
      timestamp_tick[f] = f + 1;
      timestamp_tick[f + frames] = data->ticks[0]; 

      // delete [] data->index;
      // delete [] data->value;
      // delete data;
      r++;
      f++;
    }
  }
 
  float* pixels_sum = filter.PixelsSum();
  for (int i = 0 ; i < pixels; i++) {
    pixels_sum[i] /= frames;
  }

  float* frames_sum = filter.FramesSum();

  xpcs::H5Result::write2DData(conf->getFilename(), 
                        "exchange", 
                        "pixelSum", 
                        conf->getFrameHeight(), 
                        conf->getFrameWidth(), 
                        pixels_sum);

  xpcs::H5Result::write2DData(conf->getFilename(), 
                              "exchange", 
                              "frameSum", 
                              2, 
                              frames, 
                              frames_sum);

  float *partitions_mean = filter.PartitionsMean();
  float *partial_part_mean = filter.PartialPartitionsMean();
  int total_static_partns = conf->getTotalStaticPartitions();
  int partitions = (int) ceil((double) frames/ swindow);
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
                        "exchange", 
                        "partition-mean-total", 
                        total_static_partns,
                        partitions_mean);

  xpcs::H5Result::write2DData(conf->getFilename(), 
                        "exchange", 
                        "partition-mean-partial", 
                        partitions, 
                        total_static_partns,
                        partial_part_mean);
                        
  xpcs::H5Result::write1DData(conf->getFilename(), 
                        "exchange", 
                        "partition_norm_factor", 
                        1, 
                        &norm_factor);

  xpcs::H5Result::write2DData(conf->getFilename(), 
                              "exchange", 
                              "timestamp_clock", 
                              2, 
                              frames, 
                              timestamp_clock);

  xpcs::H5Result::write2DData(conf->getFilename(), 
                              "exchange", 
                              "timestamp_tick", 
                              2, 
                              frames, 
                              timestamp_tick);

  float *tau = new float[delays_per_level.size()];
  for (int x = 0 ; x < delays_per_level.size(); x++)
  {   
    std::tuple<int, int> value = delays_per_level[x];
    tau[x] = std::get<1>(value);
  }

  xpcs::H5Result::write1DData(conf->getFilename(), 
                              "exchange", 
                              "tau", 
                              (int)delays_per_level.size(), 
                              tau);

  {
    xpcs::Benchmark benchmark("Computing G2");
    xpcs::Corr::multiTau2(filter.Data(), g2s, ips, ifs);
  }

  
  xpcs::Benchmark benchmark("Normalizing Data");
  Eigen::MatrixXf G2s = Eigen::Map<Eigen::MatrixXf>(g2s, pixels, delays_per_level.size());
  Eigen::MatrixXf IPs = Eigen::Map<Eigen::MatrixXf>(ips, pixels, delays_per_level.size());
  Eigen::MatrixXf IFs = Eigen::Map<Eigen::MatrixXf>(ifs, pixels, delays_per_level.size());

  xpcs::Corr::normalizeG2s(G2s, IPs, IFs);

}
