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

  xpcs::io::ImmReader reader(conf->getIMMFilePath().c_str());
  xpcs::filter::SparseFilter filter;

  int r = 0;
  if (frameFrom > 0) {
      reader.SkipFrames(frameFrom);
      r += frameFrom;
  }

  while (r <= frameTo) {
    struct xpcs::io::ImmBlock* data = reader.NextFrames();
    filter.Apply(data);
    r++;
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


  delete [] pixels_sum;
  delete [] frames_sum;
}
