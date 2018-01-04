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
#include "imm.h"
#include "corr.h"
#include "configuration.h"
#include "h5result.h"
#include "benchmark.h"

#include "gflags/gflags.h"
#include "spdlog/spdlog.h"

using namespace std;
namespace spd = spdlog; 

DEFINE_bool(g2out, false, "Write intermediate output from G2 computation");
DEFINE_string(imm, "", "The path to IMM file. By default the file specified in HDF5 metadata is used");
DEFINE_string(inpath, "", "The path prefix to replace");
DEFINE_string(outpath, "", "The path prefix to replace with");

int main(int argc, char** argv)
{

    gflags::ParseCommandLineFlags(&argc, &argv, true);

    if (argc < 2) {
        fprintf(stderr, "Please specify a HDF5 metadata file");
        return 1;
    }

    auto console = spd::stdout_color_mt("console");

    Eigen::initParallel();
    
    
    Configuration *conf = Configuration::instance();
    conf->init(argv[1]);

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

    console->info("Processing IMM file at path arg{}..", conf->getIMMFilePath().c_str());

    int* dqmap = conf->getDQMap();
    int *sqmap = conf->getSQMap();

    int frames = conf->getFrameTodoCount();
    int frameFrom = conf->getFrameStartTodo();
    int frameTo = conf->getFrameEndTodo();
    int swindow = conf->getStaticWindowSize();

    console->info("Data frames {0}..", frames);
    console->debug("Frames count={0}, from={1}, todo={2}", frames, frameFrom, frameTo);

    int pixels = conf->getFrameWidth() * conf->getFrameHeight();
    int maxLevel = Corr::calculateLevelMax(frames, 4);
    vector<std::tuple<int,int> > delays_per_level = Corr::delaysPerLevel(frames, 4, maxLevel);

    float* g2s = new float[pixels * delays_per_level.size()];
    float* ips = new float[pixels * delays_per_level.size()];
    float* ifs = new float[pixels * delays_per_level.size()];

    IMM imm(conf->getIMMFilePath().c_str(), frameFrom, frameTo, -1);

    MatrixXf G2(pixels, delays_per_level.size());
    MatrixXf IP(pixels, delays_per_level.size());
    MatrixXf IF(pixels, delays_per_level.size());

    {
        Benchmark benchmark("Writing frame-sums, pixel-sums, partition-means");
        float* fsum = imm.getFrameSums();
        H5Result::write2DData(conf->getFilename(), "exchange", "frameSum", 2, frames, fsum);
        float* psum = imm.getPixelSums();
        for (int i = 0 ; i < pixels; i++) {
            psum[i] /= frames;
        }

        H5Result::write2DData(conf->getFilename(), "exchange", "pixelSum", conf->getFrameHeight(), conf->getFrameWidth(), psum);
        int totalStaticPartns = conf->getTotalStaticPartitions();
        float* totalPartmean = imm.getTotalPartitionMean();
        float* partialPartmean = imm.getPartialPartitionMean();
        int partitions = (int) ceil((double)frames/swindow);

        H5Result::write1DData(conf->getFilename(), "exchange", "partition-mean-total", totalStaticPartns, totalPartmean);
        H5Result::write2DData(conf->getFilename(), "exchange", "partition-mean-partial", partitions, totalStaticPartns, partialPartmean);
    }

    {
        Benchmark b("Writing timestamps and taus");
        float* tclock = imm.getTimestampClock();
        float* ttick = imm.getTimestampTick();

        H5Result::write2DData(conf->getFilename(), "exchange", "timestamp_clock", 2, frames, tclock);
        H5Result::write2DData(conf->getFilename(), "exchange", "timestamp_tick", 2, frames, ttick);

        float *tau = new float[delays_per_level.size()];
        for (int x = 0 ; x < delays_per_level.size(); x++)
        {   
            std::tuple<int, int> value = delays_per_level[x];
            tau[x] = std::get<1>(value);
        }
        H5Result::write1DData(conf->getFilename(), "exchange", "tau", (int)delays_per_level.size(), tau);
    }
        
    {
        Benchmark benchmark("Computing G2");
        Corr::multiTau2(imm.getSparseData(), g2s, ips, ifs);
    }

    Benchmark benchmark("Normalizing G2");
    MatrixXf G2s = Map<MatrixXf>(g2s, pixels, delays_per_level.size());
    MatrixXf IPs = Map<MatrixXf>(ips, pixels, delays_per_level.size());
    MatrixXf IFs = Map<MatrixXf>(ifs, pixels, delays_per_level.size());

    Corr::normalizeG2s(G2s, IPs, IFs);
    
    if (FLAGS_g2out) {

        Benchmark b("Writing G2s, IPs and IFs");
        H5Result::write2DData(conf->getFilename(), "exchange", "G2", pixels, delays_per_level.size(), g2s);
        H5Result::write2DData(conf->getFilename(), "exchange", "IP", pixels, delays_per_level.size(), ips);
        H5Result::write2DData(conf->getFilename(), "exchange", "IF", pixels, delays_per_level.size(), ifs);
    }

}
