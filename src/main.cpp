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
#include "funcs.h"

#include "Eigen/Dense"
#include "Eigen/SparseCore"

#include "gflags/gflags.h"

using namespace Eigen;
using namespace std;




int main(int argc, char** argv)
{

    if (argc < 2) {
        fprintf(stderr, "Missing argument - <hdf5 file> <imm file>\n");
        return -1;
    }

    printf("File %s\n", argv[1]);
    Configuration *conf = Configuration::instance();
    conf->init(argv[1]);

    int* dqmap = conf->getDQMap();
    int *sqmap = conf->getSQMap();

    int frames = conf->getFrameTodoCount();
    int frameFrom = conf->getFrameStartTodo();
    int frameTo = conf->getFrameEndTodo();

    int pixels = conf->getFrameWidth() * conf->getFrameHeight();
    int maxLevel = Corr::calculateLevelMax(frames, 4);
    vector<std::tuple<int,int> > delays_per_level = Corr::delaysPerLevel(frames, 4, maxLevel);

    printf("Frames %d\n", frames);
    printf("Taus %d\n", delays_per_level.size());

    IMM imm(argv[2], frameFrom, frameTo, -1);

    // MatrixXf G2(pixels, delays_per_level.size());
    // MatrixXf IP(pixels, delays_per_level.size());
    // MatrixXf IF(pixels, delays_per_level.size());

    // Funcs funcs;
    // if (imm.getIsSparse()) {
    //     printf("IMM file is sparse\n");
    //     SparseMatF mat = imm.getSparsePixelData();

    //     printf ("Frame read %d\n", mat.cols());
    //     VectorXf frameSum = funcs.frameSum(mat);
    //     H5Result::writeFrameSum(conf->getFilename(), "exchange", frameSum);        

    //     VectorXf pixelSum = funcs.pixelSum(mat);
    //     H5Result::writePixelSum(conf->getFilename(), "exchange", pixelSum);        

    //     Corr::multiTauVec(mat, G2, IP, IF);
        
    // } else {
    //     MatrixXf pixelData = imm.getPixelData();
    //     VectorXf pixelSum = funcs.pixelSum(pixelData);
    //     H5Result::writePixelSum(conf->getFilename(), "exchange", pixelSum);
    //     Corr::multiTauVec(pixelData, G2, IP, IF);
    // }

    // H5Result::write2DData(conf->getFilename(), "exchange", "G2", G2);
    // H5Result::write2DData(conf->getFilename(), "exchange", "IP", IP);
    // H5Result::write2DData(conf->getFilename(), "exchange", "IF", IF);

    // Corr::normalizeG2s(G2, IP, IF);

    // Corr::twoTimesVec(imm.getPixelData());

}
