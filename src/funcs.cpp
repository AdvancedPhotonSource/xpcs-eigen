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

#include "funcs.h"
#include "configuration.h"

#include <iostream>
#include <math.h>

using namespace Eigen;
using namespace std;

VectorXf Funcs::pixelSum(Eigen::Ref<Eigen::MatrixXf> pixelData) {
    return pixelData.rowwise().sum();
}

VectorXf Funcs::pixelSum(SparseMatF pixelData) {
    Configuration *conf = Configuration::instance();
    int fcount = conf->getFrameTodoCount();

    VectorXf psum(pixelData.rows());
    psum.setZero(pixelData.rows());

     for (int i = 0; i < pixelData.cols(); i++) {
        psum += pixelData.col(i);
    }

    return psum.array() / fcount;
} 

MatrixXf Funcs::pixelWindowSum(SparseMatF pixelData) {
    Configuration *conf = Configuration::instance();
    int fcount = conf->getFrameTodoCount();
    int swindow = conf->getStaticWindowSize();
    int totalStaticPartns = conf->getTotalStaticPartitions();
    int partitions = (int) ceil((double)fcount/swindow);

    MatrixXf pixelSums(pixelData.rows(), partitions+1);
    pixelSums.setZero(pixelSums.rows(), pixelSums.cols());
    int win = 1;
    for (int i = 0; i < pixelData.cols(); i++) {
        pixelSums.col(0) += pixelData.col(i);
        pixelSums.col(win) += pixelData.col(i);
        if ( i > 0 && (i % swindow) == 0)
          win++;
    }

    return pixelSums;
}

Eigen::MatrixXf Funcs::partitionMean(Eigen::Ref<Eigen::MatrixXf> pixelSum) 
{
    Configuration *conf = Configuration::instance();
    int fcount = conf->getFrameTodoCount();
    int swindow = conf->getStaticWindowSize();
    int totalStaticPartns = conf->getTotalStaticPartitions();
    int partitions = (int) ceil((double)fcount/swindow);

    float normFactor = conf->getNormFactor();

    map<int, map<int, vector<int>> > qbins = conf->getBinMaps();

    MatrixXf means(totalStaticPartns, partitions+1);

    int pixcount = 0;
    for (map<int, map<int, vector<int>> >::const_iterator it = qbins.begin(); 
            it != qbins.end(); it++) {

        int q = it->first;
        map<int, vector<int> > values =  it->second;
        
        for (map<int, vector<int>>::const_iterator it2 =  values.begin(); it2 != values.end(); it2++) {
            int sbin = it2->first;
            pixcount = 0;
            vector<int> pixels = it2->second;
            for(vector<int>::const_iterator pind = pixels.begin(); pind != pixels.end(); pind++) {
                int p = *pind;
                means.row(sbin-1) += pixelSum.row(p);
                pixcount++;
            }

            double denom = pixcount * swindow * normFactor;
            double tmp = means(sbin-1, 0);
            means.row(sbin-1).array() /= denom;
            denom = (double)pixcount * fcount * normFactor;
            means(sbin-1, 0) = tmp / denom;
        }
    }
    
    return means;
}

VectorXf Funcs::frameSum(SparseMatF pixelData) {
    Configuration *conf = Configuration::instance();

    VectorXf fsum(pixelData.cols());
    for (int i = 0; i < pixelData.cols(); i++) {
        fsum(i) = pixelData.col(i).sum() / conf->getDetEfficiency() / 
                    conf->getDetAdhuPhot() / conf->getDetPreset();
    }
    return fsum;
}

/**
 * Compute mask of pixles from the dqmap
 */
Eigen::MatrixXf Funcs::maskFromDQmap(int* dqmap, int w, int h) {
    for (int i = 1290; i < 1300; i++) {
        for (int j = 1335; j < 1340; j++) {
            printf (" %d ", dqmap[w*i + j]);
        }
        printf("\n");
    }
}