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

#include "corr.h"

#include "configuration.h"

#include <math.h>
#include <vector>

#include <iostream>

using namespace Eigen;
using namespace std;

void Corr::multiTau(const MatrixXf &pixelData, int pix) {
    int frames = pixelData.cols();
    int maxLevel = calculateLevelMax(frames, 4);

    Configuration* conf = Configuration::instance();

    int* dqmap = conf->getDQMap();
    int* sqmap = conf->getSQMap();

    int w = conf->getFrameWidth();
    int h = conf->getFrameHeight();

    int tau = 1;
    int level = 0;

    std::vector<double*> results;

    while (level <= maxLevel) {
        int tauIncrement = (int) pow(2.0, level);
        
        //TODO smooth intensities
        int dplCount = Corr::calculateDelayCount(4, level);

        for (int delayIndex = 0; delayIndex < dplCount; delayIndex++) {

            int lastframe = tau + tauIncrement;
            int fc = frames % 2 == 0 || level == 0 ? frames : frames - 1;

            if (lastframe > fc)
                break;

            double *result = Corr::computeG2Levels(pixelData, pix, frames, tau, level);
            results.push_back(result);

            tau += tauIncrement;
        }

        level++;
    }
}

void Corr::multiTauVec(Ref<MatrixXf> pixelData,
                       Ref<MatrixXf> G2, 
                       Ref<MatrixXf> IP,
                       Ref<MatrixXf> IF)
{
    int frames = pixelData.cols();
    int pixels = pixelData.rows();

    int maxLevel = calculateLevelMax(frames, 4);

    Configuration* conf = Configuration::instance();
    int w = conf->getFrameWidth();
    int h = conf->getFrameHeight();

    vector<std::tuple<int,int> > delays_per_level = delaysPerLevel(frames, 4, maxLevel);

    //TODO asserts for G2, IP and IF sizes

    int i = 0;
    int ll = 0; // Last level

    MatrixXf c0; 
    MatrixXf c1;

    int lastframe = frames;
    for (std::vector<std::tuple<int, int> >::iterator it = delays_per_level.begin() ;
                it != delays_per_level.end(); ++it)
    {

        // View of pixel data at two different tau values. 
        std::tuple<int, int> tau_level = *it;
        int level = std::get<0>(tau_level);
        int tau = std::get<1>(tau_level);

        if (ll != level)
        {

            // level change, smooth out intensities.
            if (lastframe % 2)
                lastframe -= 1;

            Map<MatrixXf, 0, OuterStride<> > t0(pixelData.data(), pixels, lastframe/2, OuterStride<>(2*pixels));
            Map<MatrixXf, 0, OuterStride<> > t1(pixelData.data() + pixels, pixels, lastframe/2, OuterStride<>(2*pixels));
            Map<MatrixXf> t2(pixelData.data(), pixels, lastframe/2);

            t2 = 0.5 * (t0 + t1);

            lastframe = lastframe / 2;
        }

        if (level > 0)
            tau = tau / pow(2, level) - 1;



        c0 = Map<MatrixXf>(pixelData.data(), pixels, lastframe-tau);
        c1 = Map<MatrixXf>(pixelData.data() + (tau*pixels), pixels, lastframe-tau);

        G2.col(i) = (c0.array() * c1.array()).rowwise().mean();
        IP.col(i) = c0.array().rowwise().mean();
        IF.col(i) = c1.array().rowwise().mean();

        i++;
        ll = level;
    }
}

void Corr::multiTauVec(SparseMatrix<float>& pixelData,
                       Ref<MatrixXf> G2, 
                       Ref<MatrixXf> IP,
                       Ref<MatrixXf> IF)
{
    int frames = pixelData.cols();
    int pixels = pixelData.rows();

    int maxLevel = calculateLevelMax(frames, 4);

    Configuration* conf = Configuration::instance();
    int w = conf->getFrameWidth();
    int h = conf->getFrameHeight();

    vector<std::tuple<int,int> > delays_per_level = delaysPerLevel(frames, 4, maxLevel);

    SparseMatrix<float> c0, c1;
    //TODO asserts for G2, IP and IF sizes
    int i = 0;
    int ll = 0; // Last level

    cout<<pixelData<<endl;
    int lastframe = frames;
    for (std::vector<std::tuple<int, int> >::iterator it = delays_per_level.begin() ;
                it != delays_per_level.end(); ++it)
    {

        // View of pixel data at two different tau values. 
        std::tuple<int, int> tau_level = *it;
        int level = std::get<0>(tau_level);
        int tau = std::get<1>(tau_level);

        if (ll != level)
        {

            // level change, smooth out intensities.
            if (lastframe % 2)
                lastframe -= 1;

            for (int i = 0; i < lastframe/2; i++) {
                pixelData.col(i) = 0.5 * (pixelData.col(i) + pixelData.col(i+1));
            }

            lastframe = lastframe / 2;
        }

        if (level > 0)
            tau = tau / pow(2, level) - 1;


        c0 = pixelData.middleCols(0, lastframe-tau);
        c1 = pixelData.middleCols(tau, lastframe-tau);

        G2.col(i) = c0.cwiseProduct(c1) * (VectorXf::Ones(c0.cols()) * 1.0/c0.cols());
        IP.col(i) = c0 * (VectorXf::Ones(c0.cols()) * 1.0/c0.cols());
        IF.col(i) = c1 * (VectorXf::Ones(c0.cols()) * 1.0/c0.cols());

        i++;
        ll = level;
    }   
}

//TODO: Refactor this function and possibly break into sub function for the unit-tests. 
void Corr::normalizeG2s( Ref<MatrixXf> G2,
                   Ref<MatrixXf> IP, 
                   Ref<MatrixXf> IF)
{

    Configuration* conf = Configuration::instance();
    int w = conf->getFrameWidth();
    int h = conf->getFrameHeight();

    map<int, map<int, vector<int>> > qbins = conf->getBinMaps();

    int totalStaticPartns = conf->getTotalStaticPartitions();
    int totalDynamicPartns = conf->getTotalDynamicPartitions();

    // Final result out of this function. 
    printf("qbin size - %d\n", qbins.size());
    MatrixXf g2(qbins.size(), G2.cols());
    MatrixXf stdError(qbins.size(), G2.cols());

    // A matrix of the form total_static partitions over tau values.
    // Each row is a static partition and column is different tau value.  
    MatrixXf g2Sums(totalStaticPartns, G2.cols());
    MatrixXf ipSums(totalStaticPartns, G2.cols());
    MatrixXf ifSums(totalStaticPartns, G2.cols());

    // Count of pixels for each static partition
    VectorXi partitionPixelCounts(totalStaticPartns);

    // Sum pixels for each static partitions. 
    for (map<int, map<int, vector<int>> >::const_iterator it = qbins.begin(); 
            it != qbins.end(); it++) {

        int q = it->first;
        map<int, vector<int> > values =  it->second;

        for (map<int, vector<int>>::const_iterator it2 =  values.begin(); it2 != values.end(); it2++) {
            int sbin = it2->first;

            vector<int> pixels = it2->second;

            for(vector<int>::const_iterator pind = pixels.begin(); pind != pixels.end(); pind++) {
                int p = *pind;
                
                g2Sums.row(sbin - 1) = g2Sums.row(sbin - 1) + G2.row(p);
                ipSums.row(sbin - 1) = ipSums.row(sbin - 1) + IP.row(p);
                ifSums.row(sbin - 1) = ifSums.row(sbin - 1) + IF.row(p);

                partitionPixelCounts(sbin-1) =+ 1;
            }
        }
    }

    // Compute averag of each static parition
    VectorXi sbinCounts(totalDynamicPartns);
    for (map<int, map<int, vector<int>> >::const_iterator it = qbins.begin(); 
            it != qbins.end(); it++) {

        int q = it->first;
        map<int, vector<int> > values =  it->second;

        for (map<int, vector<int>>::const_iterator it2 =  values.begin(); it2 != values.end(); it2++) {
            int sbin = it2->first;

            g2Sums.row(sbin-1) = g2Sums.row(sbin-1) / partitionPixelCounts(sbin-1);
            ipSums.row(sbin-1) = ipSums.row(sbin-1) / partitionPixelCounts(sbin-1);
            ifSums.row(sbin-1) = ifSums.row(sbin-1) / partitionPixelCounts(sbin-1);

            // // Normalize per bin static partition . 
            g2Sums.row(sbin-1) = g2Sums.row(sbin-1).array() / (ipSums.row(sbin-1).array() * ifSums.row(sbin-1).array());
        }
    }

    // Compute the mean of normalized g2 values. 
    for (map<int, map<int, vector<int>> >::const_iterator it = qbins.begin(); 
            it != qbins.end(); it++) {

        int q = it->first;
        map<int, vector<int> > values =  it->second;

        int count = 0;
        for (map<int, vector<int>>::const_iterator it2 =  values.begin(); it2 != values.end(); it2++) {
            int sbin = it2->first;
            g2.row(q - 1) += g2Sums.row(sbin-1);
            count++;
        }

        // Mean of each G2 across tau values. This is our final g2 Matrix. 
        g2.row(q - 1).array() = g2.row(q - 1).array() / count;
    }


    //TODO: return value + standard-error.  
}

double* Corr::computeG2Levels(const Eigen::MatrixXf &pixelData, 
                                int pixel,
                                int frameCount, 
                                int tau, 
                                int level) {

    int tauIncrement = (level == 0) ? 1 : (int) pow(2.0, level);
    double numerator = 0.0d;
    double sumPast = 0.0d;
    double sumFuture = 0.0d;

    int count = 0;

    int requiredFrames = tau + tauIncrement - 1;
    int remainingFrames = frameCount - 1;

    int frameIndex = 0;

    while (remainingFrames >= requiredFrames) {

        numerator += pixelData(pixel, frameIndex) * pixelData(pixel, frameIndex + tau);
        sumPast += pixelData(pixel, frameIndex);
        sumFuture += pixelData(pixel, frameIndex + tau);

        count++;

        remainingFrames -= tauIncrement;
        frameIndex += tauIncrement;
    }

    double *result = new double[3];

    result[0] = numerator/count;
    result[1] =  sumFuture/count;
    result[2] = sumPast/count;

    return result;
}

vector<tuple<int, int>> Corr::delaysPerLevel(int frameCount, int dpl, int maxDelay)
{
    vector< tuple<int, int> > result;

    int ll_dpl = 0;
    for (int i=0; i<= maxDelay; i++)
    {
        int step = (int) pow(2.0, i);

        int dpll = i == 0 ? dpl * 2 : dpl;

        for (int j=0; j<dpll; j++) {
            
            if ( (ll_dpl+step) >= frameCount) break;


            result.push_back(std::make_tuple(i, ll_dpl + step));
            ll_dpl = ll_dpl + step;
        }
    }   

    return result;
}

int Corr::calculateLevelMax(int frameCount, int dpl) {
    if (frameCount < dpl * 2) return 0;

    return (int) (floor(log2(frameCount) - log2(1.0 + 1.0/(double)(dpl))) - log2(dpl));
}

int Corr::calculateDelayCount(int dpl, int level) {
    return level == 0 ? (2 * dpl -1) : dpl;
}


