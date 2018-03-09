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

#include <math.h>
#include <vector>
#include <stdio.h>
#include <iostream>

#include <omp.h>

#include "configuration.h"
#include "h5_result.h"

#include "xpcs/data_structure/sparse_data.h"


namespace xpcs {

using std::vector;
using std::string;
using std::tuple;
using std::get;
using std::min;
using std::map;

using Eigen::VectorXi;
using Eigen::MatrixXf;
using Eigen::SparseMatrix;
using Eigen::VectorXf;
using Eigen::RowMajor;
using Eigen::Ref;
using Eigen::Map;
using Eigen::OuterStride;

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

    vector<double*> results;

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
                       Ref<MatrixXf> IF) {
    int frames = pixelData.cols();
    int pixels = pixelData.rows();

    int maxLevel = calculateLevelMax(frames, 4);

    Configuration* conf = Configuration::instance();
    int w = conf->getFrameWidth();
    int h = conf->getFrameHeight();

    vector<tuple<int,int> > delays_per_level = Corr::delaysPerLevel(frames, 4, maxLevel);

    //TODO asserts for G2, IP and IF sizes
    int i = 0;
    int ll = 0; // Last level

    ::Eigen::MatrixXf c0; 
    ::Eigen::MatrixXf c1;

    int lastframe = frames;
    for (vector<tuple<int, int> >::iterator it = delays_per_level.begin() ;
                it != delays_per_level.end(); ++it)
    {

        // View of pixel data at two different tau values. 
        tuple<int, int> tau_level = *it;
        int level = get<0>(tau_level);
        int tau = get<1>(tau_level);

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

void Corr::multiTauVec(SparseRMatF& pixelData,
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
    int fcount = conf->getFrameTodoCount();

    vector<tuple<int,int> > delays_per_level = delaysPerLevel(frames, 4, maxLevel);

    // SparseMatrix<float, RowMajor> c0, c1;

    //TODO asserts for G2, IP and IF sizes
    int i = 0;
    int ll = 0;

    int lastframe = frames;
    for (vector<tuple<int, int> >::iterator it = delays_per_level.begin() ;
                it != delays_per_level.end(); ++it)
    {

        // View of pixel data at two different tau values. 
        tuple<int, int> tau_level = *it;
        int level = get<0>(tau_level);
        int tau = get<1>(tau_level);

        if (ll != level)
        {

            // level change, smooth out intensities.
            if (lastframe % 2)
                lastframe -= 1;

            lastframe = lastframe / 2;

            printf("Last frame %d\n", lastframe);
            // int ij = 0;

            // for (int k = 0; k < lastframe; k+=2) {
            //     pixelData.row(ij) = 0.5 * (pixelData.row(k) + pixelData.row(k+1));
            //     ij++;
            // }
            float* vals = pixelData.valuePtr();
            int *ind = pixelData.innerIndexPtr();
            int *outer = pixelData.outerIndexPtr();

            for (int r = 0; r < pixels; r++)
            {
                int nonzeros = outer[r+1] - outer[r];

                if (!nonzeros) continue;

                int *iptr = ind + outer[r];
                float *vptr = vals + outer[r];

                int *ciptr = iptr;
                float *cvptr = vptr;

                int index = 0;
                int cnt = 0;
                int i0, i1;
                i0 = *iptr/2;

                if (nonzeros == 1) {
                    iptr[index] = i0;
                    vptr[index] = (0.0 + *vptr)/2;
                    continue;
                }

                int inc;
                while (i0 < lastframe && cnt < (nonzeros-1))
                {   
                    i0 = *iptr/2;
                    i1 = *(iptr+1)/2;
                    
                    inc = 0;
                    if (i0 == i1)
                    {
                        ciptr[index] = i0;
                        cvptr[index] = (*vptr + *(vptr+1)) / 2.0f;
                        inc = 2;
                    }
                    else
                    {
                        ciptr[index] = i0;
                        cvptr[index] = (*vptr + 0.0f) / 2.0f;
                        inc = 1;
                    }

                    iptr += inc; vptr += inc; cnt += inc;
                    index++; 
                }

                // In case,we are left with one off element
                if (i0 != i1 && i1 < lastframe && cnt < nonzeros)
                {
                    ciptr[index] = i1;
                    cvptr[index] = (0.0f + *vptr) / 2.0f;
                    cnt++;
                } 

            }


        }

        if (level > 0)
            tau = tau / pow(2, level);

        // cout<<pixelData.row(20).nonZeros()<<endl;
        // Get lastframe-tau number of cols starting at 0 for c0 and tau for c1. 
        const SparseMatrix<float, RowMajor> c0 = pixelData.block(0, 0, pixels, lastframe-tau);
        const SparseMatrix<float, RowMajor> c1 = pixelData.block(0, tau, pixels, lastframe-tau);

        printf("%d - %d\n", c0.cols(), c1.cols());

        // G2.col(i) = c0.cwiseProduct(c1) * (VectorXf::Ones(c0.cols()) * 1.0/c0.cols());
        // IP.col(i) = c0 * (VectorXf::Ones(c0.cols()) * 1.0/c0.cols());
        // IF.col(i) = c1 * (VectorXf::Ones(c1.cols()) * 1.0/c1.cols());

        i++;
        ll = level;
    }   

}

void Corr::multiTau2(data_structure::SparseData* data, float* G2s, float* IPs, float* IFs)
{
    Configuration* conf = Configuration::instance();
    int w = conf->getFrameWidth();
    int h = conf->getFrameHeight();
    int frames = conf->getFrameTodoCount();
    int pixels = w * h;

    int maxLevel = calculateLevelMax(frames, 4);

    vector<int> validPixels = data->ValidPixels();

    vector<tuple<int,int> > delays_per_level = delaysPerLevel(frames, 4, maxLevel);

    int pix = 355517;

    #pragma omp parallel for default(none) schedule(dynamic) shared(validPixels, delays_per_level, frames, pixels, G2s, IFs, IPs, data)
    for (int i = 0; i < validPixels.size(); i++)
    {
        data_structure::Row *row = data->Pixel(validPixels.at(i));

        int ll = 0;

        int lastframe = frames;
        int lastIndex = row->indxPtr.size();
        int tauIndex = 0;
        int g2Index = 0;
        
        // if (validPixels.at(i) == pix) {
        //     printf("Last Index %d\n", lastIndex);
        // }

        for (vector<tuple<int, int> >::iterator it = delays_per_level.begin() ;
                it != delays_per_level.end(); ++it)
        {

            // View of pixel data at two different tau values. 
            tuple<int, int> tau_level = *it;
            int level = get<0>(tau_level);
            int tau = get<1>(tau_level);
            // int lastIndex;

            // if (level == 0)
            //     lastIndex = row->indxPtr.size();

            if (ll != level)
            {   
                // if (level == 7 and validPixels.at(i) == pix)
                // {
                //     printf ("...");
                // }
                if (lastframe % 2)
                    lastframe -= 1;

                lastframe = lastframe / 2;

                int index = 0;
                int cnt = 0;
                
                int i0, i1;
                i0 = row->indxPtr[cnt] / 2.0;

                if (lastIndex > 0 && i0 < lastframe) {
                    row->indxPtr[index] = i0;
                    cnt = 1;
                }

                while (cnt < lastIndex) 
                {
                    i1 = row->indxPtr[cnt] / 2.0;
                    if (i1 >= lastframe) break;

                    if (i1 == i0) 
                    {
                        row->valPtr[index] += row->valPtr[cnt]; 
                    }
                    else 
                    {
                        row->indxPtr[++index] = i1;
                        row->valPtr[index] = row->valPtr[cnt];
                    }

                    i0 = i1;
                    cnt++;
                }

                lastIndex = index+1;

                for (int i = 0 ; i < lastIndex; i++) 
                    row->valPtr[i] /= 2.0f;



                // if (row->indxPtr.size() == 1) {
                //     row->indxPtr[index] = i0;
                //     row->valPtr[index] = row->valPtr[index]/2.0f;
                // }

                // int inc;
                // while (i1 < lastframe && cnt < lastIndex-1)
                // {   
                //     i0 = row->indxPtr[cnt] / 2;
                //     i1 = row->indxPtr[cnt+1] / 2;
                    
                //     inc = 0;
                //     if (i0 == i1)
                //     {
                //         row->indxPtr[index] = i0;
                //         row->valPtr[index] = (row->valPtr[cnt] + row->valPtr[cnt+1]) / 2.0f;
                //         inc = 2;
                //     }
                //     else
                //     {
                //         row->indxPtr[index] = i0;
                //         row->valPtr[index] = row->valPtr[cnt] / 2.0f;
                //         inc = 1;
                //     }

                //     cnt += inc;
                //     index++; 

                // }

                // // In case,we are left with one off element
                // if (i0 != i1 && i1 < lastframe && cnt < lastIndex)
                // {
                //     row->indxPtr[index] = i1;
                //     row->valPtr[index] = row->valPtr[cnt] / 2.0f;
                //     index++;
                // }

                // lastIndex = index;
                

                // if (validPixels.at(i) == pix) {
                //     printf("level = %d, lastframe = %d, lastIndex = %d\n", level, lastframe, lastIndex);
                //     for (int ii = 0; ii < lastIndex; ii++) {
                //         printf("%d -> %f\n", row->indxPtr[ii], row->valPtr[ii]);
                //     }
                //     printf("\n");
                // }
            }

            if (level > 0)
                tau = tau / pow(2, level);

            g2Index = tauIndex * pixels + validPixels[i]; // * delays_per_level.size() + tauIndex;
            //G2s[g2Index]  = 0.0;
            //IFs[g2Index]  = 0.0;
            //IPs[g2Index]  = 0.0;

            // if (validPixels.at(i) == pix && lastframe == 7) {
            //     printf("last index %d\n", lastIndex);
            // }

            for (int r = 0; r < lastIndex; r++)
            {
                int src = row->indxPtr.at(r);
                int dst = src;
                
                if (src < (lastframe-tau)) {
                    IPs[g2Index] += row->valPtr.at(r);
                    int limit = min(lastIndex, src+tau+1);
                    
                    for (int j = r+1; j < limit; j++)
                    {
                        dst = row->indxPtr.at(j);
                        if (dst == (src+tau)) {
                            G2s[g2Index] += row->valPtr.at(r) * row->valPtr.at(j);
                        }
                    }
                }

                if (src >= tau)
                    IFs[g2Index] += row->valPtr.at(r);

            }

            if ( (lastframe - tau) > 0) {
                G2s[g2Index] /= (lastframe-tau);
                IPs[g2Index] /= (lastframe-tau);
                IFs[g2Index] /= (lastframe-tau);
            }
            
            // if (validPixels.at(i) >= 94045 and validPixels.at(i) <= 94048) {
            //     printf("%d, g2Index=%d, tau=%d, IP=%f\n", validPixels[i], g2Index, tauIndex, IPs[g2Index]);
            // }

            ll = level;
            tauIndex++;
        }
    }
}

void Corr::twoTimesVec(Eigen::Ref<Eigen::MatrixXf> pixelData)
{

    Configuration* conf = Configuration::instance();
    int frames = conf->getFrameTodoCount();

    MatrixXf g2(pixelData.rows(), int(0.5 * frames * (frames - 1)) );

    int index = 0;
    for (int i = 0; i < frames; i++)
    {
        for (int j = i+1; j < frames; j++)
        {
            // g2.col(index) = pixelData.col(i) * pixelData.col(j);
            index++;
        }
    }

    // For self dot product of two times, create an upper triangular view of the original data
    // MatrixXf pixelData2 = MatrixXf(pixelData.triangularView<Eigen::Uppper>());
    // MatrixXf prod = pixelData.array() * pixelData2.array();

    // Sum products based on the q-bin pixels

    // map<int, map<int, vector<int>> > qbins = conf->getBinMaps();


}

//TODO: Refactor this function and possibly break into sub function for the unit-tests. 
void Corr::normalizeG2s(Eigen::Ref<Eigen::MatrixXf> G2,
                   Eigen::Ref<Eigen::MatrixXf> IP, 
                   Eigen::Ref<Eigen::MatrixXf> IF)
{

    Configuration* conf = Configuration::instance();
    int w = conf->getFrameWidth();
    int h = conf->getFrameHeight();

    std::map<int, std::map<int, vector<int>> > qbins = conf->getBinMaps();

    int totalStaticPartns = conf->getTotalStaticPartitions();
    int totalDynamicPartns = conf->getTotalDynamicPartitions();

    // Final result out of this function. 
    MatrixXf g2(qbins.size(), G2.cols());
    MatrixXf stdError(qbins.size(), G2.cols());

    g2.setZero(g2.rows(), g2.cols());
    stdError.setZero(stdError.rows(), stdError.cols());

    // A matrix of the form total_static partitions over tau values.
    // Each row is a static partition and column is different tau value.  
    MatrixXf g2Sums(totalStaticPartns, G2.cols());
    MatrixXf ipSums(totalStaticPartns, G2.cols());
    MatrixXf ifSums(totalStaticPartns, G2.cols());

    g2Sums.setZero(g2Sums.rows(), g2Sums.cols());
    ipSums.setZero(ipSums.rows(), ipSums.cols());
    ifSums.setZero(ifSums.rows(), ifSums.cols());

    // Count of pixels for each static partition
    VectorXi partitionPixelCounts(totalStaticPartns);
    partitionPixelCounts.setZero(totalStaticPartns);

    float sum = 0.0;
    float sum2 = 0.0;

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
                
                g2Sums.row(sbin - 1) += G2.row(p);
                ipSums.row(sbin - 1) += IP.row(p);
                ifSums.row(sbin - 1) += IF.row(p);

                if ( (sbin-1) == 0 ) {
                    sum += IP(p, 0);
                    sum2 += IP(p, 1);
                }

                partitionPixelCounts(sbin-1) += 1;
            }
        }
    }

    // Compute averag of each static parition
    VectorXi sbinCounts(totalDynamicPartns);
    sbinCounts.setZero(totalDynamicPartns);

    for (map<int, map<int, vector<int>> >::const_iterator it = qbins.begin(); 
            it != qbins.end(); it++) {

        int q = it->first;
        map<int, vector<int> > values =  it->second;

        for (map<int, vector<int>>::const_iterator it2 =  values.begin(); it2 != values.end(); it2++) {
            int sbin = it2->first;

            g2Sums.row(sbin-1).array() /= partitionPixelCounts(sbin-1);
            ipSums.row(sbin-1).array() /= partitionPixelCounts(sbin-1);
            ifSums.row(sbin-1).array() /= partitionPixelCounts(sbin-1);

            // Normalize per bin static partition.

            g2Sums.row(sbin-1) = g2Sums.row(sbin-1).array() / (ipSums.row(sbin-1).array() * ifSums.row(sbin-1).array());
        }
    }

    // // Compute the mean of normalized g2 values per dynamic bin. 
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
        g2.row(q - 1).array() = g2.row(q - 1).array() / (float)count;
    }

    // Compute the standard error.
    for (map<int, map<int, vector<int>> >::const_iterator it = qbins.begin(); 
            it != qbins.end(); it++) {

        int q = it->first;
        map<int, vector<int> > values =  it->second;

        VectorXf avgG2(G2.cols());
        avgG2.setZero();
        VectorXf tmpAvg(G2.cols());
        tmpAvg.setZero();
        VectorXf stdG2(G2.cols());
        stdG2.setZero();
        VectorXf samples(G2.cols());
        samples.setZero();
        VectorXf incrs(G2.cols());
        incrs.setOnes();

        for (map<int, vector<int>>::const_iterator it2 =  values.begin(); it2 != values.end(); it2++) {
            int sbin = it2->first;

            vector<int> pixels = it2->second;
            for(vector<int>::const_iterator pind = pixels.begin(); pind != pixels.end(); pind++) {
                int p = *pind;
                tmpAvg = avgG2;

                VectorXf normalizedG2 =  G2.row(p).array() / (IP.row(p).array() * IF.row(p).array());
                samples += incrs;
                normalizedG2 = Eigen::isinf(normalizedG2.array()).select(0, normalizedG2);
                samples += Eigen::isinf(normalizedG2.array()).select(0, incrs);

                VectorXf tmp0 = (normalizedG2 - tmpAvg);
                VectorXf tmp00 = tmp0.array() / samples.array();
                avgG2 += tmp00;
                VectorXf tmp1 = (normalizedG2 - tmpAvg);
                VectorXf tmp2 = (normalizedG2 - avgG2);
                VectorXf tmp3 = tmp1.array() * tmp2.array();

                stdG2 += tmp3;
            }
        }

        VectorXf stdNorm = stdG2.array() / samples.array();
        samples = 1.0f / samples.array();
        // std::cout<<samples.array().sqrt()<<std::endl;
        stdError.row(q - 1).array() = samples.array().sqrt() * stdNorm.array().sqrt();
    }

    H5Result::write2DData(conf->getFilename(), "exchange", "norm-0-g2", g2);
    H5Result::write2DData(conf->getFilename(), "exchange", "norm-0-stderr", stdError);

}

double* Corr::computeG2Levels(const Eigen::MatrixXf &pixelData, 
                                int pixel,
                                int frameCount, 
                                int tau, 
                                int level) 
{

    int tauIncrement = (level == 0) ? 1 : (int) pow(2.0, level);
    double numerator = 0.0;
    double sumPast = 0.0;
    double sumFuture = 0.0;

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

} //namespace xpcs
