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
#ifndef CORR_H
#define CORR_H

#include <vector>
#include <tuple>

#include "Eigen/Dense"
#include "Eigen/SparseCore"

#include "imm.h"
#include "sparsedata.h"


class Corr {

public:
    static int calculateDelayCount(int dpl, int level);

    static int calculateLevelMax(int frameCount, int dpl);

    static std::vector< std::tuple<int, int> > delaysPerLevel(int frameCount, int dpl, int maxDelay);


    /**
     * Non vectorized version of G2 calcuation, here only for reference. 
     * I will delete it from the source soon. 
     */
    static void multiTau(const Eigen::MatrixXf &pixelData, int pix);

    /**
     *  Compute G2,IP, and IF from a Dense matrix. 
     */
    static void multiTauVec(Eigen::Ref<Eigen::MatrixXf> pixelData, 
                            Eigen::Ref<Eigen::MatrixXf> G2, 
                            Eigen::Ref<Eigen::MatrixXf> IP, 
                            Eigen::Ref<Eigen::MatrixXf> IF);

    /**
     *  Compute G2,IP, and IF from a Sparse Matrix. 
     */
    static void multiTauVec(SparseRMatF& pixelData,
                            Eigen::Ref<Eigen::MatrixXf> G2, 
                            Eigen::Ref<Eigen::MatrixXf> IP, 
                            Eigen::Ref<Eigen::MatrixXf> IF);

    static void multiTau2(SparseData *data, float* G2, float* IP, float* IF);

    static void twoTimesVec(Eigen::Ref<Eigen::MatrixXf> pixelData);

    static void normalizeG2s(Eigen::Ref<Eigen::MatrixXf> g2,
                      Eigen::Ref<Eigen::MatrixXf> IP, 
                      Eigen::Ref<Eigen::MatrixXf> IF);

    static double* computeG2Levels(const Eigen::MatrixXf &pixelData, 
                                int pixel,
                                int frameCount, 
                                int tau, 
                                int level);

};

#endif