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
#include <omp.h>

#include "hdf5.h"
#include "imm.h"
#include "corr.h"
#include "configuration.h"
#include "h5result.h"
#include "funcs.h"
#include "benchmark.h"

#include "Eigen/Dense"
#include "Eigen/SparseCore"

#include "gflags/gflags.h"
#include "spdlog/spdlog.h"

using namespace Eigen;
using namespace std;

typedef Eigen::Triplet<double> T;

void correlaton1()
{
       // for (int tau = 1; tau <= 2; tau++)
       //  {
       //      const SparseMatrix<float, RowMajor> c0 = mat.block(0, 0, rows, cols-tau);
       //      const SparseMatrix<float, RowMajor> c1 = mat.block(0, tau, rows, cols-tau);
       //      result.col(tau-1) = c0.cwiseProduct(c1) * (VectorXf::Ones(c0.cols()) * 1.0/c0.cols());
       //  }
}

int main(int argc, char** argv)
{
    srand (time(NULL));

    int rows = 10;
    int cols = 5;
    // int rowno = atoi(argv[1]);
    // printf("Printing row # %d\n", rowno);

    auto console = spd::stdout_color_mt("console");

    SparseMatrix<float, RowMajor> mat(rows, cols);
    std::vector<T> tripletList;
    {
        Benchmark b("Generating 25% sparse matrix");
        for (int i = 0; i < rows; i++)
        {
            for (int j = 0; j < cols; j++)
            {
                int r = rand() % 100;
                if ( r > 75 )
                {
                    tripletList.push_back(T(i, j, r));
                }
            }
        }

        mat.setFromTriplets(tripletList.begin(), tripletList.end());
        mat.makeCompressed();
    }


    cout<<mat<<endl;
    // MatrixXf result(rows, 2);
    {
        Benchmark b("Averaging out matrix");
        float* vals = mat.valuePtr();
        int *ind = mat.innerIndexPtr();
        int *outer = mat.outerIndexPtr();

        #pragma omp parallel for default(none) shared(mat, ind, vals, outer, rows, cols)
        for (int r = 0; r < rows; r++)
        {
            int nonzeros = outer[r+1] - outer[r];

            // printf("Number of nonzeros at row %d are %d\n", r, nonzeros);

            if (!nonzeros) continue;

            int *iptr = ind + outer[r];
            float *vptr = vals + outer[r];

            int *ciptr = iptr;
            float *cvptr = vptr;

            // Limit is number of frames we need to have for smoothing. 
            //   This gets halved after each iteration (level). 
            int limit = cols/2;
            if (limit % 2)
                limit -= 1;

            
            int index = 0;
            int cnt = 0;
            int i0, i1;
            i0 = 0;

            if (nonzeros == 1) {
                iptr[index] = i0;
                vptr[index] = (0.0 + *vptr)/2;
                continue;
            }

            int inc;
            while (i0 < limit && cnt < (nonzeros-1))
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
            if (i0 != i1 && i1 < limit)
            {
                ciptr[index] = i1;
                cvptr[index] = (0.0f + *vptr) / 2.0f;
            } 

        }
             
    }    
    cout<<mat<<endl;

    SparseMatrix<float, RowMajor> c0 = mat.block(0, 0, rows, cols-3);
    cout<<c0<<endl;
       //      const SparseMatrix<float, RowMajor> c1 = mat.block(0, tau, rows, cols-tau);

}
