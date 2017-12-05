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

int main(int argc, char** argv)
{
    srand (time(NULL));

    int rows = 10;
    int cols = 5;
    int rowno = atoi(argv[1]);
    printf("Printing row # %d\n", rowno);

    auto console = spd::stdout_color_mt("console");

    SparseMatrix<float, RowMajor> mat(rows, cols);
    std::vector<T> tripletList;
    {
        Benchmark b("Generating 25% sparse matrix");
        for (int i = 0; i < rows; i++)
        {
            for (int j = 0; j < cols; j++)
            {
                if ((rand() % 100) > 50)
                {
                    tripletList.push_back(T(i, j, 1+random()%30));
                }
            }
        }

        mat.setFromTriplets(tripletList.begin(), tripletList.end());
        mat.makeCompressed();
    }


    cout<<mat<<endl;
    {
        Benchmark b("Averaging out matrix");
        float* vals = mat.valuePtr();
        int *ind = mat.innerIndexPtr();
        int *outer = mat.outerIndexPtr();
        int nonzeros = outer[rowno+1] - outer[rowno];

        ind += outer[rowno];
        vals += outer[rowno];

        printf("R%d :", rowno);

        for (int i = 0; i < nonzeros; i++)
        {
            printf(" %f ", *vals);   
            vals++;
            ind++;
        }

        printf("\n");

        
        // #pragma omp parallel for num_threads(5) default(none) shared(mat, rows, cols)
        // for (int r = 0; r < rows; r++)
        // {

        //     // SparseMatrix<float, RowMajor> c0 = mat.block(r, 0, 1, cols);
        //     // SparseMatrix<float, RowMajor> c1 = mat.block(r, 1, 1, cols-1);
        //     // printf("Thread # %d and sum %f, %f \n", omp_get_thread_num(), c0.sum(), c1.sum());
        //     mat.row(r) = 0.5 * c0;
        // }

    }    
    // cout<<mat<<endl;
}
