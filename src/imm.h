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


#ifndef IMM_H
#define IMM_H

#include "immHeader.h"
#include "sparsedata.h"
#include <stdio.h>
#include <stdlib.h>
#include <string>

#include "Eigen/Dense"
#include "Eigen/SparseCore"
#include "spdlog/spdlog.h"

using namespace Eigen;

typedef Eigen::SparseMatrix<float> SparseMatF;
typedef Eigen::SparseMatrix<float, RowMajor> SparseRMatF;

class IMM
{
public:
   IMM(const char* filename, int frameFrom, int frameTo, int pixels_per_frame);
   
   ~IMM();

   Eigen::MatrixXf getPixelData();
   
   SparseRMatF getSparsePixelData();

   bool getIsSparse();

   float* getTimestampClock();
   float* getTimestampTick();

   float* getFrameSums();
   float* getPixelSums();

   float* getTotalPartitionMean();
   float* getPartialPartitionMean();

   SparseData* getSparseData();

private:

    long m_frameStartTodo;
    long m_frameEndTodo;
    long m_pixelsPerFrame;
    long m_frames;
    
    const char* m_filename;
    IMMHeader *m_ptrHeader;

    FILE *m_ptrFile;

    // Internal data pointer for storing pixels. 
    float *m_data;
    float *m_timestampClock;
    float *m_timestampTick;

    MatrixXf m_pixelData;
    SparseRMatF m_sparsePixelData;
    
    bool m_isSparse;

    // Initialize the file ptr and read-in file header. 
    void init();

    // Loads the sparse IMM file 
    void load_sparse();

    // Loads the sparse IMM to internanl structures. Unlinke the load_sparse method
    // it doesn't generate a matrix. 
    void load_sparse2();

    // Load non-sparse data. 
    void load_nonsprase();
    void load_nonsparse2();

    SparseData *m_sdata;

    std::shared_ptr<spdlog::logger> _logger;

    float* m_pixelSums;
    float* m_frameSums;
    float* m_partialPartitionMean;
    float* m_totalPartitionMean;
};

#endif