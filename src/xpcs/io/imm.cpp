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
#include "imm.h"
// #include "reader.h"

#include <stdio.h>
#include <iostream>

#include "xpcs/configuration.h"
#include "xpcs/data_structure/dark_image.h"

namespace xpcs {
namespace io {

Imm::Imm(const std::string& filename)  {
    file_ = fopen(filename.c_str(), "rb");

    if (file_ == NULL) return ; //TODO handle error

    header_ = new Header();
    fread(header_, 1024, 1, file_);
    compression_ = header_->compression != 0 ? true : false;
    rewind(file_);
}

ImmBlock* Imm::NextFrames(int count) {
    int **index = new int*[count];
    float **value = new float*[count];
    double *clock = new double[count];
    double *ticks = new double[count];

    std::vector<int> ppf;

    int done = 0, pxs = 0;
    while (done < count) {
        // printf("inloop ftell(file_) %ld\n", ftell(file_));
        fread(header_, 1024, 1, file_);
        pxs = header_->dlen;
        // printf("Buffer # = %ld, pxs = %d\n", header_->buffer_number, pxs);

        index[done] = new int[pxs];
        value[done] = new float[pxs];
        short *tmpmem = new short[pxs];
        
        if (compression_) {
            fread(index[done], pxs * 4, 1, file_);
        } 
        
        // else {
        //     // for (int i = 0 ; i < pxs; i++)
        //     //     index[done][i] = i;
        // }

        fread(tmpmem, pxs * 2, 1, file_);
        std::copy(tmpmem, tmpmem + pxs, value[done]);
        delete [] tmpmem;
        ppf.push_back(pxs);

        clock[done] = header_->elapsed;
        ticks[done] = header_->corecotick;
        done++;
    }

    struct ImmBlock *ret = new ImmBlock;
    ret->index = index;
    ret->value = value;
    ret->frames = count;
    ret->pixels_per_frame = ppf;
    ret->clock = clock;
    ret->ticks = ticks;
    ret->id = 0;

    return ret;
}

void Imm::SkipFrames(int count) {
    int done = 0;
    int image_bytes = compression_ ? 6 : 2;

    while (done < count) {
        fread(header_, 1024, 1, file_);
        fseek(file_, header_->dlen * image_bytes, SEEK_CUR);
        done++;
    }
}

void Imm::Reset() {
    rewind(file_);
}

bool Imm::compression() { return compression_; }

} // namespace io
} // namespace xpcs
