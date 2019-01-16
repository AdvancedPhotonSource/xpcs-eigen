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
#include "ufxc_reader.h"

#include <stdio.h>
#include <iostream>
#include <vector>
#include <iterator>

#include "xpcs/configuration.h"

namespace xpcs {
namespace io {

UfxcReader::UfxcReader(const std::string& filename) {
    file_ = fopen(filename.c_str(), "rb");
    if (file_ == NULL) return ; //TODO handle error

    std::vector<uint> data;

    uint buffer[4096];
    size_t read = fread(&buffer, sizeof(uint), 4096, file_);
    while (read) {
        for (int i = 0; i < read; i++) {
            data.push_back(buffer[i]);
        }
        read = fread(&buffer, sizeof(uint), 4096, file_);
    }

    auto it = data.begin();
    uint value = *it >> 21;
    uint f0 = value;
    int counter = 0;
    int idx = 1;
    data_frames_[value - f0] = std::vector<uint>();
    data_frames_[value-f0].push_back(*it);

    ++it;
    int bf = 0;
    int ff = 0; // frame number
    int tmp_count = 0;
    for(; it != data.end(); ++it) {
        int diff = (*it >> 21) - value;
        if (diff < -2000) {
            bf += 2048;
        } else if (diff > 2000) {
            bf -= 2048;
        } 
        ff = (*it >> 21) + bf - f0;
        if (data_frames_.find(ff) == data_frames_.end()) { 
            data_frames_[ff] = std::vector<uint>();
	}
        data_frames_[ff].push_back(*it);
        value = *it >> 21;
    }
   
    last_frame_index = 0;

}

UfxcReader::~UfxcReader() {
}

ImmBlock* UfxcReader::NextFrames(int count) {
    int **index = new int*[count];
    float **value = new float*[count];
    double *clock = new double[count];
    double *ticks = new double[count];

    std::vector<int> ppf;


    int done = 0, pxs = 0;

    while (done < count) {
        clock[done] = last_frame_index;
        ticks[done] = last_frame_index;

        if (data_frames_.find(last_frame_index) == data_frames_.end()) {
            index[done] = new int[0];
            value[done] = new float[0];

            ppf.push_back(0);
            last_frame_index++;
            done++;
            continue;
        }

        std::vector<uint> frame = data_frames_[last_frame_index];
        index[done] = new int[frame.size()];
        value[done] = new float[frame.size()];
        ppf.push_back(frame.size());

        int idx = 0;
        for (auto& it : frame) {
            int pix = it & 0x7fff;
            float val = (it >> 15) & 0x3;
            index[done][idx] = pix;
            value[done][idx] = val;
	    idx++;
        }
        done++;
        last_frame_index++;
    }
   
    ImmBlock *ret = new ImmBlock;
    ret->index = index;
    ret->value = value;
    ret->frames = count;
    ret->pixels_per_frame = ppf;
    ret->clock = clock;
    ret->ticks = ticks;
    ret->id = 1;

    return ret;
}

void UfxcReader::SkipFrames(int count) {
    int done = 0;
}

void UfxcReader::Reset() {
    rewind(file_);
}

bool UfxcReader::compression() { return true; }

} // namespace io
} // namespace xpcs
