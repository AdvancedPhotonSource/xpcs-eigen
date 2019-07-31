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

#include "sparse_filter_burst.h"

#include <math.h>
#include <set>

#include <stdio.h>
#include <iostream>

#include "xpcs/configuration.h"
#include "xpcs/io/reader.h"
#include "xpcs/data_structure/sparse_data.h"

namespace xpcs {
namespace filter {


SparseFilterBurst::SparseFilterBurst() {
  Configuration *conf = Configuration::instance();

  pixel_mask_ = conf->getPixelMask();
  sbin_mask_ = conf->getSbinMask();
  flatfield_ = conf->getFlatField();
  frames_todo_ = conf->getFrameTodoCount();
  frame_width_ = conf->getFrameWidth();
  frame_height_ = conf->getFrameHeight();
  stride_size_ = conf->FrameStride();
  average_size_ = conf->FrameAverage();
  normalizedByFramesum_ = conf->IsNormalizedByFramesum();
  burst_size_ = conf->NumberOfBursts();

  pixels_value_ = new float[frame_width_ * frame_height_];
  sparse_map_ = new short[frame_width_ * frame_height_];
  real_frames_todo_ = conf->getRealFrameTodoCount();
  
  frame_index_ = 0;
  
  data_ = new xpcs::data_structure::SparseData(frame_width_ * frame_height_);
}

SparseFilterBurst::~SparseFilterBurst() {

}

void SparseFilterBurst::Apply(xpcs::io::ImmBlock* blk) {
  int **indx = blk->index;
  float **val = blk->value;
  int frames = blk->frames;
  int pix_cnt = frame_width_ * frame_height_;
  std::vector<int> ppf = blk->pixels_per_frame;

  for (int j = 0; j < pix_cnt; j++) {
    pixels_value_[j] = 0.0f;
    sparse_map_[j] = 0;
  }

  if (burst_data_) {
    //TODO
  }

  burst_data_ = new xpcs::data_structure::SparseData(frame_width_ * frame_height_);
  assert(frames == burst_size_);

  // Keep track of pixels that were part of any of the frame. 
  for (int i = 0; i < burst_size_; i++) {
    int pixels = ppf[i];
    int *index = indx[i];
    float *value = val[i];

    for (int j = 0; j < pixels; j++) {

      if (pixel_mask_[index[j]] != 0) {
        int pix = index[j];
        float v = value[j] * flatfield_[pix];

        xpcs::data_structure::Row *row = burst_data_->Pixel(pix);
        row->indxPtr.push_back(pix);
        row->valPtr.push_back(v);

        pixels_value_[pix] += v;
        sparse_map_[pix] = 1; 
      }
    }

  }

  int ind = 0;
  for (int i = 0 ; i < pix_cnt; i++) {
    if (sparse_map_[i] == 0) continue;

    int pix = i;

    float v = pixels_value_[pix];

    xpcs::data_structure::Row *row = data_->Pixel(pix);
    row->indxPtr.push_back(frame_index_);
    row->valPtr.push_back(v);

  }
  frame_index_++;

}

xpcs::data_structure::SparseData* SparseFilterBurst::Data() {
  return data_;
}

xpcs::data_structure::SparseData* SparseFilterBurst::BurstData() {
  return burst_data_;
}

} // namespace io
} // namespace xpcs
