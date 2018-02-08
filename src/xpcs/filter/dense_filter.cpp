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

#include "dense_filter.h"

#include <math.h>

#include <stdio.h>
#include <iostream>

#include "xpcs/configuration.h"
#include "xpcs/io/imm_reader.h"
#include "xpcs/data_structure/sparse_data.h"

namespace xpcs {
namespace filter {


DenseFilter::DenseFilter() {
  Configuration *conf = Configuration::instance();

  pixel_mask_ = conf->getPixelMask();
  sbin_mask_ = conf->getSbinMask();
  flatfield_ = conf->getFlatField();
  eff_ = conf->getDetEfficiency();
  det_adhu_ = conf->getDetAdhuPhot();
  preset_ =  conf->getDetPreset();
  frames_todo_ = conf->getFrameTodoCount();
  frame_width_ = conf->getFrameWidth();
  frame_height_ = conf->getFrameHeight();
  normFactor_ = conf->getNormFactor();
  static_window_ = conf->getStaticWindowSize();
  total_static_partns_ = conf->getTotalStaticPartitions();
  swindow_ = conf->getStaticWindowSize();
  int partitions = (int) ceil((double)frames_todo_/swindow_);
  partial_partitions_mean_ = new float[total_static_partns_ * partitions];
  partitions_mean_ = new float[total_static_partns_];
  pixels_sum_ = new float[frame_width_ * frame_height_];
  frames_sum_ =  new float[2 * conf->getFrameTodoCount()];

  frame_index_ = 0;
  partition_no_ = 0;

  for (int i = 0; i < total_static_partns_; i++) 
    partitions_mean_[i] = 0.0;

  for (int i = 0; i < total_static_partns_ * partitions; i++)
    partial_partitions_mean_[i] = 0.0;

  for (int i = 0; i < (frame_width_ * frame_height_); i++)
    pixels_sum_[i] = 0.0f;

  data_ = new xpcs::data_structure::SparseData(frame_width_ * frame_height_);

}

DenseFilter::~DenseFilter() {

}

void DenseFilter::Apply(struct xpcs::io::ImmBlock* blk) {
  int **indx = blk->index;
  short **val = blk->value;
  int frames = blk->frames;
  std::vector<int> ppf = blk->pixels_per_frame;

  for (int i = 0; i < frames; i++) {
    int pixels = ppf[i];
    int *index = indx[i];
    short *value = val[i];
    float f_sum = 0.0f;
    int pix_cnt = 0;
    int sbin = 0;

    if (frame_index_ > 0 && (frame_index_ % static_window_) == 0) {
      partition_no_++;
    }

    for (int j = 0; j < pixels; j++) {

      if (pixel_mask_[index[j]] != 0) {
        int pix = index[j];
        float v = (float)value[j] * flatfield_[pix];

        pixels_sum_[pix] += v;
        f_sum += v;

        sbin = sbin_mask_[pix] - 1;

        xpcs::data_structure::Row *row = data_->Pixel(pix);
        row->indxPtr.push_back(frame_index_);
        row->valPtr.push_back(v);

        partitions_mean_[sbin] += v;
        partial_partitions_mean_[partition_no_ * total_static_partns_ + sbin ] += v;
        pix_cnt++;
      }
    }

    frames_sum_[frame_index_] = frame_index_ + 1.0;
    frames_sum_[frame_index_ + frames_todo_] = f_sum / pix_cnt;
    frame_index_++;
  }

}

float* DenseFilter::PixelsSum() {
  return pixels_sum_;
}

float* DenseFilter::FramesSum() {
  return frames_sum_;
}

float* DenseFilter::PartitionsMean() {
  return partitions_mean_;
}

float* DenseFilter::PartialPartitionsMean() {
  return partial_partitions_mean_;
}


xpcs::data_structure::SparseData* DenseFilter::Data() {
  return data_;
}

} // namespace io
} // namespace xpcs