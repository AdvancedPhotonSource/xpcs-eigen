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
#include "rigaku.h"

#include <stdio.h>
#include <iostream>
#include <vector>
#include <iterator>

#include "xpcs/configuration.h"
#include "xpcs/data_structure/sparse_data.h"
#include "xpcs/filter/filter.h"

#include "../benchmark.h"

namespace xpcs {
namespace io {

Rigaku::Rigaku(const std::string& filename, xpcs::filter::Filter *filter) {
  xpcs::Benchmark benchmark("Reading rigaku file");

  Configuration* conf = Configuration::instance();
  int frame_width = conf->getFrameWidth();
  int frame_height = conf->getFrameHeight();
  short *pixel_mask = conf->getPixelMask();
  double *flatfield = conf->getFlatField();
  int static_window = conf->getStaticWindowSize();
  int *sbin_mask = conf->getSbinMask();
  int total_static_partns = conf->getTotalStaticPartitions();
  int swindow = conf->getStaticWindowSize();
  int frames = conf->getFrameTodoCount();
  int partitions = (int) ceil((double)frames/swindow);
  float *partial_partitions_mean = new float[total_static_partns * partitions];
  float *partitions_mean = new float[total_static_partns];

  int framesize = frame_width * frame_height;
  
  float *pixels_sum = new float[framesize];
  float *frames_sum =  new float[2 * conf->getFrameTodoCount()];
    
  xpcs::data_structure::SparseData *data = new xpcs::data_structure::SparseData(framesize);

  file_ = fopen(filename.c_str(), "rb");
  if (file_ == NULL) return ; //TODO handle error

  long buffer_size = 4096 * 10;
  long long buffer[buffer_size];
  size_t read = fread(&buffer, sizeof(long long), buffer_size, file_);

  int *ppf = new int[frames];

  for (int i = 0; i < frames; i++) {
      ppf[i] = 0;
  }

  for (int i = 0; i < total_static_partns; i++) 
    partitions_mean[i] = 0.0;

  for (int i = 0; i < total_static_partns * partitions; i++)
    partial_partitions_mean[i] = 0.0;

  for (int i = 0; i < framesize; i++)
    pixels_sum[i] = 0.0f;

  uint previous_frame = 0;    
  uint frame = 0;
  int pixcount = 0;
  int totalpix = 0;

  int sbin = 0;
  int partition_no = 0;

  float fsum = 0.0;

  while (read) 
  {
    for (int i = 0; i < read; i++) 
    {
        frame = buffer[i] >> 40;    
        if (frame != previous_frame) 
        {
            ppf[previous_frame] = pixcount;
            pixcount = 0;
            previous_frame = frame;

            frames_sum[frame] = frame + 1.0;
            frames_sum[frame + frames] = fsum / (float)framesize;

            if (frame > 0 && (frame % static_window) == 0) 
            {
              partition_no++;
            }

            fsum = 0.0;
        }
        
        uint pix = (buffer[i] >> 16) & 0xFFFFF;
          
        if (pixel_mask[pix] == 0) 
          continue;

          
        float val = buffer[i] & 0x7FF;
        val *= flatfield[pix];

        xpcs::data_structure::Row *row = data->Pixel(pix);

        row->indxPtr.push_back(frame);
        row->valPtr.push_back(val);

        sbin = sbin_mask[pix] - 1;
    
        partitions_mean[sbin] += val;

        pixels_sum[pix] += val;

        fsum += val;

        pixcount++;
        totalpix++;
    }
    
    read = fread(&buffer, sizeof(long long), buffer_size, file_);
  }

  filter->FramesSum(frames_sum);
  filter->PixelsSum(pixels_sum);
  filter->PartitionsMean(partitions_mean);
  filter->PartialPartitionsMean(partial_partitions_mean);
  filter->Data(data);

}

Rigaku::~Rigaku() 
{
}

void Rigaku::Process() 
{
  
}

ImmBlock* Rigaku::NextFrames(int count) {
    // int **index = new int*[count];
    // float **value = new float*[count];
    // double *clock = new double[count];
    // double *ticks = new double[count];

    // std::vector<int> ppf;
    // int done = 0, pxs = 0;

    // while (done < count) {
    //     clock[done] = last_frame_index_;
    //     ticks[done] = last_frame_index_;

    //     if (ppf_[last_frame_index_] == 0) {
    //         index[done] = new int[0];
    //         value[done] = new float[0];

    //         ppf.push_back(0);
    //         last_frame_index_++;
    //         done++;
    //         continue;
    //     }
        
    //     int pix_cnt = ppf_[last_frame_index_];
    //     index[done] = new int[pix_cnt];
    //     value[done] = new float[pix_cnt];
    //     ppf.push_back(pix_cnt);

    //     uint idx = 0;
    //     for (int i = last_pixel_index_; i < (last_pixel_index_ + pix_cnt); i++) {
    //         long long it = data_[i];

    //         uint pix = (it >> 16) & 0xFFFFF;
            
    //         int row = pix % frame_height_;
    //         int col = pix / frame_height_;

    //         float val = it & 0x7FF;

    //         index[done][idx] = pix;
    //         value[done][idx] = val * flatfield[pix];;
    //         idx++;
    //     }
    //     done++;
    //     last_frame_index_++;
    //     last_pixel_index_ += ppf[last_frame_index_];
    // }
   
    ImmBlock *ret = new ImmBlock;
    ret->index = NULL;
    ret->value = NULL;
    ret->frames = 0;
    ret->pixels_per_frame = std::vector<int>();
    ret->clock = NULL;
    ret->ticks = NULL;
    ret->id = -1;

    return ret;
}

void Rigaku::SkipFrames(int count) {
    int done = 0;
}

void Rigaku::Reset() {
    rewind(file_);
}

bool Rigaku::compression() { return true; }

} // namespace io
} // namespace xpcs
