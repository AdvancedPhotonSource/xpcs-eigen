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
  uint frames_start_todo = conf->getFrameStartTodo();
  uint frames = conf->getFrameTodoCount();
  int partitions = (int) ceil((double)frames/swindow);
  float *partial_partitions_mean = new float[total_static_partns * partitions];
  float *partitions_mean = new float[total_static_partns];

  int stride_factor = conf->FrameStride();
  int average_factor = conf->FrameAverage();

  int framesize = frame_width * frame_height;
  float *pixels_sum = new float[framesize];
  float *frames_sum =  new float[2 * conf->getFrameTodoCount()];

  double *clock = new double[frames];
  double *ticks = new double[frames];
    
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
  int sbin = 0;
  int partition_no = 0;

  int read_in_count = stride_factor > 1 ? stride_factor : average_factor;
  if (stride_factor > 1 && average_factor > 1)
    read_in_count = stride_factor * average_factor;

  std::map<int, float> pixel_key_value;
  // std::vector<int> sparse_pixel_mask;
  // int sparse_pixel_mask[framesize];
  // float pixel_values[framesize];

  // for (int i = 0 ; i < framesize; i++) {
  //   sparse_pixel_mask[i] = 0;
  //   pixel_values[i] = 0.0;
  // }
  
  int real_frame_index = 0;
  int next_expected_frame = frames_start_todo + read_in_count;
  previous_frame = frames_start_todo + 1;
  bool not_frame_endtodo = true;

  while (read && not_frame_endtodo) 
  {
    for (int i = 0; i < read; i++) 
    {
        frame = (buffer[i] >> 40);

        if (frame <= frames_start_todo) {
          continue;
        }

        // We reached to end of frames required for this read.
        if (real_frame_index >= frames) {
          not_frame_endtodo = false;
          break;
        }
        
        if (stride_factor > 1 && (frame != 0 && (frame % stride_factor != 0)))
          continue;

        // Two set of conditions, 1) averaging is on and we reached the boundary of frames to be averaged
        // 2) averaging is off and we reached the end of frame. 
        if ( (average_factor > 1 && frame > next_expected_frame) ||
             (average_factor == 1 && frame != previous_frame) ) 
        {
          int pixcount = 0;
          float fsum = 0.0;

          for (auto it = pixel_key_value.begin(); it != pixel_key_value.end(); ++it)
          {
            xpcs::data_structure::Row *row = data->Pixel(it->first);
            float value = it->second / average_factor;

            row->indxPtr.push_back(real_frame_index);
            row->valPtr.push_back(value);

            sbin = sbin_mask[it->first] - 1;
            
            // printf("Total static partitions %d, sbin %d\n", total_static_partns, sbin);
            partitions_mean[sbin] += value;
            partial_partitions_mean[partition_no * total_static_partns + sbin ] += value;
            pixels_sum[it->first] += value;
            fsum += value;

            pixcount++;
          }

          ppf[real_frame_index] = pixcount;

          frames_sum[real_frame_index] = real_frame_index + 1.0;
          frames_sum[real_frame_index + frames] = fsum / (float)framesize;
          clock[real_frame_index] = real_frame_index;
          ticks[real_frame_index] = real_frame_index;
          
          // printf("%d %d %f %f\n", frame, real_frame_index, fsum, frames_sum[real_frame_index + frames]);


          if (real_frame_index > 0 && (real_frame_index % static_window) == 0) 
          {
            partition_no++;
          }

          previous_frame = frame;
          real_frame_index++;
          next_expected_frame += read_in_count;

          pixel_key_value.clear();
        }

        

        
        uint pix = (buffer[i] >> 16) & 0xFFFFF;
        
        pix = (pix % frame_height) * frame_width + (pix / frame_height);

        if (pixel_mask[pix] == 0) 
          continue;
          
        float val = buffer[i] & 0x7FF;
        val *= flatfield[pix];

        auto key_exists = pixel_key_value.find(pix);
        if (key_exists != pixel_key_value.end()) {
          pixel_key_value[pix] += val;
        } else {
          pixel_key_value[pix] = val;
        }

    }
    read = fread(&buffer, sizeof(long long), buffer_size, file_);

  }

  // Edge cases when we don't have even number of frames / average*stride_factor.
  // drain sparse_pixel_mask.  
  if (pixel_key_value.size() > 1)
  {
    int pixcount = 0;
    float fsum = 0.0;

    for (auto it = pixel_key_value.begin(); it != pixel_key_value.end(); ++it)
    {
      xpcs::data_structure::Row *row = data->Pixel(it->first);
      float value = it->second / average_factor;

      row->indxPtr.push_back(real_frame_index);
      row->valPtr.push_back(value);

      sbin = sbin_mask[it->first] - 1;
            
      // printf("Total static partitions %d, sbin %d\n", total_static_partns, sbin);
      partitions_mean[sbin] += value;
      partial_partitions_mean[partition_no * total_static_partns + sbin ] += value;
      pixels_sum[it->first] += value;
      fsum += value;
      pixcount++;
    }

    ppf[real_frame_index] = pixcount;

    frames_sum[real_frame_index] = real_frame_index + 1.0;
    frames_sum[real_frame_index + frames] = fsum / (float)framesize;
    clock[real_frame_index] = real_frame_index;
    ticks[real_frame_index] = real_frame_index;

  }

  filter->FramesSum(frames_sum);
  filter->PixelsSum(pixels_sum);
  filter->PartitionsMean(partitions_mean);
  filter->PartialPartitionsMean(partial_partitions_mean);
  filter->TimestampClock(clock);
  filter->TimestampTicks(ticks);
  filter->Data(data);
}


Rigaku::~Rigaku() 
{
}

void Rigaku::Process() 
{
  
}

ImmBlock* Rigaku::NextFrames(int count) {

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
