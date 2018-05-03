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
#ifndef XPCS_DENSE_FILTER_H
#define XPCS_DENSE_FILTER_H

#include "filter.h"

namespace xpcs {

namespace data_structure {
  class SparseData;
  class DarkImage;
}

namespace io {
  struct ImmBlock;
}

namespace filter {  

class DenseFilter : public Filter  {

public:
  DenseFilter(xpcs::data_structure::DarkImage* dark_image);

  ~DenseFilter();

  void Apply(struct xpcs::io::ImmBlock* data);

  float* PixelsSum();

  float* FramesSum();

  float* PartitionsMean();

  float* PartialPartitionsMean();

  float* TimestampClock();

  float* TimestampTicks();

  xpcs::data_structure::SparseData* Data();

private:

  xpcs::data_structure::SparseData *data_;

  xpcs::data_structure::DarkImage *dark_image_;

  short *pixel_mask_;

  int *sbin_mask_;

  float sigma_;

  float threshold_;

  float *pixels_value_;

  float *pixels_sum_;

  float *frames_sum_;

  float *partial_partitions_mean_;

  float *partitions_mean_;

  float *timestamp_clock_;

  float *timestamp_ticks_;

  double *flatfield_;

  int frame_width_;

  int frame_height_; 

  int static_window_; 

  int total_static_partns_; 

  int frame_index_;

  int global_frame_index_;

  int frames_todo_;

  int real_frames_todo_;
  
  int swindow_;

  int partition_no_;

  int stride_size_;

  int average_size_;

  float eff_; 

  float det_adhu_; 

  float preset_; 

  float normFactor_; 

};

} //namespace filter
} //namespace xpcs

#endif