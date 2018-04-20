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

#include "dense_average.h"

#include <set>
#include <math.h>
#include <stdio.h>
#include <iostream>

#include "xpcs/configuration.h"
#include "xpcs/io/imm_reader.h"
#include "xpcs/data_structure/sparse_data.h"

namespace xpcs {
namespace filter {

DenseAverage::DenseAverage() {
  Configuration *conf = Configuration::instance();
  pixels_ = conf->getFrameWidth() * conf->getFrameHeight();
  average_size_ = conf->FrameAverage();
  stride_size_ = conf->FrameStride();
  pixels_value_ = new float[pixels_];
}

DenseAverage::~DenseAverage() {

}

void DenseAverage::Apply(struct xpcs::io::ImmBlock* blk) {
  int **indx = blk->index;
  float **val = blk->value;
  int frames = blk->frames;
  std::vector<int> ppf = blk->pixels_per_frame;

  if (frames < 2) return;

  for (int j = 0; j < ppf[0]; j++) {
    int px = j;
    float v = val[0][j];
    pixels_value_[px] = v; 
  }

  for (int i = stride_size_ ; i < frames; i+=stride_size_) {
    for (int j = 0; j < ppf[i]; j++) {
      int px = j;
      float v = val[i][j];
      pixels_value_[px] += v; 
      // if (v > 0)
      //   printf("%d - %f\n", px, v);
    }
  }

  int **new_index = new int*[1];
  float **new_val = new float*[1];
  int new_frames = 1;
  new_index[0] = new int[ppf[0]];
  new_val[0] = new float[ppf[0]];
  std::vector<int> new_ppf = {ppf[0]};

  int ind = 0;
  for (int i = 0; i < ppf[0]; i++) {
    int px = i;
    new_index[0][ind] = px;
    new_val[0][ind] = pixels_value_[px] / average_size_;
    ind++;
  }

  blk->index = new_index;
  blk->value = new_val;
  blk->frames = new_frames;
  blk->pixels_per_frame = new_ppf;

  //TODO smart pointers to handle memory
}

} // namespace io
} // namespace xpcs