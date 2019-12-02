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
#include "hdf5.h"

#include <stdio.h>
#include <iostream>
#include <vector>
#include <iterator>

#include "xpcs/configuration.h"

namespace xpcs {
namespace io {

Hdf5::Hdf5(const std::string& filename) {
    xpcs::Configuration *conf = xpcs::Configuration::instance();
    int fw = conf->getFrameWidth();
    int fh = conf->getFrameHeight();

    hsize_t count[3] = {1, fw, fh};
    hid_t file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    hid_t gid = H5Gopen2(file, "/entry/data", H5P_DEFAULT);
    dataset_id_ = H5Dopen1(gid, "data");
    space_id_ = H5Dget_space(dataset_id_);
    memspace_id_ = H5Screate_simple(3, count, NULL);

    buffer_ = new unsigned short[fw * fh];
}

Hdf5::~Hdf5() {
}

ImmBlock* Hdf5::NextFrames(int count) {
    xpcs::Configuration *conf = xpcs::Configuration::instance();
    int fw = conf->getFrameWidth();
    int fh = conf->getFrameHeight();

    int **index = new int*[count];
    float **value = new float*[count];
    double *clock = new double[count];
    double *ticks = new double[count];

    int done = 0;
    std::vector<int> ppf;

    hsize_t hdf_count[3] = {1, fw, fh};

    while (done < count) {
      hsize_t offset[3] = {last_frame_index_, 0, 0};
      herr_t errstatus = H5Sselect_hyperslab(space_id_, H5S_SELECT_SET, offset, NULL, hdf_count, NULL);
      hid_t status = H5Dread(dataset_id_, H5T_NATIVE_UINT16, memspace_id_, space_id_, H5P_DEFAULT, buffer_);

      // count total non-zero values.
      int sparse = 0;
      for (int i = 0; i < (fw * fh); i++) {
        if (buffer_[i] != 0) sparse++;
      }

      index[done] = new int[sparse];
      value[done] = new float[sparse];

      int idx = 0;
      for (int i = 0; i < (fw * fh); i++) {
        if (buffer_[i] != 0) {
          index[done][idx] = i;
          value[done][idx] = buffer_[i];

          idx++;
        }
      }
      ppf.push_back(sparse);

      clock[done] = last_frame_index_;
      ticks[done] = last_frame_index_;
      
      last_frame_index_++;
      done++;
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

void Hdf5::SkipFrames(int count) {
  last_frame_index_ += count;
}

void Hdf5::Reset() {
  last_frame_index_ = 0;
}

bool Hdf5::compression() { return false; }

} // namespace io
} // namespace xpcs
