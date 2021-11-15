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
#include<fstream>

#include "xpcs/configuration.h"

using namespace std;

namespace xpcs {
namespace io {

Hdf5::Hdf5(const std::string& filename, bool tranposed) : tranposed_(tranposed) {
    xpcs::Configuration *conf = xpcs::Configuration::instance();
    int fw = conf->getFrameWidth();
    int fh = conf->getFrameHeight();
    hsize_t count[3] = {1, fw, fh};
    hid_t file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    hid_t gid = H5Gopen2(file, "/entry/data", H5P_DEFAULT);
    
    if (gid < 0) {
      std::cout<<"Failed to open group /entry/data" <<std::endl;
      return;
    }

    dataset_id_ = H5Dopen2(gid, "data", H5P_DEFAULT);
    
    if (dataset_id_ < 0) 
    {
      std::cout<<"Failed to open dataset data"<<std::endl;
      return;
    }

    datatype_id_ = H5Dget_type(dataset_id_);
    // std::cout<<"Data type " << datatype_id_ << " " <<  H5T_NATIVE_UINT16 <<std::endl;
    space_id_ = H5Dget_space(dataset_id_);
    memspace_id_ = H5Screate_simple(3, count, NULL);
    int rank = H5Sget_simple_extent_ndims(space_id_);
    //TODO: Check return values.

    hsize_t maxdims[3];
    H5Sget_simple_extent_dims(space_id_, dims_, maxdims);

    buffer_ = new unsigned short[fw * fh];

    // std::cout<<file<<std::endl;
    // std::cout<<space_id_<<std::endl;
    // std::cout<<memspace_id_<<std::endl;
    // std::cout<<dataset_id_<<std::endl;

    // unsigned int flags;
    // size_t nelmts = 1;
    // unsigned int values_out[1] = {99};
    // char filter_name[80];

    // hid_t dcpl = H5Dget_create_plist (dataset_id_);
    // H5Pget_filter2 (dcpl, (unsigned) 0, &flags, &nelmts, values_out, sizeof(filter_name), filter_name, NULL);

    // std::cout<<filter_name<<std::endl;

    // std::cout<<"Rank : " <<rank<<std::endl;

    last_frame_index_ = 0;
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

    hsize_t hdf_count[3] = {1, dims_[1], dims_[2]};

    while (done < count) {
      hsize_t offset[3] = {last_frame_index_, 0, 0};
      herr_t errstatus = H5Sselect_hyperslab(space_id_, H5S_SELECT_SET, offset, NULL, hdf_count, NULL);
      hid_t status = H5Dread(dataset_id_, datatype_id_, memspace_id_, space_id_, H5P_DEFAULT, buffer_);

      // unsigned short* buffer2 = new unsigned short[fw * fh];
      // for (int i = 0; i < (fw * fh); i++) {
      //   int new_index = (i % fh) * fh + (i / fh);
      //   buffer2[new_index] = buffer_[i];
      // }
      
      // ofstream wf("frame2.txt");
      // wf.write((char*)buffer2, (fw*fh)*sizeof(unsigned short));
      // wf.close();
      // exit(1);

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
          int pixel_index = i;
          // if (tranposed_) {
          //   int row = i % fh;
          //   int col = i / fh;
          //   pixel_index = row * fw + col;
          // }
          index[done][idx] = pixel_index;
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

bool Hdf5::compression() { return true; }

} // namespace io
} // namespace xpcs
