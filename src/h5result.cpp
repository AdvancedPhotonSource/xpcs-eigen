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

#include "h5result.h"

#include "hdf5.h"

void H5Result::writeG2(const std::string &file, 
                       const std::string &grpname,
                       Eigen::Ref<Eigen::MatrixXf> g2) {

    hid_t file_id, exchange_grp_id, dataset_id, dataspace_id;
    hsize_t dims[3];

    file_id = H5Fopen(file.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    // Disable hdf std err printout for opening an non-existing group. 
    H5Eset_auto(H5E_DEFAULT, NULL, NULL);
    exchange_grp_id = H5Gopen2(file_id, grpname.c_str(), H5P_DEFAULT);
    if (exchange_grp_id < 0) {
        exchange_grp_id = H5Gcreate(file_id, grpname.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    }

    //TODO :: Log error if the grp creation fail. 

    dataset_id = H5Dopen2(exchange_grp_id, "g2", H5P_DEFAULT);

    if (dataset_id < 0) {

        dims[0] = g2.rows();
        dims[1] = g2.cols();
        dims[2] = 1;

        dataspace_id = H5Screate_simple(3, dims, NULL);
        dataset_id = H5Dcreate(exchange_grp_id, "g2", H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    }

    hid_t stats = H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, g2.data());

    H5Dclose(dataset_id);
    if (dataspace_id) H5Sclose(dataspace_id);
    H5Gclose(exchange_grp_id);
    H5Fclose(file_id);
}