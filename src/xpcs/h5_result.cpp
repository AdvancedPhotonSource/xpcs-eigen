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

#include "h5_result.h"

#include "hdf5.h"

#include "configuration.h"

namespace xpcs {

void H5Result::write2DData(const std::string &file, 
                           const std::string &grpname,
                           const std::string &nodename,
                           Eigen::Ref<Eigen::MatrixXf> mat,
                           bool compression)
{
    hid_t file_id, exchange_grp_id, dataset_id, dataspace_id;
    hid_t plist_id;
    hsize_t dims[2];
    hsize_t cdims[2];

    file_id = H5Fopen(file.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    // Disable hdf std err printout for opening an non-existing group.
    H5Eset_auto(H5E_DEFAULT, NULL, NULL);
    exchange_grp_id = H5Gopen2(file_id, grpname.c_str(), H5P_DEFAULT);
    if (exchange_grp_id < 0) {
        exchange_grp_id = H5Gcreate(file_id, grpname.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    }

    //TODO :: Log error if the grp creation fail. 
    dataset_id = H5Dopen2(exchange_grp_id, nodename.c_str(), H5P_DEFAULT);

    if (dataset_id < 0) {

        dims[0] = mat.cols();
        dims[1] = mat.rows();

        dataspace_id = H5Screate_simple(2, dims, NULL);

        hid_t plist_id = H5Pcreate (H5P_DATASET_CREATE);

        cdims[0] = 20;
        cdims[0] = 20;
        H5Pset_chunk(plist_id, 2, cdims);

        H5Pset_deflate(plist_id, 6);

        dataset_id = H5Dcreate2(exchange_grp_id, nodename.c_str(), H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT, plist_id, H5P_DEFAULT);
    }

    hid_t stats = H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, mat.data());

    H5Dclose(dataset_id);
    if (dataspace_id) H5Sclose(dataspace_id);
    H5Gclose(exchange_grp_id);
    H5Fclose(file_id);    
}

void H5Result::write2DData(const std::string &file, 
                           const std::string &grpname,
                           const std::string &nodename,
                           int size0,
                           int size1,
                           float* data,
                           bool compression
                           )
{
    hid_t file_id, exchange_grp_id, dataset_id, dataspace_id;
    hsize_t dims[2];
    hsize_t cdims[2];

    file_id = H5Fopen(file.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    // Disable hdf std err printout for opening an non-existing group.
    H5Eset_auto(H5E_DEFAULT, NULL, NULL);
    exchange_grp_id = H5Gopen2(file_id, grpname.c_str(), H5P_DEFAULT);
    if (exchange_grp_id < 0) {
        exchange_grp_id = H5Gcreate(file_id, grpname.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    }

    //TODO :: Log error if the grp creation fail. 

    dataset_id = H5Dopen2(exchange_grp_id, nodename.c_str(), H5P_DEFAULT);

    if (dataset_id < 0) {

        printf("Compresing data\n");
        
        dims[0] = size0;
        dims[1] = size1;

        dataspace_id = H5Screate_simple(2, dims, NULL);
         
        hid_t plist_id = H5Pcreate (H5P_DATASET_CREATE);

        cdims[0] = 20;
        cdims[1] = 20;
        H5Pset_chunk(plist_id, 2, cdims);

        H5Pset_deflate(plist_id, 6);

        dataset_id = H5Dcreate2(exchange_grp_id, nodename.c_str(), H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT, plist_id, H5P_DEFAULT);
    }

    hid_t stats = H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);

    H5Dclose(dataset_id);
    if (dataspace_id) H5Sclose(dataspace_id);
    H5Gclose(exchange_grp_id);
    H5Fclose(file_id);    
}


void H5Result::write3DData(const std::string &file, 
                           const std::string &grpname,
                           const std::string &nodename,
                           int size0,
                           int size1,
                           int size2,
                           float* data)
{
    hid_t file_id, exchange_grp_id, dataset_id, dataspace_id;
    hsize_t dims[3];
    hsize_t start_3d[3];
    hsize_t stride_3d[3];
    hsize_t count_3d[3];


    file_id = H5Fopen(file.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    // Disable hdf std err printout for opening an non-existing group.
    H5Eset_auto(H5E_DEFAULT, NULL, NULL);
    exchange_grp_id = H5Gopen2(file_id, grpname.c_str(), H5P_DEFAULT);
    if (exchange_grp_id < 0) {
        exchange_grp_id = H5Gcreate(file_id, grpname.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    }

    //TODO :: Log error if the grp creation fail. 

    dataset_id = H5Dopen2(exchange_grp_id, nodename.c_str(), H5P_DEFAULT);

    if (dataset_id < 0) {

        dims[0] = size0;
        dims[1] = size1;
        dims[2] = size2;

        dataspace_id = H5Screate_simple(3, dims, NULL);

        start_3d[0] = 0;
        start_3d[1] = 0;
        start_3d[2] = 0;

        stride_3d[0] = 1;
        stride_3d[1] = 1;
        stride_3d[2] = 1;

        count_3d[0] = size0;
        count_3d[1] = size1;
        count_3d[2] = size2;

        H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, start_3d, stride_3d, count_3d, NULL);

        dataset_id = H5Dcreate(exchange_grp_id, nodename.c_str(), H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    }

    hid_t stats = H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);

    H5Dclose(dataset_id);
    if (dataspace_id) H5Sclose(dataspace_id);
    H5Gclose(exchange_grp_id);
    H5Fclose(file_id);    
}


void H5Result::write2DData(const std::string &file, 
                           const std::string &grpname,
                           const std::string &nodename,
                           int size0,
                           int size1,
                           double* data,
                           bool compression)
{
    hid_t file_id, exchange_grp_id, dataset_id, dataspace_id;
    hsize_t dims[2];

    file_id = H5Fopen(file.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    // Disable hdf std err printout for opening an non-existing group.
    H5Eset_auto(H5E_DEFAULT, NULL, NULL);
    exchange_grp_id = H5Gopen2(file_id, grpname.c_str(), H5P_DEFAULT);
    if (exchange_grp_id < 0) {
        exchange_grp_id = H5Gcreate(file_id, grpname.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    }

    //TODO :: Log error if the grp creation fail. 

    dataset_id = H5Dopen2(exchange_grp_id, nodename.c_str(), H5P_DEFAULT);

    if (dataset_id < 0) {

        dims[0] = size0;
        dims[1] = size1;

        dataspace_id = H5Screate_simple(2, dims, NULL);
        dataset_id = H5Dcreate(exchange_grp_id, nodename.c_str(), H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    }

    hid_t stats = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);

    H5Dclose(dataset_id);
    if (dataspace_id) H5Sclose(dataspace_id);
    H5Gclose(exchange_grp_id);
    H5Fclose(file_id);    
}




void H5Result::write1DData(const std::string &file, 
                           const std::string &grpname,
                           const std::string &nodename,
                           Eigen::Ref<Eigen::VectorXf> vec)
{
    hid_t file_id, exchange_grp_id, dataset_id, dataspace_id;
    hsize_t dims[2];

    file_id = H5Fopen(file.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    // Disable hdf std err printout for opening an non-existing group.
    H5Eset_auto(H5E_DEFAULT, NULL, NULL);
    exchange_grp_id = H5Gopen2(file_id, grpname.c_str(), H5P_DEFAULT);
    if (exchange_grp_id < 0) {
        exchange_grp_id = H5Gcreate(file_id, grpname.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    }

    //TODO :: Log error if the grp creation fail. 

    dataset_id = H5Dopen2(exchange_grp_id, nodename.c_str(), H5P_DEFAULT);

    if (dataset_id < 0) {

        dims[0] = 1;
        dims[1] = vec.rows();

        dataspace_id = H5Screate_simple(2, dims, NULL);
        dataset_id = H5Dcreate(exchange_grp_id, nodename.c_str(), H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    }

    hid_t stats = H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, vec.data());

    H5Dclose(dataset_id);
    if (dataspace_id) H5Sclose(dataspace_id);
    H5Gclose(exchange_grp_id);
    H5Fclose(file_id);    
}


void H5Result::write1DData(const std::string &file, 
                           const std::string &grpname,
                           const std::string &nodename,
                           int size,
                           float* data)
{
    hid_t file_id, exchange_grp_id, dataset_id, dataspace_id;
    hsize_t dims[2];

    file_id = H5Fopen(file.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    // Disable hdf std err printout for opening an non-existing group.
    H5Eset_auto(H5E_DEFAULT, NULL, NULL);
    exchange_grp_id = H5Gopen2(file_id, grpname.c_str(), H5P_DEFAULT);
    if (exchange_grp_id < 0) {
        exchange_grp_id = H5Gcreate(file_id, grpname.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    }

    //TODO :: Log error if the grp creation fail. 

    dataset_id = H5Dopen2(exchange_grp_id, nodename.c_str(), H5P_DEFAULT);

    if (dataset_id < 0) {

        dims[0] = 1;
        dims[1] = size;

        dataspace_id = H5Screate_simple(2, dims, NULL);
        dataset_id = H5Dcreate(exchange_grp_id, nodename.c_str(), H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    }

    hid_t stats = H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);

    H5Dclose(dataset_id);
    if (dataspace_id) H5Sclose(dataspace_id);
    H5Gclose(exchange_grp_id);
    H5Fclose(file_id);    
}


void H5Result::writePixelSum(const std::string &file, 
                   const std::string &grpname,
                   Eigen::Ref<Eigen::VectorXf> pixelSum) {

    hid_t file_id, exchange_grp_id, dataset_id, dataspace_id;
    hsize_t dims[2];

    file_id = H5Fopen(file.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    // Disable hdf std err printout for opening an non-existing group.
    H5Eset_auto(H5E_DEFAULT, NULL, NULL);
    exchange_grp_id = H5Gopen2(file_id, grpname.c_str(), H5P_DEFAULT);
    if (exchange_grp_id < 0) {
        exchange_grp_id = H5Gcreate(file_id, grpname.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    }

    //TODO :: Log error if the grp creation fail. 
    dataset_id = H5Dopen2(exchange_grp_id, "pixelSum", H5P_DEFAULT);

    Configuration *conf = Configuration::instance();
    int w = conf->getFrameWidth();
    int h = conf->getFrameHeight();

    if (dataset_id < 0) {

        dims[0] = h;
        dims[1] = w;

        dataspace_id = H5Screate_simple(2, dims, NULL);
        dataset_id = H5Dcreate(exchange_grp_id, "pixelSum", H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    }

    hid_t stats = H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, pixelSum.data());

    H5Dclose(dataset_id);
    if (dataspace_id) H5Sclose(dataspace_id);
    H5Gclose(exchange_grp_id);
    H5Fclose(file_id);
}


void H5Result::writeFrameSum(const std::string &file, 
                   const std::string &grpname,
                   Eigen::Ref<Eigen::VectorXf> frameSum) {

    hid_t file_id, exchange_grp_id, dataset_id, dataspace_id;
    hsize_t dims[2];

    file_id = H5Fopen(file.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    // Disable hdf std err printout for opening an non-existing group.
    H5Eset_auto(H5E_DEFAULT, NULL, NULL);
    exchange_grp_id = H5Gopen2(file_id, grpname.c_str(), H5P_DEFAULT);
    if (exchange_grp_id < 0) {
        exchange_grp_id = H5Gcreate(file_id, grpname.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    }

    //TODO :: Log error if the grp creation fail. 
    dataset_id = H5Dopen2(exchange_grp_id, "frameSum", H5P_DEFAULT);

    Configuration *conf = Configuration::instance();

    if (dataset_id < 0) {

        dims[0] = 1;
        dims[1] = conf->getFrameTodoCount();

        dataspace_id = H5Screate_simple(2, dims, NULL);
        dataset_id = H5Dcreate(exchange_grp_id, "frameSum", H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    }

    hid_t stats = H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, frameSum.data());

    H5Dclose(dataset_id);
    if (dataspace_id) H5Sclose(dataspace_id);
    H5Gclose(exchange_grp_id);
    H5Fclose(file_id);
}

}
