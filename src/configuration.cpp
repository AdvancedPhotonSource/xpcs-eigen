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

#include "configuration.h"
#include "hdf5.h"

#include <iostream>

#include <boost/algorithm/string.hpp>

using namespace std;

Configuration *Configuration::s_instance = 0;


Configuration* Configuration::instance()
{
    if (!s_instance)
        s_instance = new Configuration();

    return s_instance;
}

Configuration::Configuration() 
{

}

Configuration::~Configuration()
{
}

void Configuration::init(const string &path)
{
    //TODO: Force initialization to once per instantiation
    m_filename = path;

    file_id = H5Fopen(path.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    std::string value = getString("/xpcs/compression");

    if (boost::iequals(value, "ENABLED"))
        this->compression = true;
    
    this->dqmap = get2DTable("/xpcs/dqmap");
    this->sqmap = get2DTable("/xpcs/sqmap");

    this->xdim = getInteger("/measurement/instrument/detector/x_dimension");
    this->ydim = getInteger("/measurement/instrument/detector/y_dimension");
    this->frameCount = getInteger("/xpcs/data_end_todo");


    this->m_totalStaticPartitions = 0;
    this->m_totalDynamicPartitions = 0;


    this->m_validPixelMask = new short[this->xdim * this->ydim];
    // Build mapping from q-bin to s-bins and from s-bins to pixels. 
    for (int i=0; i<(xdim*ydim); i++) {
        if (dqmap[i] < 1 || sqmap[i] < 1) continue;

        // Mark this pixel to be part of a valid mask. 
        m_validPixelMask[i] = 1;

        map<int, map<int, vector<int>> >::iterator it = m_mapping.find(dqmap[i]);

        // Highest q-map value is equal to total number of dynamic partitions. 
        if (this->m_totalDynamicPartitions < dqmap[i])
            this->m_totalDynamicPartitions = dqmap[i];

        if (it != m_mapping.end()) {
            map<int, vector<int> > &v = it->second;
            
            map<int, vector<int>>::iterator it2 = v.find(sqmap[i]);
            
            if (this->m_totalStaticPartitions < sqmap[i])
                this->m_totalStaticPartitions = sqmap[i];

            if (it2 != v.end()) {
                vector<int> &v2 = it2->second;
                v2.push_back(i);
            } else {
                vector<int> data;
                data.push_back(i);
                v[sqmap[i]] = data;
            }
        }
        else {
            map<int, vector<int> > mapping;
            vector<int> data;
            data.push_back(i);
            mapping[sqmap[i]] = data;
            m_mapping[dqmap[i]] = mapping;
        }
    }

    H5Fclose(file_id);
}

/**
 * Data is stored in row major format. 
 */
int* Configuration::get2DTable(const std::string &path)
{
    hid_t dataset_id;

    dataset_id = H5Dopen(this->file_id, path.c_str(), H5P_DEFAULT);
    hid_t space = H5Dget_space(dataset_id);

    hsize_t dims[2] = {0, 0};
    int ndims = H5Sget_simple_extent_dims (space, dims, NULL);

    int *data = new int[dims[0] * dims[1]];

    H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);

    H5Dclose(dataset_id);

    return data;
}

std::string Configuration::getString(const std::string &path)
{
    hid_t  dataset_id;
    herr_t status;

    dataset_id = H5Dopen(this->file_id, path.c_str(), H5P_DEFAULT);
    
    hid_t dtype = H5Dget_type(dataset_id);
    hid_t size = H5Dget_storage_size(dataset_id);

    char * str = new char[size + 1];
    str[size] = '\0';

    status = H5Dread(dataset_id, dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, str);

    H5Dclose(dataset_id);

    return std::string(str);
}

int Configuration::getInteger(const std::string &path)
{
    hid_t  dataset_id;
    herr_t status;

    dataset_id = H5Dopen(this->file_id, path.c_str(), H5P_DEFAULT);
    
    hid_t dtype = H5Dget_type(dataset_id);
    hid_t size = H5Dget_storage_size(dataset_id);

    int value = 0;

    status = H5Dread(dataset_id, dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &value);

    H5Dclose(dataset_id);

    return value;
}

int* Configuration::getDQMap()
{
    return this->dqmap;
}

int* Configuration::getSQMap()
{
    return this->sqmap;
}

int Configuration::getFrameWidth()
{
    return this->xdim;
}

int Configuration::getFrameHeight()
{
    return this->ydim;
}

int Configuration::getFrameCount()
{
    return this->frameCount;
}

std::map<int, std::map<int, std::vector<int>> > Configuration::getBinMaps()
{
    return this->m_mapping;
}


int Configuration::getTotalStaticPartitions()
{
    return this->m_totalStaticPartitions;
}

int Configuration::getTotalDynamicPartitions()
{
    return this->m_totalDynamicPartitions;
}

std::string Configuration::getFilename()
{
    return this->m_filename;
}

short* Configuration::getPixelMask()
{
    return this->m_validPixelMask;
}
