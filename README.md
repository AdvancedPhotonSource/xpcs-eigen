[![Build Status](https://travis-ci.org/AdvancedPhotonSource/xpcs-eigen.svg?branch=master)](https://travis-ci.org/AdvancedPhotonSource/xpcs-eigen)


## X-Ray Photon Correlation Spectroscopy.

### Summary
---

XPCS (X-ray Photon Correlation Spectroscopy) is experimental technique used at
Sector 8 of APS. The main principle behind the technique is to study the
dynamics in materials at nano-scale by studying correlations in time series of
area detector images. These correlations involve analyzing the pixel-by-pixel
correlations for different time intervals. A multi-tau algorithm is commonly
used that defines the time intervals in a logarithmic manner. The current
state-of-the-art detector at Sector 8 is capable of acquiring megapixels of
frames at a rate of 60 Hz (120 MB/sec). A suite of CCD detectors are used with
data rates ranging from 4 MB/sec to 120 MB/sec. The next generation Fast CCD2
detector being built in collaboration between ANL and LBL is projected to run
at 200 Hz generating 400 MB/sec of data. The challenge is to compute these
correlations at a decent rate that comes within a small factor of the rate of
the data acquisition. Doing the analysis at (near) real time saves a lot of
effort and time of the scientists.

### Implementation
---

This C++ implementation , currently a work in progress, is based on the Eigen
(http://eigen.tuxfamily.org/index.php?title=Main_Page) library. We are handling
both sparse and non-sparse data using the equivalent constructs from the Eigen
library. 

### Compiling
---

You will need to have a latest (>= 3.5) version of```cmake``` installed in
order to compile this project. Additionally, you need the HDF5 C library.
Please follwing the instruction here
(https://support.hdfgroup.org/HDF5/release/obtain518.html) for installing it. 

Next, follow these steps:


* Install using ```apt```, ```yum```, ```pacman``` or your favourite package manager:
    1. eigen3
    2. hdf5
    3. gflags
    4. spdlog
* Create a build directory under the root of the project 
  ```mkdir build && cd build```
* ```cmake ../``` 
* ```make -j```


This will generate the binary executable ```corr```. 

### Running
----

The ```corr`` executable accepts two input arguments. The HDF5 containing the
configuration for the analysis job and an `IMM` file containing the binary
data. 

```
./corr configuration.hdf5 data.imm
```

