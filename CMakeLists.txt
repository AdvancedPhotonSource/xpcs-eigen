project(xpcs-eigen)

cmake_minimum_required(VERSION 3.5)

add_definitions(-std=c++0x)

set(sources
    src/configuration.cpp
    src/corr.cpp
    src/imm.cpp
    src/main.cpp)

# TODO. Fix the hard-coded path. 
include_directories("/anaconda/pkgs/hdf5-1.8.17-9/include" "/Users/faisal/Development/xpcs-eigen/boost_1_61_0/")
link_directories("/anaconda/pkgs/hdf5-1.8.17-9/lib")

add_executable(corr ${sources})

target_link_libraries(corr hdf5)