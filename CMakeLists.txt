project(xpcs-eigen)

cmake_minimum_required(VERSION 3.5)

#set(CMAKE_BUILD_TYPE RELEASE)
set(CMAKE_BUILD_TYPE DEBUG)

add_definitions(-std=c++0x)

set(headers
    src/xpcs/benchmark.h
    src/xpcs/configuration.h
    src/xpcs/corr.h
    src/xpcs/dark_image.h
    src/xpcs/funcs.h
    src/xpcs/h5_result.h
    src/xpcs/ds/row.h
    src/xpcs/ds/sparse_data.h
    src/xpcs/filter/filter.h
    src/xpcs/filter/sparse_filter.h
    src/xpcs/io/imm_header.h
   
)

set(sources
    src/xpcs/configuration.cpp
    src/xpcs/h5_result.cpp
    src/xpcs/corr.cpp
    src/xpcs/funcs.cpp
    src/xpcs/dark_image.cpp
    src/xpcs/ds/sparse_data.cpp
    src/xpcs/filter/sparse_filter.cpp
    src/xpcs/io/imm_reader.cpp
    src/xpcs/main2.cpp
)

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()


# TODO. Fix the hard-coded path. 
include_directories("/usr/local/hdf5/include" 
                    "src" 
                    "src/external")

link_directories("/usr/local/hdf5/lib")

add_subdirectory("src/external/gflags")
add_executable(corr ${sources})

target_link_libraries(corr hdf5 gflags)
