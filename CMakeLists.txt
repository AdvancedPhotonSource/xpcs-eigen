project(xpcs-eigen)

cmake_minimum_required(VERSION 3.5)

#set(CMAKE_BUILD_TYPE RELEASE)
set(CMAKE_BUILD_TYPE DEBUG)

add_definitions(-fno-stack-protector)
add_definitions(-std=c++0x)
add_definitions(-fno-stack-protector)


set(headers
    src/xpcs/benchmark.h
    src/xpcs/configuration.h
    src/xpcs/corr.h
    src/xpcs/data_structure/dark_image.h
    src/xpcs/funcs.h
    src/xpcs/h5_result.h
    src/xpcs/data_structure/row.h
    src/xpcs/data_structure/sparse_data.h
    src/xpcs/io/reader.h
    src/xpcs/io/imm.h
    src/xpcs/io/ufxc.h
    src/xpcs/io/rigaku.h
    src/xpcs/io/hdf5.h
    src/xpcs/filter/filter.h
    src/xpcs/filter/sparse_filter.h
    src/xpcs/filter/dense_filter.h
    src/xpcs/filter/stride.h
    src/xpcs/filter/average.h
    src/xpcs/filter/dense_average.h
)

set(sources
    src/xpcs/configuration.cpp
    src/xpcs/h5_result.cpp
    src/xpcs/corr.cpp
    src/xpcs/funcs.cpp
    src/xpcs/data_structure/dark_image.cpp
    src/xpcs/data_structure/sparse_data.cpp
    src/xpcs/io/imm.cpp
    src/xpcs/io/ufxc.cpp
    src/xpcs/io/rigaku.cpp
    src/xpcs/io/hdf5.cpp
    src/xpcs/filter/sparse_filter.cpp
    src/xpcs/filter/dense_filter.cpp
    src/xpcs/filter/stride.cpp
    src/xpcs/filter/average.cpp
    src/xpcs/filter/dense_average.cpp
    src/xpcs/main.cpp
)

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()


find_package(HDF5 REQUIRED COMPONENTS CXX HL)
include_directories( ${HDF5_INCLUDE_DIR} )

# TODO. Fix the hard-coded path. 
include_directories( "src" "src/external")

add_subdirectory("src/external/gflags")

add_executable(corr ${sources})

target_link_libraries(corr ${HDF5_LIBRARIES} gflags)

