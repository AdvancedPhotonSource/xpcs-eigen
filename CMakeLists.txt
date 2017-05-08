project(xpcs-eigen)

cmake_minimum_required(VERSION 3.5)

add_definitions(-std=c++0x)

set(sources
    src/configuration.cpp
    src/corr.cpp
    src/imm.cpp
    src/main.cpp)

# TODO. Fix the hard-coded path. 
include_directories("/local/fkhan/libexec/include")
link_directories("/local/fkhan/libexec/lib")

add_executable(corr ${sources})

target_link_libraries(corr h