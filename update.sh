#!/bin/sh

git pull
[ -d build ] && rm -rf build
mkdir build
cd build
/local/cmake-3.9.0-rc4-Linux-x86_64/bin/cmake ../
make clean
make -j 4 VERBOSE=1

