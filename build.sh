#!/bin/bash

mkdir -p build
cd build
cmake -DDLONG=OFF -DBUILD_SHARED_LIBS=ON ..
make
cd ..
