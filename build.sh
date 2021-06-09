#!/bin/bash

mkdir -p build
cd build
cmake -DBUILD_SHARED_LIBS=ON ..
make
cd ..
