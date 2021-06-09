#!/bin/bash

mkdir -p build
cd build
cmake -DCMAKE_BUILD_TYPE=$1 -DQPOASES_SCHUR_COMPLEMENT=$2 ..
make
cd ..
