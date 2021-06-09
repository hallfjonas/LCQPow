#!/bin/bash

mkdir -p build
cd build
cmake -DSOLVER_MA57 ..
make
cd ..
