#!/bin/bash

mkdir -p build
cd build
rm -r *
cmake ..
make
./RunUnitTests
cd ..
