#!/bin/bash

# Break if any command fails
set -e

# Configure and build
mkdir -p build
cd build
cmake ..
make
cd ..
