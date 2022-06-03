#!/bin/bash

# Break if any command fails
set -e

# Configure and build
mkdir -p build
cd build
cmake ..
make
cd ..

# Run all test examples
for f in ./build/bin/tests/*; do
    echo "Running $f"
    $f
done

# Run unit tests
build/bin/tests/RunUnitTests
