#!/bin/bash

# Break if any command fails
set -e

# Configure and build
mkdir -p build
cd build
cmake ..
make
cd ..

# Run all examples
for f in ./build/bin/examples/*; do
    echo "Running $f"
    valgrind --leak-check=full --error-exitcode=1 $f >> test/mem_examples.log
done

# Run all test examples
for f in ./build/bin/tests/*_TEST; do
    echo "Running $f"
    valgrind --leak-check=full --error-exitcode=1 $f >> test/mem_tests.log
done
