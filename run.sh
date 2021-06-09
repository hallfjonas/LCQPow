#!/bin/bash

# Break if any command fails
set -e

# Run unit tests
# build/bin/tests/RunUnitTests
for f in ./build/bin/examples/*; do
    echo "Running $f"
    $f
done
