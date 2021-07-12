#!/bin/bash

# Break if any command fails
set -e

# Run examples
for f in ./build/bin/examples/*; do
    echo "Running $f"
    $f
done

# Run unit tests
build/bin/tests/RunUnitTests
