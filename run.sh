#!/bin/bash

# Break if any command fails
set -e

# Run unit tests
# build/bin/tests/RunUnitTests

# Run all examples
build/bin/examples/warm_up
build/bin/examples/warm_up_schur
build/bin/examples/warm_up_osqp
build/bin/examples/warm_up_w_A
build/bin/examples/OptimizeOnCircle
