#	LCQPow -- A Solver for Quadratic Programs with Linear Complementarity Constraints.

LCQPow is a open-source solver for Quadratic Programs with Complementarity Constraints. The approach is based on a standard penalty homotopy reformulated using sequential convex programming. The convex sequence derives from linearizing the (necessarily) nonconvex penalty function. This leads to a constant objective Hessian matrix throughout all iterates, and thus enables us to solve the linear complementarity quadratic program with a single factorization of the KKT matrix (by using qpOASES).

The entire strategy is presented in detail in [this paper](https://ieeexplore.ieee.org/abstract/document/9439931).

## Requirements
* Build process is currently only tested on Ubuntu 20.04
* CMake version >= VERSION 3.20.3
* All dependencies are linked automatically. These are
   * [qpOASES](https://github.com/coin-or/qpOASES)
   * [OSQP](https://github.com/osqp/osqp)
   * [googletest](https://github.com/google/googletest)

## GETTING STARTED (#getting-started)
1. **Clone the repository**, and recursively initialize the submodules

```
$ git clone https://github.com/hallfjonas/LCQPow.git
$ cd LCQPow
$ git submodule update --init --recursive
```

2. **Configure and build** the project
```
$ mkdir build
$ cd build
$ cmake ..
$ make
```

3. To test the build you can **run the test examples**
```
$ bin/examples/warm_up
$ bin/examples/warm_up_sparse
$ bin/examples/warm_up_osqp
$ bin/examples/OptimizeOnCircle
```

## MATLAB Interface
The **MATLAB interface** is built automatically if matlab is successfully detected by CMake. Make sure that your **linker can locate the created libraries** by exporting the library path in **the same shell as the one you call matlab in**:
```
$ export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:<LCQPow-dir>/build/lib
$ matlab
```
where `<LCQPow-dir>` represents the base directory of your local LCQPow repository.

Whenever **using the MATLAB interface** make sure that you have built the project according to the [getting started section](#getting-started), and further let MATLAB locate the interface by executing
```
addpath(<LCQPow-dir>/build/lib)
```

Navigate your MATLAB editor to the directory `<LCQPow-dir>/interfaces/matlab/examples` and play with any of the provided codes.

Type `help LCQPow` in order to obtain the documentation of our MATLAB interface.

## License
The file LICENSE contains a copy of the GNU Lesser General Public License (v2.1). Please read it carefully before using LCQPow!

## CONTACT THE AUTHORS
If you have got questions, remarks or comments on LCQPow, it is strongly encouraged to report them by creating a new issue on this github page.

Finally, you may contact the main author directly:
        Jonas Hall,  hall.f.jonas@gmail.com

Also bug reports, source code enhancements or success stories are most welcome!
