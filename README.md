#	LCQPow - A Solver for Quadratic Programs with Linear Complementarity Constraints

LCQPow is a open-source solver for Quadratic Programs with Complementarity Constraints. The approach is based on a standard penalty homotopy reformulated using sequential convex programming. The convex sequence derives from linearizing the (necessarily) nonconvex penalty function. This leads to a constant objective Hessian matrix throughout all iterates, and thus enables us to solve the linear complementarity quadratic program with a single factorization of the KKT matrix (by using qpOASES).

The entire strategy is presented in detail in [this paper](https://ieeexplore.ieee.org/abstract/document/9439931).

## Requirements
* Build process is currently only tested on Ubuntu >= 18.04
* CMake version >= 3.13.0
* This project depends on a few external repos. These are included and linked automatically:
   * [qpOASES](https://github.com/coin-or/qpOASES)
   * [OSQP](https://github.com/osqp/osqp)
   * [googletest](https://github.com/google/googletest)
   * [pybind11](https://github.com/pybind/pybind11)

## GETTING STARTED
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
You may specify the following flags (bracket represents default):
```
BUILD TYPE         [Release] / Debug
EXAMPLES           [ON] /  OFF
MATLAB INTERFACE   [ON] /  OFF
PYTHON INTERFACE   [ON] /  OFF
DOCUMENTATION      [ON] /  OFF
UNIT_TESTS         [ON] /  OFF
PROFILING           ON  / [OFF]
QPOASES_SCHUR       ON  / [OFF]
```

3. To test the build you can **run the examples** in the directory `<LCQPow-dir>/examples`.

4. Even more examples, in particular variations of the options, can be found in `<LCQPow-dir>/test/examples` directory. Those are not included in the examples directory in order to keep the example set neatly arranged.

## MATLAB Interface
The **MATLAB interface** is built automatically if matlab is successfully detected by CMake. Make sure that your **linker can locate the created libraries**, e.g. by exporting the library path in **the same shell as the one you call matlab in**:
```
$ export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:<LCQPow-dir>/build/lib
$ matlab
```
where `<LCQPow-dir>` represents the base directory of your local LCQPow repository.

Whenever **using the MATLAB interface** make sure that **MATLAB is able to locate the interface** by executing
```
addpath("<LCQPow-dir>/build/lib")
```

Navigate your MATLAB editor to the directory `<LCQPow-dir>/interfaces/matlab/examples` and play with any of the provided codes.

Type `help LCQPow` in order to obtain the documentation of our MATLAB interface.

## Python Interface
Thanks to a contribution by [Sotaro Katayama](https://github.com/mayataka) you can call the solver through its python interface. Doing so can be en/disabled by setting the respective cmake flag. The python interface requires the **Python 3 development** package as well as **Eigen3**. Make sure these are installed. On Ubuntu this can be achieved by running
```
apt install libeigen3-dev
apt-get install python3-dev
```

Remark: unlike the matlab interface this is more in an experimental stage.

## Sparse vs Dense
The most tested version of LCQPow uses qpOASES with dense linear algebra. There exist two alternatives:
  - using OSQP, which exploits sparsity naturally,
  - or using qpOASES Schur complement method (uses sparse linear solver MA57).

Usage of the qpOASES sparse method relies on some Matlab libraries (libmwma57.so, libmwlapack.so, libmwblas.so, libmwmetis.so), which are automatically detected and linked if they exist.

## License
The file LICENSE contains a copy of the GNU Lesser General Public License (v2.1). Please read it carefully before using LCQPow!

## CONTACT THE AUTHORS
If you have got questions, remarks or comments on LCQPow, it is strongly encouraged to report them by creating a new issue on this github page.

Finally, you may contact the main author directly:
        Jonas Hall,  hall.f.jonas@gmail.com

Also bug reports, source code enhancements or success stories are most welcome!


## Credits
The design of this software project is in large parts inspired by that of [qpOASES](https://github.com/coin-or/qpOASES). This includes the object-oriented design, and the specific classes introduces (though each of the classes are quite different from qpOASES itself). Some code snippets may have more similarity to the code of qpOASES than others (in particular constructors and destructors may be very similar). Additionally, some files contain one-to-one code snippets copied from qpOASES (e.g., the matlab interface contains code that parses qpOASES options). In such files the copyright header is included explicitly.

## Logo Design
Thank you Johanna Schmidt for designing this logo!