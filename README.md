#	LCQPow -- A Solver for Quadratic Programs with Linear Complementarity Constraints.

LCQPow is a open-source solver for Quadratic Programs with Complementarity Constraints. The approach is based on a standard penalty homotopy reformulated using sequential convex programming. The convex sequence derives from linearizing the (necessarily) nonconvex penalty function. This leads to a constant objective Hessian matrix throughout all iterates, and thus enables us to solve the linear complementarity quadratic program with a single factorization of the KKT matrix (by using qpOASES).

The entire strategy is presented in detail in [this paper](https://ieeexplore.ieee.org/abstract/document/9439931).

## Requirements
- Build process is currently only tested on Ubuntu 20.04
- CMake version >= VERSION 3.20.3
- all library dependencies are linked automatically (qpOASES, OSQP, googletest)
-

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

3. The **MATLAB interface** is built automatically if matlab is successfully detected by CMake.

## License
The file LICENSE contains a copy of the GNU Lesser General Public License (v2.1). Please read it carefully before using LCQPow!

## CONTACT THE AUTHORS
If you have got questions, remarks or comments on LCQPow, it is strongly encouraged to report them by creating a new issue on this github page.

Finally, you may contact the main author directly:
        Jonas Hall,  hall.f.jonas@gmail.com

Also bug reports, source code enhancements or success stories are most welcome!
