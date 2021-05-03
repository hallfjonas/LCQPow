##
##	This file is part of lcqpOASES.
##
##	lcqpOASES -- A Solver for Quadratic Programs with Commplementarity Constraints.
##	Copyright (C) 2020 - 2021 by Jonas Hall et al.
##
##	lcqpOASES is free software; you can redistribute it and/or
##	modify it under the terms of the GNU Lesser General Public
##	License as published by the Free Software Foundation; either
##	version 2.1 of the License, or (at your option) any later version.
##
##	lcqpOASES is distributed in the hope that it will be useful,
##	but WITHOUT ANY WARRANTY; without even the implied warranty of
##	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
##	See the GNU Lesser General Public License for more details.
##
##	You should have received a copy of the GNU Lesser General Public
##	License along with lcqpOASES; if not, write to the Free Software
##	Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
##

## User Configuration (adapt to your local settings)
## NOTE: We assume all external libraries to be placed into /usr/local/lib, such that no further linking is required.
##    	If you desire to link an external library from elsewhere,
## 		please do so by adding the path via -L/path-to-external-lib-dir -lexternal
## 		In this case, make sure that your system can find this library (e.g. by setting the LD_LIBRARY_PATH)

## Do not touch this
# include directories, relative
TOP = $(realpath $(dir $(lastword $(MAKEFILE_LIST))))
IDIR =   ${TOP}/include
SRCDIR = ${TOP}/src
BINDIR = ${TOP}/bin
LIBDIR = ${TOP}/lib
DEBUGDIR = ${TOP}/debug
BUILDDIR = ${TOP}/build

# qpOASES
QPOASES_IDIR   = ${TOP}/external/qpOASES/include
QPOASES_LINK   = -lqpOASES

# OSQP
OSQP_IDIR   = ${TOP}/external/osqp/include
OSQP_LINK   = -losqp

# Compiler flags
CPP = g++
CC  = gcc
AR  = ar
RM  = rm
F77 = gfortran
ECHO = echo
CD = cd
CP = cp
CMAKE = cmake
MKDIR = mkdir
CPPFLAGS = -Wall -pedantic -Wshadow -Wfloat-equal -O3 -Wconversion -Wsign-conversion -fPIC -DLINUX -D__USE_LONG_INTEGERS__ -D__NO_COPYRIGHT__

# file extensions
OBJEXT = o
LIBEXT = a
DLLEXT = so
EXE =
MEXOCTEXT = mex
DEF_TARGET = -o $@
SHARED = -shared

# Links to libraries
LCQPOASES_LIB = ${LIBDIR}/liblcqpOASES.${LIBEXT} ${LIBDIR}/liblcqpOASES.${DLLEXT}
LCQPOASES_LINK = -L${LIBDIR} -llcqpOASES

LINK_LIBRARIES = ${QPOASES_LINK} ${OSQP_LINK}

##
##	end of file
##
