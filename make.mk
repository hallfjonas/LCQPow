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

## user configuration (adapt to your local settings)
# 1) qpOASES
QPOASES_IDIR   = /usr/local/qpOASES/include
QPOASES_LIB_DIR = /usr/local/qpOASES/bin
QPOASES_LINK   = -L${QPOASES_LIB_DIR} -Wl,-rpath=${QPOASES_LIB_DIR} -lqpOASES 

# 2) CasADi
CASADI_IDIR   = /usr/local/casadi
CASADI_LINK = -L/usr/local/casadi/build/lib -lcasadi


## Do not touch this
# include directories, relative
TOP = $(realpath $(dir $(lastword $(MAKEFILE_LIST))))
IDIR =   ${TOP}/include
SRCDIR = ${TOP}/src
BINDIR = ${TOP}/bin
DEBUGDIR = ${TOP}/debug

# Compiler flags
CPP = g++
CC  = gcc
AR  = ar
RM  = rm
F77 = gfortran
ECHO = echo
CD = cd
CP = cp
CPPFLAGS = -Wall -pedantic -Wshadow -Wfloat-equal -O3 -Wconversion -Wsign-conversion -fPIC -DLINUX -D__USE_LONG_INTEGERS__ -D__USE_LONG_FINTS__ -D${DEF_SOLVER} -D__NO_COPYRIGHT__

# file extensions
OBJEXT = o
LIBEXT = a
DLLEXT = so
EXE =
MEXOCTEXT = mex
DEF_TARGET = -o $@
SHARED = -shared

# Links to libraries
LINK_LIBRARIES = ${QPOASES_LINK} ${CASADI_LINK}
LINK_DEPENDS = ${BINDIR}/liblcqpOASES.${LIBEXT} ${BINDIR}/liblcqpOASES.${DLLEXT}


##
##	end of file
##
