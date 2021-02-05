##
##	This file is part of qpOASES.
##
##	qpOASES -- An Implementation of the Online Active Set Strategy.
##	Copyright (C) 2007-2017 by Hans Joachim Ferreau, Andreas Potschka,
##	Christian Kirches et al. All rights reserved.
##
##	qpOASES is free software; you can redistribute it and/or
##	modify it under the terms of the GNU Lesser General Public
##	License as published by the Free Software Foundation; either
##	version 2.1 of the License, or (at your option) any later version.
##
##	qpOASES is distributed in the hope that it will be useful,
##	but WITHOUT ANY WARRANTY; without even the implied warranty of
##	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
##	See the GNU Lesser General Public License for more details.
##
##	You should have received a copy of the GNU Lesser General Public
##	License along with qpOASES; if not, write to the Free Software
##	Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
##



##
##	Filename:  make_linux.mk
##	Author:    Hans Joachim Ferreau, Andreas Potschka, Christian Kirches
##	Version:   3.2
##	Date:      2007-2017
##

################################################################################
# user configuration

# include directories, relative
IDIR =   ${TOP}/include
SRCDIR = ${TOP}/src
BINDIR = ${TOP}/bin
DEBUGDIR = ${TOP}/debug

# Include directories (ADAPT TO YOUR LOCAL SETTINGS!)
# 1) qpOASES
QPOASES_IDIR   = /usr/local/qpOASES/include
QPOASES_LIB_DIR = /usr/local/qpOASES/bin
QPOASES_LINK   = -L${QPOASES_LIB_DIR} -Wl,-rpath=${QPOASES_LIB_DIR} -lqpOASES 

# 2) CasADi
CASADI_IDIR   = /usr/local/casadi/build
CASADI_LIBDIR = /usr/local/casadi/build/lib

################################################################################
# do not touch this

CPP = g++
CC  = gcc
AR  = ar
RM  = rm
F77 = gfortran
ECHO = echo
CD = cd
CP = cp

# file extensions
OBJEXT = o
LIBEXT = a
DLLEXT = so
EXE =
MEXOCTEXT = mex
DEF_TARGET = -o $@
SHARED = -shared

CPPFLAGS = -Wall -pedantic -Wshadow -Wfloat-equal -O3 -Wconversion -Wsign-conversion -fPIC -DLINUX -D__USE_LONG_INTEGERS__ -D__USE_LONG_FINTS__ -D${DEF_SOLVER} -D__NO_COPYRIGHT__
#          -g -D__DEBUG__ -D__NO_COPYRIGHT__ -D__SUPPRESSANYOUTPUT__ -D__USE_SINGLE_PRECISION__

LINK_LIBRARIES = ${QPOASES_LINK}
LINK_DEPENDS = ${BINDIR}/liblcqpOASES.${LIBEXT} ${BINDIR}/liblcqpOASES.${DLLEXT}


##
##	end of file
##
