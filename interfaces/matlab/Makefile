##
##	This file is part of LCQPow.
##
##	LCQPow -- A Solver for Quadratic Programs with Commplementarity Constraints.
##	Copyright (C) 2020 - 2021 by Jonas Hall et al.
##
##	LCQPow is free software; you can redistribute it and/or
##	modify it under the terms of the GNU Lesser General Public
##	License as published by the Free Software Foundation; either
##	version 2.1 of the License, or (at your option) any later version.
##
##	LCQPow is distributed in the hope that it will be useful,
##	but WITHOUT ANY WARRANTY; without even the implied warranty of
##	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
##	See the GNU Lesser General Public License for more details.
##
##	You should have received a copy of the GNU Lesser General Public
##	License along with LCQPow; if not, write to the Free Software
##	Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
##

include ../../make.mk

##
##	flags
##

IDIR   = ../../include

IFLAGS      =  -I. \
               -I${IDIR} \
               -I${MATLAB_IDIR} \
               -I${QPOASES_IDIR} \
               -I${OSQP_IDIR}

LINK_MATLAB_LIBS = -L${MATLAB_LIBDIR} ${QPOASES_LINK} ${OSQP_LINK} -lmex -lmat -lmx

EXTRAFLAGS = -largeArrayDims -O -D__cpluplus -D__MATLAB__ -fexceptions

LCQPow_OBJECT = LCQPow.${OBJEXT}

##
##	targets
##

all: ${BINDIR}/LCQPow.${MEXEXT} copyMs

${BINDIR}/LCQPow.${MEXEXT}: ${LCQPow_OBJECT}
	@${ECHO} "Creating" $@
	@${CPP} -shared -o $@ $? ${LCQPow_LINK} ${LINK_LIBRARIES} ${LINK_MATLAB_LIBS}

copyMs:
	@${ECHO} "Copying M files ..."
	@${CP} LCQPow.m ${BINDIR}

clean:
	@${RM} -f *.${OBJEXT} *.mex*

clobber: clean

%.${OBJEXT}: %.cpp
	@${ECHO} "Creating" $@
	@${CPP} ${DEF_TARGET} -c ${IFLAGS} ${EXTRAFLAGS} ${CPPFLAGS} $<
