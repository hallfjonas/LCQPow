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

include ../../make.mk

##
##	flags
##

IDIR   = ../../include

IFLAGS      =  -I. \
               -I${IDIR} \
               -I${MATLAB_IDIR}
			   
LINK_MATLAB_LIBS = -L${MATLAB_LIBDIR} -lmex -lmat -lmx

EXTRAFLAGS = -largeArrayDims -O -D__cpluplus -D__MATLAB__ -fexceptions

LCQPOASES_OBJECT = LCQPanther.${OBJEXT}

##
##	targets
##

all: ${BINDIR}/LCQPanther.${MEXEXT} copyMs

${BINDIR}/qpOASES.${MEXEXT}: ${LCQPOASES_OBJECT} 
	@${ECHO} "Creating" $@
	@${CPP} -shared -o $@ $? ${QPOASES_LINK} ${LINK_LIBRARIES} ${LINK_MATLAB_LIBS}

copyMs:
	@${ECHO} "Copying M files ..."
	@${CP} LCQPanther.m ${BINDIR}
	
clean:
	@${RM} -f *.${OBJEXT} *.mex*

clobber: clean

%.${OBJEXT}: %.cpp 
	@${ECHO} "Creating" $@
	@${CPP} ${DEF_TARGET} -c ${IFLAGS} ${EXTRAFLAGS} ${CPPFLAGS} $<
