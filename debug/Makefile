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
##	Filename:  examples/Makefile
##	Author:    Hans Joachim Ferreau
##	Version:   3.2
##	Date:      2007-2017
##

$(info ************  DEBUG  ************)

include ../make.mk

##
##	flags
##

IFLAGS      =  -I. \
               -I${IDIR} \
			   -I/home/jonas/Documents/casadi

LFLAGS = -L/home/jonas/Documents/casadi/lib -lcasadi

DEBUGFLAGS = -g

##
##  files
##

QPOASES_SRC = $(wildcard ${SRCDIR}/*.cpp)

WARM_UP = ${BINDIR}/warm_up${EXE}

MOVING_MASSES = ${BINDIR}/moving_masses_condensed${EXE}

LCQP_FROM_FILE = ${BINDIR}/solve_lcqp_from_file${EXE}

##
##	targets
##

all: ${MOVING_MASSES} ${WARM_UP} ${LCQP_FROM_FILE}

.PHONY: moving_masses warm_up lcqp_from_file

moving_masses: clean ${MOVING_MASSES}
warm_up: clean ${WARM_UP}
lcqp_from_file: clean ${LCQP_FROM_FILE}

${BINDIR}/%${EXE}:
	@${ECHO} "Building debug..." $@
	@${CPP} ${DEBUGFLAGS} $(notdir $@).cpp ${QPOASES_SRC} ${DEF_TARGET} ${IFLAGS} ${LFLAGS}

clean:
	@${ECHO} "Cleaning up (examples)"
	@${RM} -f *.${OBJEXT} ${WARM_UP} ${MOVING_MASSES} ${LCQP_FROM_FILE}

clobber: clean



##
##	end of file
##
