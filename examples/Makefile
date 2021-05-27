##
##	This file is part of LCQPanther.
##
##	LCQPanther -- A Solver for Quadratic Programs with Commplementarity Constraints.
##	Copyright (C) 2020 - 2021 by Jonas Hall et al.
##
##	LCQPanther is free software; you can redistribute it and/or
##	modify it under the terms of the GNU Lesser General Public
##	License as published by the Free Software Foundation; either
##	version 2.1 of the License, or (at your option) any later version.
##
##	LCQPanther is distributed in the hope that it will be useful,
##	but WITHOUT ANY WARRANTY; without even the implied warranty of
##	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
##	See the GNU Lesser General Public License for more details.
##
##	You should have received a copy of the GNU Lesser General Public
##	License along with LCQPanther; if not, write to the Free Software
##	Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
##

include ../make.mk

##
##	flags
##

IFLAGS      =  -I. \
               -I${IDIR} \
			   -isystem${QPOASES_IDIR}\
			   -isystem${OSQP_IDIR}


LFLAGS 		= ${QPOASES_LINK} ${OSQP_LINK}

DEBUGFLAGS = -g -Wall -pedantic -Wshadow -Wfloat-equal

QPOASES_SRC = $(wildcard ${QPOASESDIR}/src/*.cpp)

##
##  files
##

LCQPanther_EXAMPLES = \
	${BINDIR}/warm_up${EXE} \
	${BINDIR}/warm_up_w_A${EXE} \
	${BINDIR}/warm_up_sparse${EXE} \
	${BINDIR}/solve_lcqp_from_file${EXE} \
	${BINDIR}/OptimizeOnCircle${EXE}

LCQPanther_DEBUG = \
	${DEBUGDIR}/warm_up${EXE} \
	${DEBUGDIR}/warm_up_w_A${EXE} \
	${DEBUGDIR}/warm_up_sparse${EXE} \
	${DEBUGDIR}/solve_lcqp_from_file${EXE} \
	${DEBUGDIR}/OptimizeOnCircle${EXE}

LCQPanther_SRC = $(wildcard ${SRCDIR}/*.cpp)
OSQP_SRC = $(wildcard /usr/local/osqp/src/*.c)

##
##	targets
##

all: clean ${LCQPanther_DEBUG} ${LCQPanther_EXAMPLES}

.PHONY: debug examples

debug: clean ${LCQPanther_DEBUG}

examples: clean ${LCQPanther_EXAMPLES}

moving_masses: clean ${MOVING_MASSES}
warm_up: clean ${WARM_UP}
lcqp_from_file: clean ${LCQP_FROM_FILE}

${BINDIR}/%${EXE}: %.${OBJEXT} ${LCQPanther_LIBS}
	@mkdir -p ${BINDIR}
	@${ECHO} "Creating" $@
	@${CPP} ${DEF_TARGET} ${CPPFLAGS} $< ${LCQPanther_LINK} ${LFLAGS}

# If we need to debug and step into qpOASES: Replace LFLAGS with QPOASES_SRC and move to front
${DEBUGDIR}/%${EXE}:
	@mkdir -p ${DEBUGDIR}
	@${ECHO} "Building debug..." $@
	${CPP} ${DEBUGFLAGS} ${PLOTLIBFLAGS} $(notdir $@).cpp ${LCQPanther_SRC} ${DEF_TARGET} ${IFLAGS} ${LFLAGS}
	@${CPP} ${DEBUGFLAGS} ${PLOTLIBFLAGS} $(notdir $@).cpp ${LCQPanther_SRC} ${DEF_TARGET} ${IFLAGS} ${LFLAGS}

clean:
	@${ECHO} "Cleaning up (examples)"
	@${RM} -f *.${OBJEXT} ${LCQPanther_EXAMPLES} ${LCQPanther_DEBUG}

clobber: clean

${LCQPanther_LIBS}:
	@cd ..; ${MAKE} -s src

%.${OBJEXT}: %.cpp
	@${ECHO} "Creating" $@
	@${CPP} ${DEF_TARGET} -c ${IFLAGS} ${CPPFLAGS} $<


##
##	end of file
##
