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

include ../make.mk

##
##	flags
##

IFLAGS      =  -I. \
               -I${IDIR} \
			   -I${QPOASES_IDIR}

LFLAGS 		= ${QPOASES_LINK}

DEBUGFLAGS = -g -D__USE_LONG_INTEGERS__ -Wall -pedantic -Wshadow -Wfloat-equal

QPOASES_SRC = $(wildcard ${SRCDIR}/qpOASES/*.cpp)

##
##  files
##

LCQPOASES_SRC = $(wildcard ${SRCDIR}/*.cpp)
QPOASES_EXM = $(wildcard ./*.cpp)


QPOASES_EXES = \
	${BINDIR}/warm_up${EXE} \
	${BINDIR}/warm_up_from_file${EXE} \
	${BINDIR}/initial_value_OCP${EXE} \
	${BINDIR}/michael_schumacher_problem${EXE} \
	${BINDIR}/solve_lcqp_from_file${EXE}
	
QPOASES_DEBUG = $(wildcard ${DEBUGDIR}/*${EXE})

WARM_UP = ${DEBUGDIR}/warm_up${EXE}

MOVING_MASSES = ${DEBUGDIR}/moving_masses_condensed${EXE}

LCQP_FROM_FILE = ${DEBUGDIR}/solve_lcqp_from_file${EXE}

PLOTS = ${BINDIR}/minimal_matplotlib

##
##	targets
##

all: ${QPOASES_DEBUG} ${QPOASES_EXES}

.PHONY: moving_masses warm_up lcqp_from_file

moving_masses: clean ${MOVING_MASSES}
warm_up: clean ${WARM_UP}
lcqp_from_file: clean ${LCQP_FROM_FILE}

${BINDIR}/%${EXE}: %.${OBJEXT} ${LINK_DEPENDS}
	@mkdir -p ${BINDIR}
	@${ECHO} "Creating" $@
	${CPP} ${DEF_TARGET} ${CPPFLAGS} $< ${QPOASES_LINK} ${LINK_LIBRARIES}
	@${CPP} ${DEF_TARGET} ${CPPFLAGS} $< ${QPOASES_LINK} ${LINK_LIBRARIES}

# If we need to debug and step into qpOASES: Replace LFLAGS with QPOASES_SRC and move to front
${DEBUGDIR}/%${EXE}:
	@mkdir -p ${DEBUGDIR}
	@${ECHO} "Building debug..." $@
	${CPP} ${DEBUGFLAGS} $(notdir $@).cpp ${LCQPOASES_SRC} ${PLOT_SRC} ${DEF_TARGET} ${IFLAGS} ${LFLAGS}
	@${CPP} ${DEBUGFLAGS} $(notdir $@).cpp ${LCQPOASES_SRC} ${PLOT_SRC} ${DEF_TARGET} ${IFLAGS} ${LFLAGS}

clean:
	@${ECHO} "Cleaning up (examples)"
	@${RM} -f *.${OBJEXT} ${QPOASES_EXES} ${QPOASES_DEBUG}

clobber: clean

${LINK_DEPENDS}:
	@cd ..; ${MAKE} -s src

%.${OBJEXT}: %.cpp
	@${ECHO} "Creating" $@
	${CPP} ${DEF_TARGET} -c ${IFLAGS}/usr/local/qpOASES/${CPPFLAGS} $<
	@${CPP} ${DEF_TARGET} -c ${IFLAGS} ${CPPFLAGS} $<


##
##	end of file
##
