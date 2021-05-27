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
			   -I${QPOASES_IDIR} \
			   -I${OSQP_IDIR}

LCQPanther_DEPENDS = \
	${IDIR}/LCQProblem.hpp \
	${IDIR}/Utilities.hpp \
	${IDIR}/Subsolver.hpp \
	${IDIR}/SubsolverOSQP.hpp \
	${IDIR}/SubsolverQPOASES.hpp \
	${IDIR}/SubsolverBase.hpp

LCQPanther_OBJECTS = \
	LCQProblem.${OBJEXT} \
	Utilities.${OBJEXT} \
	Subsolver.${OBJEXT} \
	SubsolverOSQP.${OBJEXT} \
	SubsolverQPOASES.${OBJEXT}

##
##	targets
##

all: ${LCQPanther_LIBS}


${LIBDIR}/libLCQPanther.${LIBEXT}: ${LCQPanther_OBJECTS}
	@mkdir -p ${LIBDIR}
	@${ECHO} "Creating static lib" $@
	@${AR} r $@ $^

${LIBDIR}/libLCQPanther.${DLLEXT}: ${LCQPanther_OBJECTS}
	@mkdir -p ${LIBDIR}
	@${ECHO} "Creating shared lib" $@
	@${CPP} ${DEF_TARGET} ${SHARED} $^ ${LINK_LIBRARIES}

clean:
	@${ECHO} "Cleaning up (src)"
	@${RM} -f *.${OBJEXT} *.${LIBEXT} *.${DLLEXT}

clobber: clean

%.${OBJEXT}: %.cpp ${LCQPanther_DEPENDS}
	@${ECHO} "Creating" $@
	@${ECHO} ${CPP} ${DEF_TARGET} -c ${IFLAGS} ${CPPFLAGS} $<
	@${CPP} ${DEF_TARGET} -c ${IFLAGS} ${CPPFLAGS} $<

##
##	end of file
##
