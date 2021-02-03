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
##	Filename:  src/Makefile
##	Author:    Hans Joachim Ferreau, Andreas Potschka, Christian Kirches
##	Version:   3.2
##	Date:      2007-2017
##

include ../make.mk

##
##	flags
##

IFLAGS      =  -I. \
               -I${IDIR}

LCQPOASES_OBJECTS = \
	LCQProblem.${OBJEXT} \
	Utilities.${OBJEXT}

QPOASES_DEPENDS = \
	${IDIR}/LCQProblem.hpp \
	${IDIR}/Utilities.hpp


##
##	targets
##

all: ${LINK_DEPENDS}


${BINDIR}/libqpOASES.${LIBEXT}: ${QPOASES_OBJECTS} ${QPOASES_EXTRAS_OBJECTS}
	@${ECHO} "Creating static lib" $@
	@${AR} r $@ $^

${BINDIR}/libqpOASES.${DLLEXT}: ${QPOASES_OBJECTS} ${QPOASES_EXTRAS_OBJECTS}
	@${ECHO} "Creating shared lib" $@
	@${CPP} ${DEF_TARGET} ${SHARED} $^ ${LINK_LIBRARIES}

clean:
	@${ECHO} "Cleaning up (src)"
	@${RM} -f *.${OBJEXT} *.${LIBEXT} *.${DLLEXT}

clobber: clean

%.${OBJEXT}: %.cpp ${QPOASES_DEPENDS}
	@${ECHO} "Creating" $@
	@${CPP} ${DEF_TARGET} -c ${IFLAGS} ${CPPFLAGS} $<

##
##	end of file
##
