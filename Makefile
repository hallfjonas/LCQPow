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


include make.mk

##
##	targets
##


all: src examples
#src_aw testing

src:
	@cd $@; ${MAKE} -s

examples: src
	@cd $@; ${MAKE} examples -s

doc:
	@cd $@; ${MAKE} -s

test:
	@mkdir -p ${BUILDDIR}
	@cd build && ${CMAKE} ..
	@cd build && ${MAKE}
	build/TestUtils

debugging:
	@cd $@; ${MAKE} -s

clean:
	@cd src               			&& ${MAKE} -s clean
	@cd examples          			&& ${MAKE} -s clean
	@${ECHO} "Cleaning up (debug)"  && ${RM} -rf debug
	@${ECHO} "Cleaning up (bin)"  	&& ${RM} -rf bin
	@${ECHO} "Cleaning up (build)"  && ${RM} -rf build
	@${ECHO} "Cleaning up (lib)"  	&& ${RM} -rf lib

clobber: clean

install:
	@${CP} ${LIBDIR}/* ${INSTALLDIR}

.PHONY : all src examples doc testing debugging clean clobber dependencies


##
##   end of file
##
