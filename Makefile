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


include make.mk

##
##	targets
##


all: src examples
#src_aw testing

src:
	@cd $@; ${MAKE} -s 

examples: src
	@cd $@; ${MAKE} -s

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

dependencies: 
	@${ECHO} "Building dependency qpOASES at ${TOP}/external/qpOASES"
	@${CD} ${TOP}/external/qpOASES; ${MKDIR} -p bin; ${MAKE} all; ${CP} bin/libqpOASES.* /usr/local/lib; ${MAKE} clean

	@${ECHO} "Building dependency OSQP at ${TOP}/external/osqp"
	@${RM} -r ${TOP}/external/osqp/build; ${MKDIR} -p ${TOP}/external/osqp/build
	@${CD} ${TOP}/external/osqp/build; ${CMAKE} -G "Unix Makefiles" -DDLONG=OFF ..; ${CMAKE} --build .; ${CMAKE} --build . --target install;

.PHONY : all src examples doc testing debugging clean clobber dependencies


##
##   end of file
##
