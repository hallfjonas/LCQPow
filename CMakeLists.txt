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

cmake_minimum_required(VERSION 3.0) # version can be different

project(LCQPanther) #name of your project

add_subdirectory(UnitTesting/googletest) # add googletest subdirectory

AUX_SOURCE_DIRECTORY(src SRC_FILES)

add_compile_options(-Wall -pedantic -Wshadow -Wfloat-equal -O3 -Wconversion -Wsign-conversion -fPIC -DLINUX -D__USE_LONG_INTEGERS__ -D__USE_LONG_FINTS__)

include_directories(SYSTEM /usr/local/include/osqp include /usr/local/qpOASES/include /usr/local/matplotlib-cpp /usr/include/python3.6m) # this is so we can #include <gtest/gtest.h>

add_executable(TestUtils UnitTesting/TestUtils.cpp ${SRC_FILES}) # add this executable

target_link_libraries(TestUtils PRIVATE gtest qpOASES osqp python3.6m) # link google test to this executable