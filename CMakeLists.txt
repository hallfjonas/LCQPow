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

# cmake minimum required version
cmake_minimum_required(VERSION 3.10.2)

# Project name
project(lcqpanther)

# Save auxiliar source files to variable
aux_source_directory(src SRC_FILES)

## Options ------------------------------------------------------------------------------
option(
    BUILD_EXAMPLES
    "Option to build examples."
    ON
)

option(
    BUILD_DEBUG
    "Option to compile with the -g flag."
    OFF
)

option(
    RUN_GPROF_ANALYSIS
    "Option to compile with the -pg flag and run an analysis on all examples."
    OFF
)

option(
    RUN_VALGRIND_ANALYSIS
    "Option to compile with the -g flag and run a valgrind analysis on all examples."
    OFF
)

## DOWNLOAD AND CHECK ALL SUBMODULES ---- from Code, Tech, and Tutorials ----------------
find_package(Git QUIET)
if (GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
# Update submodules as needed
    option(GIT_SUBMODULE "Check submodules during build" ON)
    if (GIT_SUBMODULE)
        message(STATUS "Submodule update")
        execute_process(
            COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
            RESULT_VARIABLE GIT_SUBMOD_RESULT
        )

        if (NOT GIT_SUBMOD_RESULT EQUAL "0")
            message(
                FATAL_ERROR
                "git submodule update --init --recursive failed with ${GIT_SUBMOD_RESULT}."
            )
        endif()
    endif()
endif()

if (NOT EXISTS "${PROJECT_SOURCE_DIR}/external/osqp/CMakeLists.txt")
    message(
        FATAL_ERROR
        "The submodule OSQP does not exist. GIT_SUBMODULE was turned off or failed."
    )
endif()

if (NOT EXISTS "${PROJECT_SOURCE_DIR}/external/qpOASES/CMakeLists.txt")
    message(
        FATAL_ERROR
        "The submodule qpOASES does not exist. GIT_SUBMODULE was turned off or failed."
    )
endif()


if (NOT EXISTS "${PROJECT_SOURCE_DIR}/external/googletest/CMakeLists.txt")
    message(
        FATAL_ERROR
        "The submodule googletest does not exist. GIT_SUBMODULE was turned off or failed."
    )
endif()

## ADD ALL SUBMODULES -------------------------------------------------------------------
# OSQP
add_subdirectory(external/osqp)

# qpOASES
add_subdirectory(external/qpOASES)

# qpOASES
add_subdirectory(external/googletest)

# Add include directories of dependencies: qpOASES, OSQP
include_directories(${PROJECT_NAME}
    PRIVATE external/osqp/include
    PRIVATE external/qpOASES/include
)

## CREATE LCQPanther libraries ----------------------------------------------------------
# Include project headers
include_directories(${PROJECT_NAME} PUBLIC include)

# Add compile options
add_compile_options(
    -Wall
    -pedantic
    -Wshadow
    -Wfloat-equal
    -O3
    -Wconversion
    -Wsign-conversion
    -fPIC
    -DLINUX
)

# create static lib
add_library(${PROJECT_NAME}-static STATIC ${SRC_FILES})
set_target_properties(${PROJECT_NAME}-static PROPERTIES OUTPUT_NAME ${PROJECT_NAME})

# create shared lib
add_library(${PROJECT_NAME}-shared SHARED ${SRC_FILES})
set_target_properties(${PROJECT_NAME}-shared PROPERTIES OUTPUT_NAME ${PROJECT_NAME})

## Build debug --------------------------------------------------------------------------
if (${BUILD_DEBUG} OR ${RUN_VALGRIND_ANALYSIS} OR ${RUN_GPROF_ANALYSIS})
    add_compile_options(-g)
endif()

## Build gprof --------------------------------------------------------------------------
if (${RUN_GPROF_ANALYSIS})
    add_compile_options(-pg)
endif()

## Build examples -----------------------------------------------------------------------
if (${BUILD_EXAMPLES})
    aux_source_directory(examples EXAMPLE_FILES)

    FOREACH(ELEMENT ${EXAMPLE_FILES})
        get_filename_component(EXAMPLE_NAME ${ELEMENT} NAME_WE)
        add_executable(${EXAMPLE_NAME} ${ELEMENT})
        target_link_libraries(
            ${EXAMPLE_NAME}
            PUBLIC ${PROJECT_NAME}-shared
            PRIVATE qpOASES osqp gtest
        )
    endforeach()
endif()

## Build unit tests ---------------------------------------------------------------------
# Add Unit testing source file
add_executable(RunUnitTests test/RunUnitTests.cpp)

# Add libraries: gtest, qpOASES, OSQP
target_link_libraries(
    RunUnitTests
    PUBLIC ${PROJECT_NAME}-shared
    PRIVATE qpOASES osqp gtest
)

## End of cmake file --------------------------------------------------------------------
