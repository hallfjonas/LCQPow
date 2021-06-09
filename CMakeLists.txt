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
cmake_minimum_required(VERSION 3.20.3)

# Project name
project(lcqpanther)

## Options ------------------------------------------------------------------------------
option(
    BUILD_EXAMPLES
    "Option to build examples"
    ON
)

option(
    BUILD_DOCUMENTATION
    "Option to build documentation"
    ON
)

option(
    UNIT_TESTS
    "Option to build and run unit tests"
    ON
)

option(
    __USE_LONG_INTEGERS__
    "Use long integers"
    ON
)

option(
    __USE_LONG_FINTS__
    "Use long fortran ints"
    ON
)

## Paths --------------------------------------------------------------------------------
set(
    CMAKE_INSTALL_PREFIX
    "${CMAKE_BINARY_DIR}"
)

## ADD ALL EXTERNAL PROJECTS ------------------------------------------------------------

# external modules
include(ExternalProject)

# qpOASES
ExternalProject_Add(
    qpoases
    GIT_REPOSITORY "https://github.com/hallfjonas/qpOASES.git"
    PREFIX external
    CONFIGURE_COMMAND ""
    BUILD_IN_SOURCE ON
    BUILD_COMMAND make src;
    INSTALL_DIR ${CMAKE_INSTALL_PREFIX}
    INSTALL_COMMAND cp ${CMAKE_BINARY_DIR}/external/src/qpoases/bin/libqpOASES.so ${CMAKE_INSTALL_PREFIX}/lib/libqpOASES.so
)

set(qpoases_lib "-L${CMAKE_BINARY_DIR}/lib -lqpOASES")
set(qpoases_include "${CMAKE_BINARY_DIR}/external/src/qpoases/include")

# OSQP
ExternalProject_Add(
    osqp
    GIT_REPOSITORY "https://github.com/hallfjonas/osqp.git"
    PREFIX external
    CMAKE_ARGS
        -DDLONG=OFF
        -DBUILD_SHARED_LIBS=ON
    BUILD_COMMAND cmake --build .
    INSTALL_COMMAND cp ${CMAKE_BINARY_DIR}/external/src/osqp-build/out/libosqp.so ${CMAKE_INSTALL_PREFIX}/lib/libosqp.so
)

set(osqp_lib "-L${CMAKE_BINARY_DIR}/lib -losqp")
set(osqp_include "${CMAKE_BINARY_DIR}/external/src/osqp/include")

# googletest
ExternalProject_Add(
    gtest
    GIT_REPOSITORY "https://github.com/google/googletest.git"
    PREFIX external
    INSTALL_COMMAND cp ${CMAKE_BINARY_DIR}/external/src/gtest-build/lib/libgtest.a ${CMAKE_INSTALL_PREFIX}/lib/
)

set(gtest_lib "-L${CMAKE_BINARY_DIR}/external/src/gtest-build/out -lgtest")
set(gtest_include "${CMAKE_BINARY_DIR}/external/src/gtest/googletest/include")

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
    -D__USE_LONG_INTEGERS__
    -D__USE_LONG_FINTS__
    -DSOLVER_MA57
)

# Save auxiliar source files to variable
aux_source_directory(src SRC_FILES)

# create static lib
add_library(${PROJECT_NAME}-static STATIC ${SRC_FILES})
set_target_properties(
    ${PROJECT_NAME}-static
    PROPERTIES
    OUTPUT_NAME ${PROJECT_NAME}
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
)

# Add the external projects as depenedncies for our project
add_dependencies(
    ${PROJECT_NAME}-static
    qpoases
    osqp
)

# Add include directories of dependencies: qpOASES, OSQP
include_directories(${PROJECT_NAME}
    PRIVATE ${qpoases_include}
    PRIVATE ${osqp_include}
)

# create shared lib
add_library(${PROJECT_NAME}-shared SHARED ${SRC_FILES})
set_target_properties(
    ${PROJECT_NAME}-shared
    PROPERTIES
    OUTPUT_NAME ${PROJECT_NAME}
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
)

# Add the external projects as depenedncies for our project
add_dependencies(
    ${PROJECT_NAME}-shared
    qpoases
    osqp
)

## Build examples -----------------------------------------------------------------------
if (${BUILD_EXAMPLES})
    aux_source_directory(examples EXAMPLE_FILES)

    FOREACH(ELEMENT ${EXAMPLE_FILES})
        # get filename w/o dir and extension
        get_filename_component(EXAMPLE_NAME ${ELEMENT} NAME_WE)

        # generate executable target
        add_executable(${EXAMPLE_NAME} ${ELEMENT})

        # link libraries
        target_link_libraries(
            ${EXAMPLE_NAME}
            PUBLIC ${PROJECT_NAME}-shared
            PRIVATE ${qpoases_lib} ${osqp_lib}
        )

        # specify output directory
        set_target_properties(
            ${EXAMPLE_NAME}
            PROPERTIES
            RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin/examples"
        )
    endforeach()
endif()

# ## Build unit tests ---------------------------------------------------------------------
# if (${UNIT_TESTS})
#     # Add include directories of dependencies: qpOASES, OSQP
#     include_directories(${PROJECT_NAME}
#         PRIVATE ${gtest_include}
#         PRIVATE ${qpoases_include}
#         PRIVATE ${osqp_include}
#     )

#     # Add Unit testing source file
#     add_executable(RunUnitTests test/RunUnitTests.cpp)

#     add_dependencies(
#         RunUnitTests
#         qpoases
#         osqp
#         gtest
#     )

#     set_target_properties(
#         RunUnitTests
#         PROPERTIES
#         RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin/tests"
#     )

#     # Add libraries: gtest, qpOASES, OSQP
#     target_link_libraries(
#         RunUnitTests
#         PUBLIC ${PROJECT_NAME}-shared
#         PRIVATE ${qpoases_lib} ${osqp_lib} ${gtest_lib}
#     )
# endif()

# ## Build documentation
# if (${BUILD_DOCUMENTATION})
#     find_package(
#         Doxygen
#         REQUIRED dot
#     )

#     if (DOXYGEN_FOUND)
#         message(
#             STATUS
#             "Found doxygen: Building documentation."
#         )

#         doxygen_add_docs(
#             doxy
#             "include"
#         )

#         set(DOXYGEN_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/../docs")

#         set(DOXYGEN_GENERATE_HTML YES)
#         set(DOXYGEN_GENERATE_MAN YES)

#         set(DOXYGEN_PROJECT_BRIEF "LCQPanther -- A Linear Complementarity Quadratic Programming Solver")
#     else()
#         message(
#             WARNING
#             "Could not find doxygen: Not building documentation."
#         )
#     endif()
# endif()

## End of cmake file --------------------------------------------------------------------
