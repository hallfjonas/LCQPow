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

# qpOASES with or without Schur Complement method
option(
    QPOASES_SCHUR_COMPLEMENT
    "Whether to enable the qpOASES Schur Complement method (sparse solver must be available, e.g. through MATLAB)"
    OFF
)
if (${QPOASES_SCHUR_COMPLEMENT})
    add_compile_options(-DSOLVER_MA57)
    message("WARNING: This method requires the availability of MA57.")
endif()

set(
    qpOASES_library_link
    "-L/usr/local/lib -lqpOASES"
)

set(
    osqp_library_link
    "-L/usr/local/lib -losqp"
)

set(
    CMAKE_CXX_STANDARD_INCLUDE_DIRECTORIES
    "/usr/local/include/qpOASES"
    "/usr/local/include/osqp"
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
    -D__USE_LONG_INTEGERS__
    -D__USE_LONG_FINTS__
    -D__AVOID_LA_NAMING_CONFLICTS__
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

# Add include directories of dependencies: qpOASES, OSQP
include_directories(${PROJECT_NAME}
    PRIVATE ${qpOASES_include}
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
            PRIVATE ${qpOASES_library_link} ${osqp_library_link}
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
#         PRIVATE ${qpOASES_include}
#         PRIVATE ${osqp_include}
#     )

#     # Add Unit testing source file
#     add_executable(RunUnitTests test/RunUnitTests.cpp)

#     add_dependencies(
#         RunUnitTests
#         qpOASES
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
#         PRIVATE ${qpOASES_lib} ${osqp_lib} ${gtest_lib}
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
