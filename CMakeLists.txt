cmake_minimum_required(VERSION 3.0) # version can be different

project(lcqpOASES) #name of your project

add_subdirectory(UnitTesting/googletest) # add googletest subdirectory

add_compile_options(-Wall -pedantic -Wshadow -Wfloat-equal -O3 -Wconversion -Wsign-conversion -fPIC -DLINUX -D__USE_LONG_INTEGERS__ -D__USE_LONG_FINTS__)

include_directories(/usr/lib include /usr/local/qpOASES/include) # this is so we can #include <gtest/gtest.h>

add_executable(TestUtils UnitTesting/TestUtils.cpp src/Utilities.cpp src/LCQProblem.cpp) # add this executable

target_link_libraries(TestUtils PRIVATE gtest /usr/local/qpOASES/bin/libqpOASES.a /usr/local/qpOASES/bin/libqpOASES.so) # link google test to this executable