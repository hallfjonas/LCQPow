cmake_minimum_required(VERSION 3.0) # version can be different

project(lcqpOASES) #name of your project

add_subdirectory(UnitTesting/googletest) # add googletest subdirectory

include_directories(/usr/lib include /usr/local/qpOASES/include) # this is so we can #include <gtest/gtest.h>

add_executable(TestUtils UnitTesting/TestUtils.cpp src/Utilities.cpp) # add this executable

target_link_libraries(TestUtils PRIVATE gtest) # link google test to this executable