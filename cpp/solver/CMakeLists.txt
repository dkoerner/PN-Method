cmake_minimum_required(VERSION 2.8.8)
project("solver")

find_package( pybind11 )
if (pybind11_FOUND)
    message("-- Found pybind11 version ${pybind11_VERSION}: ${pybind11_INCLUDE_DIRS}")
else ()
    message(FATAL_ERROR "pybind11 not found")
endif ()

# required since eigen has eigen3 intermediate directory in include
include_directories( ${CMAKE_INSTALL_PREFIX}/include/eigen3 )

include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/include )
pybind11_add_module(solver src/main.cpp src/util/string.cpp src/util/dda2d.cpp)
