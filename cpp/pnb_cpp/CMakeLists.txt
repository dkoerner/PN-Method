cmake_minimum_required(VERSION 2.8.8)
project("pnb")

find_package( pybind11 )
if (pybind11_FOUND)
    message("-- Found pybind11 version ${pybind11_VERSION}: ${pybind11_INCLUDE_DIRS}")
else ()
    message(FATAL_ERROR "pybind11 not found")
endif ()

# required since eigen has eigen3 intermediate directory in include
include_directories( ${CMAKE_INSTALL_PREFIX}/include/eigen3 )

## houio
#include_directories( ${CMAKE_INSTALL_PREFIX}/include/houio )

include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/include )
#pybind11_add_module(pnb src/main.cpp src/util/string.cpp src/util/dda2d.cpp)
pybind11_add_module(pnbuilder_cpp src/main.cpp src/field/CoefficientGrid.cpp src/common/Domain.cpp src/common/GridLocation.cpp)



#target_link_libraries( pnb ${CMAKE_INSTALL_PREFIX}/lib/houio.lib )
