cmake_minimum_required(VERSION 2.8.8)
project("renderer")

find_package( pybind11 )
if (pybind11_FOUND)
    message("-- Found pybind11 version ${pybind11_VERSION}: ${pybind11_INCLUDE_DIRS}")
else ()
    message(FATAL_ERROR "pybind11 not found")
endif ()

# required since eigen has eigen3 intermediate directory in include
include_directories( ${CMAKE_INSTALL_PREFIX}/include/eigen3 )
include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/include )




pybind11_add_module(renderer src/main.cpp
	                         src/math/sph.cpp
							 src/util/string.cpp
							 src/util/threadpool.cpp
							 src/integrators/simplept.cpp
							 src/lights/pointlight.cpp
							 src/scene.cpp
							 src/volume.cpp
							 src/Image.cpp
							 src/PNSolution.cpp
	                         )
target_link_libraries(renderer c:/libs-msvs2017/lib/Half.lib
	                           c:/libs-msvs2017/lib/Iex-2_1.lib
							   c:/libs-msvs2017/lib/IexMath-2_1.lib
							   c:/libs-msvs2017/lib/IlmImf-2_1.lib
							   c:/libs-msvs2017/lib/IlmThread-2_1.lib
							   c:/libs-msvs2017/lib/Imath-2_1.lib )
target_link_libraries(renderer c:/libs-msvs2017/lib/zlibstatic.lib )
