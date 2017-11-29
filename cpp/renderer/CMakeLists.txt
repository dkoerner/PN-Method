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
include_directories( ${CMAKE_INSTALL_PREFIX}/include/OpenEXR )
include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/include )


# temp: mitsubacore
include_directories( c:/libs/mitsuba/include )
#ADD_DEFINITIONS( -DSINGLE_PRECISION )
ADD_DEFINITIONS( -DDOUBLE_PRECISION )
ADD_DEFINITIONS( -DSPECTRUM_SAMPLES=3 )
ADD_DEFINITIONS( -DMTS_BUILD_MODULE )
ADD_DEFINITIONS( -DMTS_HAS_OPENEXR )




pybind11_add_module(renderer src/main.cpp
	                         src/math/sph.cpp
							 src/math/shvector.cpp
							 src/util/string.cpp
							 src/util/threadpool.cpp
							 src/util/voxelgrid.cpp
							 src/integrators/simplept.cpp
							 src/integrators/pnispt.cpp
							 src/integrators/jispt.cpp
							 src/integrators/directpn.cpp
							 src/lights/pointlight.cpp
							 src/lights/directionallight.cpp
							 src/fields/VoxelGridField.cpp
							 src/scene.cpp
							 src/volume.cpp
							 src/Integrator.cpp
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
target_link_libraries(renderer ws2_32 psapi )
target_link_libraries(renderer c:/libs-msvs2017/lib/libboost_filesystem-vc141-mt-1_65_1.lib
	                           c:/libs-msvs2017/lib/libboost_system-vc141-mt-1_65_1.lib
							   c:/libs-msvs2017/lib/libboost_thread-vc141-mt-1_65_1.lib)

# temp: mitsubacore
target_link_libraries(renderer c:/libs/mitsuba-buildmvs17/Release/libcore.lib )

# temp: houio for bgeo loading and writing
target_link_libraries(renderer c:/libs-msvs2017/lib/houio.lib )
