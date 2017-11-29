cmake_minimum_required(VERSION 2.8.8)
project("pns")

find_package( pybind11 )
if (pybind11_FOUND)
    message("-- Found pybind11 version ${pybind11_VERSION}: ${pybind11_INCLUDE_DIRS}")
else ()
    message(FATAL_ERROR "pybind11 not found")
endif ()

# required since eigen has eigen3 intermediate directory in include
include_directories( ${CMAKE_INSTALL_PREFIX}/include/eigen3 )
include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/include )
add_definitions(-DEIGEN_INITIALIZE_MATRICES_BY_NAN)


pybind11_add_module(pnsolver src/main.cpp
	                         src/PNSystem.cpp
							 src/PNVolume.cpp
							 src/PNSolution.cpp
							 src/MultigridSolver.cpp
							 src/stencil_fopn_p0_cg.cpp
							 src/stencil_fopn_p0_sg.cpp
							 src/stencil_fopn_p1_cg.cpp
							 src/stencil_fopn_p1_sg.cpp
							 src/stencil_fopn_p2_cg.cpp
							 src/stencil_fopn_p2_sg.cpp
							 src/stencil_fopn_p3_cg.cpp
							 src/stencil_fopn_p3_sg.cpp
							 src/stencil_fopn_p4_cg.cpp
							 src/stencil_fopn_p4_sg.cpp
							 src/stencil_fopn_p5_cg.cpp
							 src/stencil_fopn_p5_sg.cpp
							 src/stencil_sopn_p0_cg.cpp
							 src/stencil_sopn_p0_sg.cpp
							 src/stencil_sopn_p1_cg.cpp
							 src/stencil_sopn_p1_sg.cpp
							 src/stencil_sopn_p2_cg.cpp
							 src/stencil_sopn_p2_sg.cpp
							 src/stencil_sopn_p3_cg.cpp
							 src/stencil_sopn_p3_sg.cpp
							 src/stencil_sopn_p4_cg.cpp
							 src/stencil_sopn_p4_sg.cpp
							 src/stencil_sopn_p5_cg.cpp
							 src/stencil_sopn_p5_sg.cpp
							 src/stencil_cda.cpp
							 src/fields/VoxelGridField.cpp
							 src/common/Domain.cpp
							 src/util/threadpool.cpp
							 src/util/voxelgrid.cpp
							 src/util/mem.cpp
							 src/math/sph.cpp
							 src/solver.cpp)
