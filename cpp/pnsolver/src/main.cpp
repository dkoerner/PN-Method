#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <complex>
#include <iostream>

#include <math/common.h>
#include <math/vector.h>
#include <math/ray.h>

#include <util/timer.h>
#include <util/mem.h>

#include<common/Domain.h>
#include<field/VoxelGridField.h>
#include<field/Constant.h>
#include<field/SHEXP.h>
#include <PNSystem.h>

#include <mgtest.h>


namespace py = pybind11;


Eigen::VectorXd to_vector( const std::vector<double>& values )
{
	Eigen::RowVectorXd x = Eigen::VectorXd( values.size() );
	for( int i=0;i<values.size();++i )
		x(i) = values[i];
	return x;
}


void buildUpAndDownsamplingMatrices( PNSystem& sys_fine, PNSystem::RealMatrix& downsampleMatrix, PNSystem::RealMatrix& upsampleMatrix )
{
	PNSystem::VoxelManager& vm_fine = sys_fine.getVoxelManager();
	PNSystem::VoxelManager vm_coarse = vm_fine.downsample();


	PNSystem::MatrixBuilderd downsampleMatrixBuilder;
	PNSystem::MatrixBuilderd upsampleMatrixBuilder;

	// 2d stamps
	std::vector<std::vector<std::pair<V3i, double>>> stamps(8);

	if( vm_fine.is2D() )
	{
		stamps[0].push_back(std::make_pair(V3i(-1, -1, 0), 0.5*0.5*1.0));
		stamps[0].push_back(std::make_pair(V3i(-1, 0, 0), 0.5*1.0*1.0));
		stamps[0].push_back(std::make_pair(V3i(-1, 1, 0), 0.5*0.5*1.0));
		stamps[0].push_back(std::make_pair(V3i(0, -1, 0), 1.0*0.5*1.0));
		stamps[0].push_back(std::make_pair(V3i(0, 0, 0), 1.0*1.0*1.0));
		stamps[0].push_back(std::make_pair(V3i(0, 1, 0), 1.0*0.5*1.0));
		stamps[0].push_back(std::make_pair(V3i(1, -1, 0), 0.5*0.5*1.0));
		stamps[0].push_back(std::make_pair(V3i(1, 0, 0), 0.5*1.0*1.0));
		stamps[0].push_back(std::make_pair(V3i(1, 1, 0), 0.5*0.5*1.0));
		stamps[1].push_back(std::make_pair(V3i(-1, -1, 0), 0.25*0.5*1.0));
		stamps[1].push_back(std::make_pair(V3i(-1, 0, 0), 0.25*1.0*1.0));
		stamps[1].push_back(std::make_pair(V3i(-1, 1, 0), 0.25*0.5*1.0));
		stamps[1].push_back(std::make_pair(V3i(0, -1, 0), 0.75*0.5*1.0));
		stamps[1].push_back(std::make_pair(V3i(0, 0, 0), 0.75*1.0*1.0));
		stamps[1].push_back(std::make_pair(V3i(0, 1, 0), 0.75*0.5*1.0));
		stamps[1].push_back(std::make_pair(V3i(1, -1, 0), 0.75*0.5*1.0));
		stamps[1].push_back(std::make_pair(V3i(1, 0, 0), 0.75*1.0*1.0));
		stamps[1].push_back(std::make_pair(V3i(1, 1, 0), 0.75*0.5*1.0));
		stamps[1].push_back(std::make_pair(V3i(2, -1, 0), 0.25*0.5*1.0));
		stamps[1].push_back(std::make_pair(V3i(2, 0, 0), 0.25*1.0*1.0));
		stamps[1].push_back(std::make_pair(V3i(2, 1, 0), 0.25*0.5*1.0));
		stamps[2].push_back(std::make_pair(V3i(-1, -1, 0), 0.25*0.25*1.0));
		stamps[2].push_back(std::make_pair(V3i(-1, 0, 0), 0.25*0.75*1.0));
		stamps[2].push_back(std::make_pair(V3i(-1, 1, 0), 0.25*0.75*1.0));
		stamps[2].push_back(std::make_pair(V3i(-1, 2, 0), 0.25*0.25*1.0));
		stamps[2].push_back(std::make_pair(V3i(0, -1, 0), 0.75*0.25*1.0));
		stamps[2].push_back(std::make_pair(V3i(0, 0, 0), 0.75*0.75*1.0));
		stamps[2].push_back(std::make_pair(V3i(0, 1, 0), 0.75*0.75*1.0));
		stamps[2].push_back(std::make_pair(V3i(0, 2, 0), 0.75*0.25*1.0));
		stamps[2].push_back(std::make_pair(V3i(1, -1, 0), 0.75*0.25*1.0));
		stamps[2].push_back(std::make_pair(V3i(1, 0, 0), 0.75*0.75*1.0));
		stamps[2].push_back(std::make_pair(V3i(1, 1, 0), 0.75*0.75*1.0));
		stamps[2].push_back(std::make_pair(V3i(1, 2, 0), 0.75*0.25*1.0));
		stamps[2].push_back(std::make_pair(V3i(2, -1, 0), 0.25*0.25*1.0));
		stamps[2].push_back(std::make_pair(V3i(2, 0, 0), 0.25*0.75*1.0));
		stamps[2].push_back(std::make_pair(V3i(2, 1, 0), 0.25*0.75*1.0));
		stamps[2].push_back(std::make_pair(V3i(2, 2, 0), 0.25*0.25*1.0));
		stamps[3].push_back(std::make_pair(V3i(-1, -1, 0), 0.5*0.25*1.0));
		stamps[3].push_back(std::make_pair(V3i(-1, 0, 0), 0.5*0.75*1.0));
		stamps[3].push_back(std::make_pair(V3i(-1, 1, 0), 0.5*0.75*1.0));
		stamps[3].push_back(std::make_pair(V3i(-1, 2, 0), 0.5*0.25*1.0));
		stamps[3].push_back(std::make_pair(V3i(0, -1, 0), 1.0*0.25*1.0));
		stamps[3].push_back(std::make_pair(V3i(0, 0, 0), 1.0*0.75*1.0));
		stamps[3].push_back(std::make_pair(V3i(0, 1, 0), 1.0*0.75*1.0));
		stamps[3].push_back(std::make_pair(V3i(0, 2, 0), 1.0*0.25*1.0));
		stamps[3].push_back(std::make_pair(V3i(1, -1, 0), 0.5*0.25*1.0));
		stamps[3].push_back(std::make_pair(V3i(1, 0, 0), 0.5*0.75*1.0));
		stamps[3].push_back(std::make_pair(V3i(1, 1, 0), 0.5*0.75*1.0));
		stamps[3].push_back(std::make_pair(V3i(1, 2, 0), 0.5*0.25*1.0));
	}else
	{
		stamps[0].push_back(std::make_pair(V3i(-1, -1, -1), 0.5*0.5*0.25));
		stamps[0].push_back(std::make_pair(V3i(-1, -1, 0), 0.5*0.5*0.75));
		stamps[0].push_back(std::make_pair(V3i(-1, -1, 1), 0.5*0.5*0.75));
		stamps[0].push_back(std::make_pair(V3i(-1, -1, 2), 0.5*0.5*0.25));
		stamps[0].push_back(std::make_pair(V3i(-1, 0, -1), 0.5*1.0*0.25));
		stamps[0].push_back(std::make_pair(V3i(-1, 0, 0), 0.5*1.0*0.75));
		stamps[0].push_back(std::make_pair(V3i(-1, 0, 1), 0.5*1.0*0.75));
		stamps[0].push_back(std::make_pair(V3i(-1, 0, 2), 0.5*1.0*0.25));
		stamps[0].push_back(std::make_pair(V3i(-1, 1, -1), 0.5*0.5*0.25));
		stamps[0].push_back(std::make_pair(V3i(-1, 1, 0), 0.5*0.5*0.75));
		stamps[0].push_back(std::make_pair(V3i(-1, 1, 1), 0.5*0.5*0.75));
		stamps[0].push_back(std::make_pair(V3i(-1, 1, 2), 0.5*0.5*0.25));
		stamps[0].push_back(std::make_pair(V3i(0, -1, -1), 1.0*0.5*0.25));
		stamps[0].push_back(std::make_pair(V3i(0, -1, 0), 1.0*0.5*0.75));
		stamps[0].push_back(std::make_pair(V3i(0, -1, 1), 1.0*0.5*0.75));
		stamps[0].push_back(std::make_pair(V3i(0, -1, 2), 1.0*0.5*0.25));
		stamps[0].push_back(std::make_pair(V3i(0, 0, -1), 1.0*1.0*0.25));
		stamps[0].push_back(std::make_pair(V3i(0, 0, 0), 1.0*1.0*0.75));
		stamps[0].push_back(std::make_pair(V3i(0, 0, 1), 1.0*1.0*0.75));
		stamps[0].push_back(std::make_pair(V3i(0, 0, 2), 1.0*1.0*0.25));
		stamps[0].push_back(std::make_pair(V3i(0, 1, -1), 1.0*0.5*0.25));
		stamps[0].push_back(std::make_pair(V3i(0, 1, 0), 1.0*0.5*0.75));
		stamps[0].push_back(std::make_pair(V3i(0, 1, 1), 1.0*0.5*0.75));
		stamps[0].push_back(std::make_pair(V3i(0, 1, 2), 1.0*0.5*0.25));
		stamps[0].push_back(std::make_pair(V3i(1, -1, -1), 0.5*0.5*0.25));
		stamps[0].push_back(std::make_pair(V3i(1, -1, 0), 0.5*0.5*0.75));
		stamps[0].push_back(std::make_pair(V3i(1, -1, 1), 0.5*0.5*0.75));
		stamps[0].push_back(std::make_pair(V3i(1, -1, 2), 0.5*0.5*0.25));
		stamps[0].push_back(std::make_pair(V3i(1, 0, -1), 0.5*1.0*0.25));
		stamps[0].push_back(std::make_pair(V3i(1, 0, 0), 0.5*1.0*0.75));
		stamps[0].push_back(std::make_pair(V3i(1, 0, 1), 0.5*1.0*0.75));
		stamps[0].push_back(std::make_pair(V3i(1, 0, 2), 0.5*1.0*0.25));
		stamps[0].push_back(std::make_pair(V3i(1, 1, -1), 0.5*0.5*0.25));
		stamps[0].push_back(std::make_pair(V3i(1, 1, 0), 0.5*0.5*0.75));
		stamps[0].push_back(std::make_pair(V3i(1, 1, 1), 0.5*0.5*0.75));
		stamps[0].push_back(std::make_pair(V3i(1, 1, 2), 0.5*0.5*0.25));
		stamps[1].push_back(std::make_pair(V3i(-1, -1, -1), 0.25*0.5*0.25));
		stamps[1].push_back(std::make_pair(V3i(-1, -1, 0), 0.25*0.5*0.75));
		stamps[1].push_back(std::make_pair(V3i(-1, -1, 1), 0.25*0.5*0.75));
		stamps[1].push_back(std::make_pair(V3i(-1, -1, 2), 0.25*0.5*0.25));
		stamps[1].push_back(std::make_pair(V3i(-1, 0, -1), 0.25*1.0*0.25));
		stamps[1].push_back(std::make_pair(V3i(-1, 0, 0), 0.25*1.0*0.75));
		stamps[1].push_back(std::make_pair(V3i(-1, 0, 1), 0.25*1.0*0.75));
		stamps[1].push_back(std::make_pair(V3i(-1, 0, 2), 0.25*1.0*0.25));
		stamps[1].push_back(std::make_pair(V3i(-1, 1, -1), 0.25*0.5*0.25));
		stamps[1].push_back(std::make_pair(V3i(-1, 1, 0), 0.25*0.5*0.75));
		stamps[1].push_back(std::make_pair(V3i(-1, 1, 1), 0.25*0.5*0.75));
		stamps[1].push_back(std::make_pair(V3i(-1, 1, 2), 0.25*0.5*0.25));
		stamps[1].push_back(std::make_pair(V3i(0, -1, -1), 0.75*0.5*0.25));
		stamps[1].push_back(std::make_pair(V3i(0, -1, 0), 0.75*0.5*0.75));
		stamps[1].push_back(std::make_pair(V3i(0, -1, 1), 0.75*0.5*0.75));
		stamps[1].push_back(std::make_pair(V3i(0, -1, 2), 0.75*0.5*0.25));
		stamps[1].push_back(std::make_pair(V3i(0, 0, -1), 0.75*1.0*0.25));
		stamps[1].push_back(std::make_pair(V3i(0, 0, 0), 0.75*1.0*0.75));
		stamps[1].push_back(std::make_pair(V3i(0, 0, 1), 0.75*1.0*0.75));
		stamps[1].push_back(std::make_pair(V3i(0, 0, 2), 0.75*1.0*0.25));
		stamps[1].push_back(std::make_pair(V3i(0, 1, -1), 0.75*0.5*0.25));
		stamps[1].push_back(std::make_pair(V3i(0, 1, 0), 0.75*0.5*0.75));
		stamps[1].push_back(std::make_pair(V3i(0, 1, 1), 0.75*0.5*0.75));
		stamps[1].push_back(std::make_pair(V3i(0, 1, 2), 0.75*0.5*0.25));
		stamps[1].push_back(std::make_pair(V3i(1, -1, -1), 0.75*0.5*0.25));
		stamps[1].push_back(std::make_pair(V3i(1, -1, 0), 0.75*0.5*0.75));
		stamps[1].push_back(std::make_pair(V3i(1, -1, 1), 0.75*0.5*0.75));
		stamps[1].push_back(std::make_pair(V3i(1, -1, 2), 0.75*0.5*0.25));
		stamps[1].push_back(std::make_pair(V3i(1, 0, -1), 0.75*1.0*0.25));
		stamps[1].push_back(std::make_pair(V3i(1, 0, 0), 0.75*1.0*0.75));
		stamps[1].push_back(std::make_pair(V3i(1, 0, 1), 0.75*1.0*0.75));
		stamps[1].push_back(std::make_pair(V3i(1, 0, 2), 0.75*1.0*0.25));
		stamps[1].push_back(std::make_pair(V3i(1, 1, -1), 0.75*0.5*0.25));
		stamps[1].push_back(std::make_pair(V3i(1, 1, 0), 0.75*0.5*0.75));
		stamps[1].push_back(std::make_pair(V3i(1, 1, 1), 0.75*0.5*0.75));
		stamps[1].push_back(std::make_pair(V3i(1, 1, 2), 0.75*0.5*0.25));
		stamps[1].push_back(std::make_pair(V3i(2, -1, -1), 0.25*0.5*0.25));
		stamps[1].push_back(std::make_pair(V3i(2, -1, 0), 0.25*0.5*0.75));
		stamps[1].push_back(std::make_pair(V3i(2, -1, 1), 0.25*0.5*0.75));
		stamps[1].push_back(std::make_pair(V3i(2, -1, 2), 0.25*0.5*0.25));
		stamps[1].push_back(std::make_pair(V3i(2, 0, -1), 0.25*1.0*0.25));
		stamps[1].push_back(std::make_pair(V3i(2, 0, 0), 0.25*1.0*0.75));
		stamps[1].push_back(std::make_pair(V3i(2, 0, 1), 0.25*1.0*0.75));
		stamps[1].push_back(std::make_pair(V3i(2, 0, 2), 0.25*1.0*0.25));
		stamps[1].push_back(std::make_pair(V3i(2, 1, -1), 0.25*0.5*0.25));
		stamps[1].push_back(std::make_pair(V3i(2, 1, 0), 0.25*0.5*0.75));
		stamps[1].push_back(std::make_pair(V3i(2, 1, 1), 0.25*0.5*0.75));
		stamps[1].push_back(std::make_pair(V3i(2, 1, 2), 0.25*0.5*0.25));
		stamps[2].push_back(std::make_pair(V3i(-1, -1, -1), 0.25*0.25*0.25));
		stamps[2].push_back(std::make_pair(V3i(-1, -1, 0), 0.25*0.25*0.75));
		stamps[2].push_back(std::make_pair(V3i(-1, -1, 1), 0.25*0.25*0.75));
		stamps[2].push_back(std::make_pair(V3i(-1, -1, 2), 0.25*0.25*0.25));
		stamps[2].push_back(std::make_pair(V3i(-1, 0, -1), 0.25*0.75*0.25));
		stamps[2].push_back(std::make_pair(V3i(-1, 0, 0), 0.25*0.75*0.75));
		stamps[2].push_back(std::make_pair(V3i(-1, 0, 1), 0.25*0.75*0.75));
		stamps[2].push_back(std::make_pair(V3i(-1, 0, 2), 0.25*0.75*0.25));
		stamps[2].push_back(std::make_pair(V3i(-1, 1, -1), 0.25*0.75*0.25));
		stamps[2].push_back(std::make_pair(V3i(-1, 1, 0), 0.25*0.75*0.75));
		stamps[2].push_back(std::make_pair(V3i(-1, 1, 1), 0.25*0.75*0.75));
		stamps[2].push_back(std::make_pair(V3i(-1, 1, 2), 0.25*0.75*0.25));
		stamps[2].push_back(std::make_pair(V3i(-1, 2, -1), 0.25*0.25*0.25));
		stamps[2].push_back(std::make_pair(V3i(-1, 2, 0), 0.25*0.25*0.75));
		stamps[2].push_back(std::make_pair(V3i(-1, 2, 1), 0.25*0.25*0.75));
		stamps[2].push_back(std::make_pair(V3i(-1, 2, 2), 0.25*0.25*0.25));
		stamps[2].push_back(std::make_pair(V3i(0, -1, -1), 0.75*0.25*0.25));
		stamps[2].push_back(std::make_pair(V3i(0, -1, 0), 0.75*0.25*0.75));
		stamps[2].push_back(std::make_pair(V3i(0, -1, 1), 0.75*0.25*0.75));
		stamps[2].push_back(std::make_pair(V3i(0, -1, 2), 0.75*0.25*0.25));
		stamps[2].push_back(std::make_pair(V3i(0, 0, -1), 0.75*0.75*0.25));
		stamps[2].push_back(std::make_pair(V3i(0, 0, 0), 0.75*0.75*0.75));
		stamps[2].push_back(std::make_pair(V3i(0, 0, 1), 0.75*0.75*0.75));
		stamps[2].push_back(std::make_pair(V3i(0, 0, 2), 0.75*0.75*0.25));
		stamps[2].push_back(std::make_pair(V3i(0, 1, -1), 0.75*0.75*0.25));
		stamps[2].push_back(std::make_pair(V3i(0, 1, 0), 0.75*0.75*0.75));
		stamps[2].push_back(std::make_pair(V3i(0, 1, 1), 0.75*0.75*0.75));
		stamps[2].push_back(std::make_pair(V3i(0, 1, 2), 0.75*0.75*0.25));
		stamps[2].push_back(std::make_pair(V3i(0, 2, -1), 0.75*0.25*0.25));
		stamps[2].push_back(std::make_pair(V3i(0, 2, 0), 0.75*0.25*0.75));
		stamps[2].push_back(std::make_pair(V3i(0, 2, 1), 0.75*0.25*0.75));
		stamps[2].push_back(std::make_pair(V3i(0, 2, 2), 0.75*0.25*0.25));
		stamps[2].push_back(std::make_pair(V3i(1, -1, -1), 0.75*0.25*0.25));
		stamps[2].push_back(std::make_pair(V3i(1, -1, 0), 0.75*0.25*0.75));
		stamps[2].push_back(std::make_pair(V3i(1, -1, 1), 0.75*0.25*0.75));
		stamps[2].push_back(std::make_pair(V3i(1, -1, 2), 0.75*0.25*0.25));
		stamps[2].push_back(std::make_pair(V3i(1, 0, -1), 0.75*0.75*0.25));
		stamps[2].push_back(std::make_pair(V3i(1, 0, 0), 0.75*0.75*0.75));
		stamps[2].push_back(std::make_pair(V3i(1, 0, 1), 0.75*0.75*0.75));
		stamps[2].push_back(std::make_pair(V3i(1, 0, 2), 0.75*0.75*0.25));
		stamps[2].push_back(std::make_pair(V3i(1, 1, -1), 0.75*0.75*0.25));
		stamps[2].push_back(std::make_pair(V3i(1, 1, 0), 0.75*0.75*0.75));
		stamps[2].push_back(std::make_pair(V3i(1, 1, 1), 0.75*0.75*0.75));
		stamps[2].push_back(std::make_pair(V3i(1, 1, 2), 0.75*0.75*0.25));
		stamps[2].push_back(std::make_pair(V3i(1, 2, -1), 0.75*0.25*0.25));
		stamps[2].push_back(std::make_pair(V3i(1, 2, 0), 0.75*0.25*0.75));
		stamps[2].push_back(std::make_pair(V3i(1, 2, 1), 0.75*0.25*0.75));
		stamps[2].push_back(std::make_pair(V3i(1, 2, 2), 0.75*0.25*0.25));
		stamps[2].push_back(std::make_pair(V3i(2, -1, -1), 0.25*0.25*0.25));
		stamps[2].push_back(std::make_pair(V3i(2, -1, 0), 0.25*0.25*0.75));
		stamps[2].push_back(std::make_pair(V3i(2, -1, 1), 0.25*0.25*0.75));
		stamps[2].push_back(std::make_pair(V3i(2, -1, 2), 0.25*0.25*0.25));
		stamps[2].push_back(std::make_pair(V3i(2, 0, -1), 0.25*0.75*0.25));
		stamps[2].push_back(std::make_pair(V3i(2, 0, 0), 0.25*0.75*0.75));
		stamps[2].push_back(std::make_pair(V3i(2, 0, 1), 0.25*0.75*0.75));
		stamps[2].push_back(std::make_pair(V3i(2, 0, 2), 0.25*0.75*0.25));
		stamps[2].push_back(std::make_pair(V3i(2, 1, -1), 0.25*0.75*0.25));
		stamps[2].push_back(std::make_pair(V3i(2, 1, 0), 0.25*0.75*0.75));
		stamps[2].push_back(std::make_pair(V3i(2, 1, 1), 0.25*0.75*0.75));
		stamps[2].push_back(std::make_pair(V3i(2, 1, 2), 0.25*0.75*0.25));
		stamps[2].push_back(std::make_pair(V3i(2, 2, -1), 0.25*0.25*0.25));
		stamps[2].push_back(std::make_pair(V3i(2, 2, 0), 0.25*0.25*0.75));
		stamps[2].push_back(std::make_pair(V3i(2, 2, 1), 0.25*0.25*0.75));
		stamps[2].push_back(std::make_pair(V3i(2, 2, 2), 0.25*0.25*0.25));
		stamps[3].push_back(std::make_pair(V3i(-1, -1, -1), 0.5*0.25*0.25));
		stamps[3].push_back(std::make_pair(V3i(-1, -1, 0), 0.5*0.25*0.75));
		stamps[3].push_back(std::make_pair(V3i(-1, -1, 1), 0.5*0.25*0.75));
		stamps[3].push_back(std::make_pair(V3i(-1, -1, 2), 0.5*0.25*0.25));
		stamps[3].push_back(std::make_pair(V3i(-1, 0, -1), 0.5*0.75*0.25));
		stamps[3].push_back(std::make_pair(V3i(-1, 0, 0), 0.5*0.75*0.75));
		stamps[3].push_back(std::make_pair(V3i(-1, 0, 1), 0.5*0.75*0.75));
		stamps[3].push_back(std::make_pair(V3i(-1, 0, 2), 0.5*0.75*0.25));
		stamps[3].push_back(std::make_pair(V3i(-1, 1, -1), 0.5*0.75*0.25));
		stamps[3].push_back(std::make_pair(V3i(-1, 1, 0), 0.5*0.75*0.75));
		stamps[3].push_back(std::make_pair(V3i(-1, 1, 1), 0.5*0.75*0.75));
		stamps[3].push_back(std::make_pair(V3i(-1, 1, 2), 0.5*0.75*0.25));
		stamps[3].push_back(std::make_pair(V3i(-1, 2, -1), 0.5*0.25*0.25));
		stamps[3].push_back(std::make_pair(V3i(-1, 2, 0), 0.5*0.25*0.75));
		stamps[3].push_back(std::make_pair(V3i(-1, 2, 1), 0.5*0.25*0.75));
		stamps[3].push_back(std::make_pair(V3i(-1, 2, 2), 0.5*0.25*0.25));
		stamps[3].push_back(std::make_pair(V3i(0, -1, -1), 1.0*0.25*0.25));
		stamps[3].push_back(std::make_pair(V3i(0, -1, 0), 1.0*0.25*0.75));
		stamps[3].push_back(std::make_pair(V3i(0, -1, 1), 1.0*0.25*0.75));
		stamps[3].push_back(std::make_pair(V3i(0, -1, 2), 1.0*0.25*0.25));
		stamps[3].push_back(std::make_pair(V3i(0, 0, -1), 1.0*0.75*0.25));
		stamps[3].push_back(std::make_pair(V3i(0, 0, 0), 1.0*0.75*0.75));
		stamps[3].push_back(std::make_pair(V3i(0, 0, 1), 1.0*0.75*0.75));
		stamps[3].push_back(std::make_pair(V3i(0, 0, 2), 1.0*0.75*0.25));
		stamps[3].push_back(std::make_pair(V3i(0, 1, -1), 1.0*0.75*0.25));
		stamps[3].push_back(std::make_pair(V3i(0, 1, 0), 1.0*0.75*0.75));
		stamps[3].push_back(std::make_pair(V3i(0, 1, 1), 1.0*0.75*0.75));
		stamps[3].push_back(std::make_pair(V3i(0, 1, 2), 1.0*0.75*0.25));
		stamps[3].push_back(std::make_pair(V3i(0, 2, -1), 1.0*0.25*0.25));
		stamps[3].push_back(std::make_pair(V3i(0, 2, 0), 1.0*0.25*0.75));
		stamps[3].push_back(std::make_pair(V3i(0, 2, 1), 1.0*0.25*0.75));
		stamps[3].push_back(std::make_pair(V3i(0, 2, 2), 1.0*0.25*0.25));
		stamps[3].push_back(std::make_pair(V3i(1, -1, -1), 0.5*0.25*0.25));
		stamps[3].push_back(std::make_pair(V3i(1, -1, 0), 0.5*0.25*0.75));
		stamps[3].push_back(std::make_pair(V3i(1, -1, 1), 0.5*0.25*0.75));
		stamps[3].push_back(std::make_pair(V3i(1, -1, 2), 0.5*0.25*0.25));
		stamps[3].push_back(std::make_pair(V3i(1, 0, -1), 0.5*0.75*0.25));
		stamps[3].push_back(std::make_pair(V3i(1, 0, 0), 0.5*0.75*0.75));
		stamps[3].push_back(std::make_pair(V3i(1, 0, 1), 0.5*0.75*0.75));
		stamps[3].push_back(std::make_pair(V3i(1, 0, 2), 0.5*0.75*0.25));
		stamps[3].push_back(std::make_pair(V3i(1, 1, -1), 0.5*0.75*0.25));
		stamps[3].push_back(std::make_pair(V3i(1, 1, 0), 0.5*0.75*0.75));
		stamps[3].push_back(std::make_pair(V3i(1, 1, 1), 0.5*0.75*0.75));
		stamps[3].push_back(std::make_pair(V3i(1, 1, 2), 0.5*0.75*0.25));
		stamps[3].push_back(std::make_pair(V3i(1, 2, -1), 0.5*0.25*0.25));
		stamps[3].push_back(std::make_pair(V3i(1, 2, 0), 0.5*0.25*0.75));
		stamps[3].push_back(std::make_pair(V3i(1, 2, 1), 0.5*0.25*0.75));
		stamps[3].push_back(std::make_pair(V3i(1, 2, 2), 0.5*0.25*0.25));
		stamps[4].push_back(std::make_pair(V3i(-1, -1, -1), 0.5*0.5*0.5));
		stamps[4].push_back(std::make_pair(V3i(-1, -1, 0), 0.5*0.5*1.0));
		stamps[4].push_back(std::make_pair(V3i(-1, -1, 1), 0.5*0.5*0.5));
		stamps[4].push_back(std::make_pair(V3i(-1, 0, -1), 0.5*1.0*0.5));
		stamps[4].push_back(std::make_pair(V3i(-1, 0, 0), 0.5*1.0*1.0));
		stamps[4].push_back(std::make_pair(V3i(-1, 0, 1), 0.5*1.0*0.5));
		stamps[4].push_back(std::make_pair(V3i(-1, 1, -1), 0.5*0.5*0.5));
		stamps[4].push_back(std::make_pair(V3i(-1, 1, 0), 0.5*0.5*1.0));
		stamps[4].push_back(std::make_pair(V3i(-1, 1, 1), 0.5*0.5*0.5));
		stamps[4].push_back(std::make_pair(V3i(0, -1, -1), 1.0*0.5*0.5));
		stamps[4].push_back(std::make_pair(V3i(0, -1, 0), 1.0*0.5*1.0));
		stamps[4].push_back(std::make_pair(V3i(0, -1, 1), 1.0*0.5*0.5));
		stamps[4].push_back(std::make_pair(V3i(0, 0, -1), 1.0*1.0*0.5));
		stamps[4].push_back(std::make_pair(V3i(0, 0, 0), 1.0*1.0*1.0));
		stamps[4].push_back(std::make_pair(V3i(0, 0, 1), 1.0*1.0*0.5));
		stamps[4].push_back(std::make_pair(V3i(0, 1, -1), 1.0*0.5*0.5));
		stamps[4].push_back(std::make_pair(V3i(0, 1, 0), 1.0*0.5*1.0));
		stamps[4].push_back(std::make_pair(V3i(0, 1, 1), 1.0*0.5*0.5));
		stamps[4].push_back(std::make_pair(V3i(1, -1, -1), 0.5*0.5*0.5));
		stamps[4].push_back(std::make_pair(V3i(1, -1, 0), 0.5*0.5*1.0));
		stamps[4].push_back(std::make_pair(V3i(1, -1, 1), 0.5*0.5*0.5));
		stamps[4].push_back(std::make_pair(V3i(1, 0, -1), 0.5*1.0*0.5));
		stamps[4].push_back(std::make_pair(V3i(1, 0, 0), 0.5*1.0*1.0));
		stamps[4].push_back(std::make_pair(V3i(1, 0, 1), 0.5*1.0*0.5));
		stamps[4].push_back(std::make_pair(V3i(1, 1, -1), 0.5*0.5*0.5));
		stamps[4].push_back(std::make_pair(V3i(1, 1, 0), 0.5*0.5*1.0));
		stamps[4].push_back(std::make_pair(V3i(1, 1, 1), 0.5*0.5*0.5));
		stamps[5].push_back(std::make_pair(V3i(-1, -1, -1), 0.25*0.5*0.5));
		stamps[5].push_back(std::make_pair(V3i(-1, -1, 0), 0.25*0.5*1.0));
		stamps[5].push_back(std::make_pair(V3i(-1, -1, 1), 0.25*0.5*0.5));
		stamps[5].push_back(std::make_pair(V3i(-1, 0, -1), 0.25*1.0*0.5));
		stamps[5].push_back(std::make_pair(V3i(-1, 0, 0), 0.25*1.0*1.0));
		stamps[5].push_back(std::make_pair(V3i(-1, 0, 1), 0.25*1.0*0.5));
		stamps[5].push_back(std::make_pair(V3i(-1, 1, -1), 0.25*0.5*0.5));
		stamps[5].push_back(std::make_pair(V3i(-1, 1, 0), 0.25*0.5*1.0));
		stamps[5].push_back(std::make_pair(V3i(-1, 1, 1), 0.25*0.5*0.5));
		stamps[5].push_back(std::make_pair(V3i(0, -1, -1), 0.75*0.5*0.5));
		stamps[5].push_back(std::make_pair(V3i(0, -1, 0), 0.75*0.5*1.0));
		stamps[5].push_back(std::make_pair(V3i(0, -1, 1), 0.75*0.5*0.5));
		stamps[5].push_back(std::make_pair(V3i(0, 0, -1), 0.75*1.0*0.5));
		stamps[5].push_back(std::make_pair(V3i(0, 0, 0), 0.75*1.0*1.0));
		stamps[5].push_back(std::make_pair(V3i(0, 0, 1), 0.75*1.0*0.5));
		stamps[5].push_back(std::make_pair(V3i(0, 1, -1), 0.75*0.5*0.5));
		stamps[5].push_back(std::make_pair(V3i(0, 1, 0), 0.75*0.5*1.0));
		stamps[5].push_back(std::make_pair(V3i(0, 1, 1), 0.75*0.5*0.5));
		stamps[5].push_back(std::make_pair(V3i(1, -1, -1), 0.75*0.5*0.5));
		stamps[5].push_back(std::make_pair(V3i(1, -1, 0), 0.75*0.5*1.0));
		stamps[5].push_back(std::make_pair(V3i(1, -1, 1), 0.75*0.5*0.5));
		stamps[5].push_back(std::make_pair(V3i(1, 0, -1), 0.75*1.0*0.5));
		stamps[5].push_back(std::make_pair(V3i(1, 0, 0), 0.75*1.0*1.0));
		stamps[5].push_back(std::make_pair(V3i(1, 0, 1), 0.75*1.0*0.5));
		stamps[5].push_back(std::make_pair(V3i(1, 1, -1), 0.75*0.5*0.5));
		stamps[5].push_back(std::make_pair(V3i(1, 1, 0), 0.75*0.5*1.0));
		stamps[5].push_back(std::make_pair(V3i(1, 1, 1), 0.75*0.5*0.5));
		stamps[5].push_back(std::make_pair(V3i(2, -1, -1), 0.25*0.5*0.5));
		stamps[5].push_back(std::make_pair(V3i(2, -1, 0), 0.25*0.5*1.0));
		stamps[5].push_back(std::make_pair(V3i(2, -1, 1), 0.25*0.5*0.5));
		stamps[5].push_back(std::make_pair(V3i(2, 0, -1), 0.25*1.0*0.5));
		stamps[5].push_back(std::make_pair(V3i(2, 0, 0), 0.25*1.0*1.0));
		stamps[5].push_back(std::make_pair(V3i(2, 0, 1), 0.25*1.0*0.5));
		stamps[5].push_back(std::make_pair(V3i(2, 1, -1), 0.25*0.5*0.5));
		stamps[5].push_back(std::make_pair(V3i(2, 1, 0), 0.25*0.5*1.0));
		stamps[5].push_back(std::make_pair(V3i(2, 1, 1), 0.25*0.5*0.5));
		stamps[6].push_back(std::make_pair(V3i(-1, -1, -1), 0.25*0.25*0.5));
		stamps[6].push_back(std::make_pair(V3i(-1, -1, 0), 0.25*0.25*1.0));
		stamps[6].push_back(std::make_pair(V3i(-1, -1, 1), 0.25*0.25*0.5));
		stamps[6].push_back(std::make_pair(V3i(-1, 0, -1), 0.25*0.75*0.5));
		stamps[6].push_back(std::make_pair(V3i(-1, 0, 0), 0.25*0.75*1.0));
		stamps[6].push_back(std::make_pair(V3i(-1, 0, 1), 0.25*0.75*0.5));
		stamps[6].push_back(std::make_pair(V3i(-1, 1, -1), 0.25*0.75*0.5));
		stamps[6].push_back(std::make_pair(V3i(-1, 1, 0), 0.25*0.75*1.0));
		stamps[6].push_back(std::make_pair(V3i(-1, 1, 1), 0.25*0.75*0.5));
		stamps[6].push_back(std::make_pair(V3i(-1, 2, -1), 0.25*0.25*0.5));
		stamps[6].push_back(std::make_pair(V3i(-1, 2, 0), 0.25*0.25*1.0));
		stamps[6].push_back(std::make_pair(V3i(-1, 2, 1), 0.25*0.25*0.5));
		stamps[6].push_back(std::make_pair(V3i(0, -1, -1), 0.75*0.25*0.5));
		stamps[6].push_back(std::make_pair(V3i(0, -1, 0), 0.75*0.25*1.0));
		stamps[6].push_back(std::make_pair(V3i(0, -1, 1), 0.75*0.25*0.5));
		stamps[6].push_back(std::make_pair(V3i(0, 0, -1), 0.75*0.75*0.5));
		stamps[6].push_back(std::make_pair(V3i(0, 0, 0), 0.75*0.75*1.0));
		stamps[6].push_back(std::make_pair(V3i(0, 0, 1), 0.75*0.75*0.5));
		stamps[6].push_back(std::make_pair(V3i(0, 1, -1), 0.75*0.75*0.5));
		stamps[6].push_back(std::make_pair(V3i(0, 1, 0), 0.75*0.75*1.0));
		stamps[6].push_back(std::make_pair(V3i(0, 1, 1), 0.75*0.75*0.5));
		stamps[6].push_back(std::make_pair(V3i(0, 2, -1), 0.75*0.25*0.5));
		stamps[6].push_back(std::make_pair(V3i(0, 2, 0), 0.75*0.25*1.0));
		stamps[6].push_back(std::make_pair(V3i(0, 2, 1), 0.75*0.25*0.5));
		stamps[6].push_back(std::make_pair(V3i(1, -1, -1), 0.75*0.25*0.5));
		stamps[6].push_back(std::make_pair(V3i(1, -1, 0), 0.75*0.25*1.0));
		stamps[6].push_back(std::make_pair(V3i(1, -1, 1), 0.75*0.25*0.5));
		stamps[6].push_back(std::make_pair(V3i(1, 0, -1), 0.75*0.75*0.5));
		stamps[6].push_back(std::make_pair(V3i(1, 0, 0), 0.75*0.75*1.0));
		stamps[6].push_back(std::make_pair(V3i(1, 0, 1), 0.75*0.75*0.5));
		stamps[6].push_back(std::make_pair(V3i(1, 1, -1), 0.75*0.75*0.5));
		stamps[6].push_back(std::make_pair(V3i(1, 1, 0), 0.75*0.75*1.0));
		stamps[6].push_back(std::make_pair(V3i(1, 1, 1), 0.75*0.75*0.5));
		stamps[6].push_back(std::make_pair(V3i(1, 2, -1), 0.75*0.25*0.5));
		stamps[6].push_back(std::make_pair(V3i(1, 2, 0), 0.75*0.25*1.0));
		stamps[6].push_back(std::make_pair(V3i(1, 2, 1), 0.75*0.25*0.5));
		stamps[6].push_back(std::make_pair(V3i(2, -1, -1), 0.25*0.25*0.5));
		stamps[6].push_back(std::make_pair(V3i(2, -1, 0), 0.25*0.25*1.0));
		stamps[6].push_back(std::make_pair(V3i(2, -1, 1), 0.25*0.25*0.5));
		stamps[6].push_back(std::make_pair(V3i(2, 0, -1), 0.25*0.75*0.5));
		stamps[6].push_back(std::make_pair(V3i(2, 0, 0), 0.25*0.75*1.0));
		stamps[6].push_back(std::make_pair(V3i(2, 0, 1), 0.25*0.75*0.5));
		stamps[6].push_back(std::make_pair(V3i(2, 1, -1), 0.25*0.75*0.5));
		stamps[6].push_back(std::make_pair(V3i(2, 1, 0), 0.25*0.75*1.0));
		stamps[6].push_back(std::make_pair(V3i(2, 1, 1), 0.25*0.75*0.5));
		stamps[6].push_back(std::make_pair(V3i(2, 2, -1), 0.25*0.25*0.5));
		stamps[6].push_back(std::make_pair(V3i(2, 2, 0), 0.25*0.25*1.0));
		stamps[6].push_back(std::make_pair(V3i(2, 2, 1), 0.25*0.25*0.5));
		stamps[7].push_back(std::make_pair(V3i(-1, -1, -1), 0.5*0.25*0.5));
		stamps[7].push_back(std::make_pair(V3i(-1, -1, 0), 0.5*0.25*1.0));
		stamps[7].push_back(std::make_pair(V3i(-1, -1, 1), 0.5*0.25*0.5));
		stamps[7].push_back(std::make_pair(V3i(-1, 0, -1), 0.5*0.75*0.5));
		stamps[7].push_back(std::make_pair(V3i(-1, 0, 0), 0.5*0.75*1.0));
		stamps[7].push_back(std::make_pair(V3i(-1, 0, 1), 0.5*0.75*0.5));
		stamps[7].push_back(std::make_pair(V3i(-1, 1, -1), 0.5*0.75*0.5));
		stamps[7].push_back(std::make_pair(V3i(-1, 1, 0), 0.5*0.75*1.0));
		stamps[7].push_back(std::make_pair(V3i(-1, 1, 1), 0.5*0.75*0.5));
		stamps[7].push_back(std::make_pair(V3i(-1, 2, -1), 0.5*0.25*0.5));
		stamps[7].push_back(std::make_pair(V3i(-1, 2, 0), 0.5*0.25*1.0));
		stamps[7].push_back(std::make_pair(V3i(-1, 2, 1), 0.5*0.25*0.5));
		stamps[7].push_back(std::make_pair(V3i(0, -1, -1), 1.0*0.25*0.5));
		stamps[7].push_back(std::make_pair(V3i(0, -1, 0), 1.0*0.25*1.0));
		stamps[7].push_back(std::make_pair(V3i(0, -1, 1), 1.0*0.25*0.5));
		stamps[7].push_back(std::make_pair(V3i(0, 0, -1), 1.0*0.75*0.5));
		stamps[7].push_back(std::make_pair(V3i(0, 0, 0), 1.0*0.75*1.0));
		stamps[7].push_back(std::make_pair(V3i(0, 0, 1), 1.0*0.75*0.5));
		stamps[7].push_back(std::make_pair(V3i(0, 1, -1), 1.0*0.75*0.5));
		stamps[7].push_back(std::make_pair(V3i(0, 1, 0), 1.0*0.75*1.0));
		stamps[7].push_back(std::make_pair(V3i(0, 1, 1), 1.0*0.75*0.5));
		stamps[7].push_back(std::make_pair(V3i(0, 2, -1), 1.0*0.25*0.5));
		stamps[7].push_back(std::make_pair(V3i(0, 2, 0), 1.0*0.25*1.0));
		stamps[7].push_back(std::make_pair(V3i(0, 2, 1), 1.0*0.25*0.5));
		stamps[7].push_back(std::make_pair(V3i(1, -1, -1), 0.5*0.25*0.5));
		stamps[7].push_back(std::make_pair(V3i(1, -1, 0), 0.5*0.25*1.0));
		stamps[7].push_back(std::make_pair(V3i(1, -1, 1), 0.5*0.25*0.5));
		stamps[7].push_back(std::make_pair(V3i(1, 0, -1), 0.5*0.75*0.5));
		stamps[7].push_back(std::make_pair(V3i(1, 0, 0), 0.5*0.75*1.0));
		stamps[7].push_back(std::make_pair(V3i(1, 0, 1), 0.5*0.75*0.5));
		stamps[7].push_back(std::make_pair(V3i(1, 1, -1), 0.5*0.75*0.5));
		stamps[7].push_back(std::make_pair(V3i(1, 1, 0), 0.5*0.75*1.0));
		stamps[7].push_back(std::make_pair(V3i(1, 1, 1), 0.5*0.75*0.5));
		stamps[7].push_back(std::make_pair(V3i(1, 2, -1), 0.5*0.25*0.5));
		stamps[7].push_back(std::make_pair(V3i(1, 2, 0), 0.5*0.25*1.0));
		stamps[7].push_back(std::make_pair(V3i(1, 2, 1), 0.5*0.25*0.5));
	}

	// iterate all coarse voxels
	std::vector<PNSystem::Voxel>& voxels_coarse = vm_coarse.getVoxels();
	for( auto&v_coarse:voxels_coarse )
	{
		// coordinate of the corresponding fine voxel
		V3i coord_fine = v_coarse.coord*2;

		// Now compute the value of the coarse voxels coefficients from gathering the values from
		// the fine voxels which are affecting the coarse voxels values. These values will be weighted,
		// as fine voxels are affected by multiple coarse voxels.

		for( int coeff_index=0;coeff_index<sys_fine.getNumCoefficients();++coeff_index )
		{
			//int coeff_index = 0;

			// get the global index of the coarse voxel/coeff
			int gi_coarse = vm_coarse.getGlobalIndex( v_coarse, coeff_index );

			// initialize the value of the current coefficient at current coarse voxel
			// this value is computed from all affected fine voxels
			double weight_sum = 0.0;

			// now iterate over all affected fine voxels. This depends on the grid location
			int grid_index = sys_fine.getStencil().getGridIndexFromCoefficient(coeff_index);
			auto& stamp = stamps[grid_index];
			std::vector<std::pair<int, double>> row_entries; // column index and value
			for( auto& tt:stamp )
			{
				V3i v_fine_coord = coord_fine + std::get<0>(tt);
				double weight = std::get<1>(tt);

				// using the stamps offset on boundary voxels may create invalid voxel coordinates
				if( !vm_fine.voxelIsValid(v_fine_coord) )
					continue;

				PNSystem::Voxel& v_fine = vm_fine.getVoxel(v_fine_coord);

				// retrieve global index of fine voxel/coeff
				int gi_fine = vm_fine.getGlobalIndex(v_fine, coeff_index);
				if( gi_fine >= 0 )
				{
					// either the fine voxel/coeff is not a boundary voxel/coeff,
					// or we have a boundary voxel with neumann BC, pointing to an
					// interior voxel
					row_entries.push_back( std::make_pair(gi_fine, weight) );
				}
				// only set weights on voxels which are actually active
				if( !vm_fine.isBoundaryCoefficient(v_fine, coeff_index)&&(gi_coarse>=0) )
					upsampleMatrixBuilder.coeff(gi_fine, gi_coarse) += weight;

				// else: value_fine remains zero, which means we have a boundary voxel with dirichlet BC
				weight_sum += weight;
			} // for each stamp entry

			// build row
			if( !vm_coarse.isBoundaryCoefficient(v_coarse, coeff_index) )
				for( auto& it:row_entries )
				{
					int gi_fine = it.first;
					double weight = it.second;
					downsampleMatrixBuilder.coeff(gi_coarse, it.first) += weight/weight_sum;
				}
		} // for each coefficient
	} // for each coarse voxel



	//return x_coarse;
	downsampleMatrix = downsampleMatrixBuilder.build( vm_coarse.getNumUnknowns(), vm_fine.getNumUnknowns() );
	upsampleMatrix = upsampleMatrixBuilder.build( vm_fine.getNumUnknowns(), vm_coarse.getNumUnknowns() );
}
/*
// runs a number of CG iterations on the given problem
// returns the square root of the residual
double run_cg_iterations( Eigen::SparseMatrix<double, Eigen::RowMajor>& A, Eigen::VectorXd& b, Eigen::VectorXd& x, Eigen::VectorXd& r, int numIterations, double tol )
{
	r = b-A*x;
	Eigen::VectorXd p = r;
	Eigen::VectorXd Ap;
	double rsold = r.squaredNorm();

	for( int i=0;i<numIterations;++i )
	{
		Ap = A*p;
		double alpha = rsold/p.dot(Ap);
		x = x + alpha*p;
		r = r - alpha*Ap;
		double rsnew = r.squaredNorm();
		if( std::sqrt(rsnew) < tol )
			return std::sqrt(rsnew);
		p = r + (rsnew/rsold)*p;
		rsold = rsnew;
	}

	return std::sqrt(rsold);
}

// runs a number of Gauss-Seidel iterations on the given problem
void run_gs_iterations( Eigen::SparseMatrix<double, Eigen::RowMajor> A, Eigen::VectorXd& b, Eigen::VectorXd& x, int numIterations )
{
	for( int k=0;k<numIterations;++k )
	{
		// iterate all rows
		for( int i=0;i<A.rows();++i )
		{
			double aii = 0.0;
			double sum = b(i);
			// iterate all non-zero column elements (this is why we need RowMajor storage order for A)
			for( Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(A, i);it;++it )
			{
				int j = it.col();
				if( i != j )
					sum -= it.value()*x.coeffRef(j);
				else
					aii = it.value();
			}
			x.coeffRef(i) = sum/aii;
		}
	}
}



std::tuple<Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd> solve_cg( PNSystem& sys )
{
	std::cout << "solve_cg: solving..." << "s\n";
	Timer timer;
	timer.start();

	std::vector<double> solve_convergence;
	std::vector<double> solve_convergence_timestamps;

	PNSystem::RealMatrix A = sys.get_A_real().transpose()*sys.get_A_real();
	Eigen::VectorXd b = sys.get_A_real().transpose()*sys.get_b_real();
	Eigen::VectorXd x = b; // initial guess
	Eigen::VectorXd r( sys.getVoxelManager().getNumUnknowns() );



	// cg solver
	double tol = 1.0e-10; // convergence error tolerance
	int numIterations = b.rows();
	{
		r = b-A*x;
		Eigen::VectorXd p = r;
		Eigen::VectorXd Ap;
		double rsold = r.squaredNorm();
		for( int i=0;i<numIterations;++i )
		{
			Ap = A*p;
			double alpha = rsold/p.dot(Ap);
			x = x + alpha*p;
			r = r - alpha*Ap;
			double rsnew = r.squaredNorm();
			solve_convergence.push_back(std::sqrt(rsnew));
			solve_convergence_timestamps.push_back(timer.elapsedSeconds());
			if( std::sqrt(rsnew) < tol )
				break;
			p = r + (rsnew/rsold)*p;
			rsold = rsnew;
		}
	} // end of cg solver

	timer.stop();
	std::cout << "solve_cg: " << timer.elapsedSeconds() << "s #iterations=" << solve_convergence.size() << "\n";
	return std::make_tuple( sys.stripBoundary(x),
							to_vector(solve_convergence),
							to_vector(solve_convergence_timestamps));
}


std::tuple<Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd> solve_sparseLU( PNSystem& sys )
{
	std::vector<double> solve_convergence;
	std::vector<double> solve_convergence_timestamps;


	Timer timer;
	timer.start();





	//Eigen::ConjugateGradient<RealMatrix> solver;
	//Eigen::BiCGSTAB<RealMatrix> solver;
	Eigen::SparseLU<PNSystem::RealMatrix> solver;
	solver.compute(sys.get_A_real());
	if(solver.info()!=Eigen::Success)
	{
		throw std::runtime_error("solve_sparseLU decomposition failed");
	}
	Eigen::VectorXd x = solver.solve(sys.get_b_real());
	if(solver.info()!=Eigen::Success)
	{
		throw std::runtime_error("solve_sparseLU solve failed");
	}


	timer.stop();
	std::cout << "solve_sparseLU: " << timer.elapsedSeconds() << "s\n";

	return std::make_tuple( sys.stripBoundary(x),
							to_vector(solve_convergence),
							to_vector(solve_convergence_timestamps));
}








std::tuple<Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd> solve_gs(PNSystem& sys)
{
	Eigen::SparseMatrix<double, Eigen::RowMajor> A = sys.get_A_real().transpose()*sys.get_A_real();
	Eigen::VectorXd b = sys.get_A_real().transpose()*sys.get_b_real();
	Eigen::VectorXd x = b; // initial guess
	Eigen::VectorXd r;

	std::vector<double> solve_convergence;
	std::vector<double> solve_convergence_timestamps;


	int numIterations = 1000;
	for( int k=0;k<numIterations;++k )
	{
		run_gs_iterations(A, b, x, 1);
		r = b - A*x;
		solve_convergence.push_back(r.norm());
		solve_convergence_timestamps.push_back(0.0);
	}

	return std::make_tuple( sys.stripBoundary(x),
							to_vector(solve_convergence),
							to_vector(solve_convergence_timestamps));

}

struct MultigridLevel
{
	MultigridLevel():
		next(0)
	{

	}

	Eigen::SparseMatrix<double, Eigen::RowMajor> A;
	Eigen::VectorXd      b;
	Eigen::VectorXd      x;
	Eigen::VectorXd      r;

	PNSystem::RealMatrix upsample; // converts x to the next coarser grid
	PNSystem::RealMatrix downsample; // converts from the next coarser grid to this grid

	MultigridLevel* next;
};

struct ProfileTimer
{
	Timer resampling;
	Timer smoothing;
	Timer lowest_level_solve;

	void print()
	{
		std::cout << "resampling time:" << resampling.elapsedSeconds() << "s\n";
		std::cout << "lowest level solve time:" << lowest_level_solve.elapsedSeconds() << "s\n";
		std::cout << "smoothing time:" << smoothing.elapsedSeconds() << "s\n";
	}
};


void multigrid_cycle( MultigridLevel* lvl_fine, ProfileTimer* timers = 0 )
{
	MultigridLevel* lvl_coarse = lvl_fine->next;

	// pre smoothing
	//run_cg_iterations( lvl_fine->A, lvl_fine->b, lvl_fine->x, lvl_fine->r, 1, 0.0 );
	if(timers)
		timers->smoothing.start();
	run_gs_iterations( lvl_fine->A, lvl_fine->b, lvl_fine->x, 5);
	if(timers)
		timers->smoothing.stop();


	// compute residual on fine level
	if(timers)
		timers->resampling.start();
	lvl_fine->r = lvl_fine->b - lvl_fine->A*lvl_fine->x;
	if(timers)
		timers->resampling.stop();
	// restriction
	lvl_coarse->b = lvl_fine->downsample*lvl_fine->r;

	// compute approximate solution to the correction equation on the coarser grid
	if( lvl_coarse->next == 0 )
	{
		// the coarse level is the last level...fully solve the thing
		if(timers)
			timers->lowest_level_solve.start();
		run_cg_iterations( lvl_coarse->A, lvl_coarse->b, lvl_coarse->x, lvl_coarse->r, 1000, 1.0e-10 );
		//run_gs_iterations( lvl_coarse->A, lvl_coarse->b, lvl_coarse->x, 1000, 1.0e-10 );
		if(timers)
			timers->lowest_level_solve.stop();
	}else
	{
		lvl_coarse->x.fill(0.0);
		for( int i=0;i<2;++i )
			multigrid_cycle( lvl_coarse, timers );
	}

	// upsample correction and apply
	if(timers)
		timers->resampling.start();
	lvl_fine->x = lvl_fine->x + lvl_fine->upsample*lvl_coarse->x;
	if(timers)
		timers->resampling.stop();

	// post smoothing
	//run_cg_iterations( lvl_fine->A, lvl_fine->b, lvl_fine->x, lvl_fine->r, 1, 0.0 );
	if(timers)
		timers->smoothing.start();
	run_gs_iterations( lvl_fine->A, lvl_fine->b, lvl_fine->x, 5);
	if(timers)
		timers->smoothing.stop();

}



std::tuple<Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd> solve_multigrid(PNSystem& sys, int numLevels)
{
	ProfileTimer timers;
	double tol = 1.0e-10; // convergence error tolerance
	//int numLevels = 9; // for res=512
	//int numLevels = 6; // for res=64
	Timer timer;
	std::vector<double> solve_convergence;
	std::vector<double> solve_convergence_timestamps;
	std::vector<MultigridLevel> levels(numLevels);

	levels[0].A = sys.get_A_real().transpose()*sys.get_A_real();
	levels[0].b = sys.get_A_real().transpose()*sys.get_b_real();
	levels[0].x = levels[0].b;
	levels[0].r = levels[0].b-levels[0].A*levels[0].x;
	buildUpAndDownsamplingMatrices(sys, levels[0].downsample, levels[0].upsample);


	PNSystem::Stencil& stencil = sys.getStencil();
	int boundaryConditions = sys.getBoundaryConditions();
	Domain domain = sys.getDomain();
	PNSystem::Fields fields = sys.getFields();

	// link up
	for( int i=0;i<levels.size()-1;++i )
	{
		domain = domain.downsample();
		fields = fields.createRestricted();
		PNSystem sys_coarse(stencil, domain, boundaryConditions);
		sys_coarse.setFields( fields );
		sys_coarse.build();

		levels[i+1].A = sys_coarse.get_A_real().transpose()*sys_coarse.get_A_real();
		levels[i+1].b = sys_coarse.get_A_real().transpose()*sys_coarse.get_b_real();
		levels[i+1].x = Eigen::VectorXd(sys_coarse.getVoxelManager().getNumUnknowns());
		levels[i+1].x.fill(0.0);
		levels[i+1].r = Eigen::VectorXd(sys_coarse.getVoxelManager().getNumUnknowns());
		buildUpAndDownsamplingMatrices(sys_coarse, levels[i+1].downsample, levels[i+1].upsample);

		levels[i].next = &levels[i+1];
	}

	timer.start();
	int maxIter = 1000;
	for( int i=0;i<maxIter;++i )
	{
		multigrid_cycle( &levels[0], &timers );
		double rmse = (levels[0].b-levels[0].A*levels[0].x).norm();
		solve_convergence.push_back( rmse );
		solve_convergence_timestamps.push_back(timer.elapsedSeconds());
		if(rmse<tol)
			break;
	}

	timer.stop();
	std::cout << "solve_multigrid_test: " << timer.elapsedSeconds() << "s #iterations=" << solve_convergence.size() << "\n";
	timers.print();
	return std::make_tuple( sys.stripBoundary(levels[0].x),
							to_vector(solve_convergence),
							to_vector(solve_convergence_timestamps));

}
*/

void setup_solver( MGTEST& mg, PNSystem& sys, int numLevels = 1 )
{
	mg = MGTEST(numLevels);

	PNSystem::Stencil& stencil = sys.getStencil();
	int boundaryConditions = sys.getBoundaryConditions();
	Domain domain = sys.getDomain();
	PNSystem::Fields fields = sys.getFields();

	for( int i=0;i<numLevels;++i )
	{
		PNSystem::RealMatrix A;
		Eigen::VectorXd x;
		PNSystem::RealMatrix downsample;
		PNSystem::RealMatrix upsample;

		PNSystem sys_level(stencil, domain, boundaryConditions);
		sys_level.setFields( fields );
		sys_level.build();
		buildUpAndDownsamplingMatrices(sys_level, downsample, upsample);

		A = sys_level.get_A_real().transpose()*sys_level.get_A_real();
		x = Eigen::VectorXd(sys_level.getVoxelManager().getNumUnknowns());
		x.fill(0.0);

		if( i==0 )
		{
			Eigen::VectorXd b = sys_level.get_A_real().transpose()*sys_level.get_b_real();
			mg.setb(b);
		}

		mg.setMultigridLevel(i, A, x, downsample, upsample);


		// downsample to next level
		domain = domain.downsample();
		fields = fields.createRestricted();
	}

	// saved for debugging
	sys.debug_downsample = mg.m_levels[0].downsample;
	sys.debug_upsample = mg.m_levels[0].upsample;
}

std::tuple<Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd> solve_multigrid(PNSystem& sys, int numLevels)
{
	MGTEST mg;
	setup_solver(mg, sys, numLevels);
	int maxIterations = 1000;
	auto result = mg.solve(maxIterations);

	sys.debug_x = std::get<0>(result);
	sys.debug_x_downsampled = sys.debug_downsample*sys.debug_x;
	sys.debug_x_up_sampled_downsampled = sys.debug_upsample*sys.debug_x_downsampled;


	return std::make_tuple(sys.stripBoundary(std::get<0>(result)), std::get<1>(result), std::get<2>(result));
}

std::tuple<Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd> solve_gs(PNSystem& sys)
{
	MGTEST mg;
	setup_solver(mg, sys);
	auto result = mg.solve_gs();
	return std::make_tuple(sys.stripBoundary(std::get<0>(result)), std::get<1>(result), std::get<2>(result));
}

std::tuple<Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd> solve_cg(PNSystem& sys)
{
	MGTEST mg;
	setup_solver(mg, sys);
	auto result = mg.solve_cg();
	return std::make_tuple(sys.stripBoundary(std::get<0>(result)), std::get<1>(result), std::get<2>(result));
}

std::tuple<Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd> solve_cg_eigen(PNSystem& sys)
{
	MGTEST mg;
	setup_solver(mg, sys);
	auto result = mg.solve_cg_eigen();
	return std::make_tuple(sys.stripBoundary(std::get<0>(result)), std::get<1>(result), std::get<2>(result));
}


PYBIND11_MODULE(pnsolver, m)
{
	m.def( "solve_multigrid", &solve_multigrid);
	m.def( "solve_cg", &solve_cg);
	m.def( "solve_cg_eigen", &solve_cg_eigen);
	m.def( "solve_gs", &solve_gs);
	//m.def( "solve_sparseLU", &solve_sparseLU);
	m.def( "physical_mem_used_by_process", &physical_mem_used_by_process);



	// MGTEST ================================
	py::class_<MGTEST> class_mgtest(m, "MGTEST");
	class_mgtest
	.def("__init__",
	[](MGTEST &m, int numLevels)
	{
		new (&m) MGTEST(numLevels);
	})
	.def("setRef",&MGTEST::setRef)
	.def("setMultigridLevel",&MGTEST::setMultigridLevel)
	.def("setb",&MGTEST::setb)
	.def("solve", &MGTEST::solve)
	.def("solve_gs", &MGTEST::solve_gs)
	.def("solve_cg", &MGTEST::solve_cg)
	.def("memoryTest", &MGTEST::memoryTest)
	;



	// PNSystem ==============================
	py::class_<PNSystem> class_pnsystem(m, "PNSystem");
	class_pnsystem
	.def("__init__",
	[](PNSystem &m, const std::string& stencil_name, Domain &domain, bool neumannBC)
	{
		PNSystem::Stencil stencil;

		if( stencil_name == "noop" )
			stencil = PNSystem::Stencil::noop();
		else
			stencil = PNSystem::findStencil(stencil_name);
		new (&m) PNSystem(stencil, domain, neumannBC);
	})
	.def("getNumCoefficients", &PNSystem::getNumCoefficients )
	.def("getResolution", &PNSystem::getResolution )
	.def("getNumVoxels", &PNSystem::getNumVoxels )
	.def("getOrder", &PNSystem::getOrder )
	.def("build", &PNSystem::build )
	.def("setField", &PNSystem::setField )
	.def("get_A_real", &PNSystem::get_A_real )
	.def("get_b_real", &PNSystem::get_b_real )
	.def("get_A_real_test", &PNSystem::get_A_real_test )
	.def("get_b_real_test", &PNSystem::get_b_real_test )
	.def("setDebugVoxel", &PNSystem::setDebugVoxel )
	.def("get_debug", &PNSystem::get_debug )

	//.def("computeGroundtruth", &PNSystem::computeGroundtruth )
	.def("getVoxelInfo", &PNSystem::getVoxelInfo )
	;

	// Domain ============================================================
	py::class_<Domain> class_domain(m, "Domain");
	class_domain
	.def("__init__",
	[](Domain &m, const Eigen::Matrix<double, 3, 1>& size, const Eigen::Matrix<int, 3, 1>& resolution, const Eigen::Matrix<double, 3, 1>& offset)
	{
		new (&m) Domain( size,
						 resolution,
						 offset);
	})
	.def("getResolution",
		[](Domain &m)
		{
			return m.getResolution();
		})
	.def("setResolution",
		[](Domain &m, const Eigen::Matrix<int, 3, 1>& resolution)
		{
			m.setResolution(resolution);
		})
	.def("getVoxelSize",
		[](Domain &m)
		{
			return m.getVoxelSize();
		})
	.def("numVoxels",
		[](Domain &m)
		{
			return m.numVoxels();
		})
	.def("worldToVoxel",
		[](Domain &m, const Eigen::Matrix<double, 3, 1>& pWS)
		{
			return m.worldToVoxel(pWS);
		})
	.def("voxelToWorld",
		[](Domain &m, const Eigen::Matrix<double, 3, 1>& pVS)
		{
			return m.voxelToWorld(pVS);
		})
	.def("getBoundMax",
		[](Domain &m)
		{
			return m.getBoundMax();
		})
	.def("getBoundMin",
		[](Domain &m)
		{
			return m.getBoundMin();
		})
	;


	// Field ============================================================
	py::class_<Field, Field::Ptr> class_field(m, "Field");
	class_field
	.def("__call__",
	[](Field &m, const Eigen::Matrix<double, 3, 1>& pWS)
	{
		return m.eval(pWS);
	})
	;

	// VoxelGrid ============================================================
	py::class_<VoxelGridField, VoxelGridField::Ptr> class_VoxelGridField(m, "VoxelGridField", class_field);
	class_VoxelGridField
	.def("__init__",
	[](VoxelGridField &m, py::array b, Domain& domain, const Eigen::Matrix<double, 3, 1>& offset)
	{
		typedef Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic> Strides;

		// Request a buffer descriptor from Python
		py::buffer_info info = b.request();

		// Some sanity checks ...
		if (info.format != py::format_descriptor<std::complex<double>>::format())
			throw std::runtime_error("Incompatible format: expected a complex array!");

		if (info.ndim != 3)
			throw std::runtime_error("Incompatible buffer dimension!");

		int res_x = int(info.shape[0]);
		int res_y = int(info.shape[1]);
		int res_z = int(info.shape[2]);

		auto data = static_cast<std::complex<double> *>(info.ptr);
		new (&m) VoxelGridField( data, domain, offset );
	})
	.def("test", &VoxelGridField::test)
	.def("getSlice", &VoxelGridField::getSlice)
	;

	// Constant ============================================================
	py::class_<Constant, Constant::Ptr> class_constant(m, "Constant", class_field);
	class_constant
	.def("__init__",
	[](Constant &m, std::complex<double> value)
	{
		new (&m) Constant(value);
	});

	// RadianceField ============================================================
	py::class_<RadianceField, RadianceField::Ptr> class_radiancefield(m, "RadianceField");
	class_radiancefield
	.def("__call__",
	[](RadianceField &m, const Eigen::Matrix<double, 3, 1>& pWS, const Eigen::Matrix<double, 3, 1>& omega)
	{
		return m.eval(pWS, omega);
	})
	;

	// SHEXP ============================================================
	py::class_<SHEXP, SHEXP::Ptr> class_shexp(m, "SHEXP", class_radiancefield);
	class_shexp
	.def("__init__",
	[](SHEXP &m, int order)
	{
		new (&m) SHEXP(order);
	})
	.def("setCoefficientField", &SHEXP::setCoefficientField)
	.def("getCoefficientField", &SHEXP::getCoefficientField)
	;
}
