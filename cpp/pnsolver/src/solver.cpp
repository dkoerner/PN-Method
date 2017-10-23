#include <solver.h>





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









