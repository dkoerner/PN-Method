#include <util/voxelgrid.h>



// needs to go to VoxelGrid.cpp for gcc compatibility...
template<> const int VoxelGrid<float>::m_dataType = 1;
template<> const int VoxelGrid<Vector3f>::m_dataType = 2;
template<> const int VoxelGrid<Vector4f>::m_dataType = 3;
template<> const int VoxelGrid<double>::m_dataType = 4;
template<> const int VoxelGrid<Vector3d>::m_dataType = 5;



