#pragma once

#include <memory>

#include <math/common.h>
#include <math/vector.h>
#include <math/bbox.h>








// datastructure for working with the solution from a PNSolver
// basically a voxelgrid with a constant number of elements per voxel
// TODO: switch everything to real only
struct PNSolution
{
	typedef std::shared_ptr<PNSolution> Ptr;

	PNSolution(const std::string& filename);
	PNSolution(int order, const V3i& resolution, const Box3d& bound, const std::complex<double> *data);

	// evaluate <- evalute SH function for given direction at given worldspace position (does interpolation of SH coefficients)
	double eval(const P3d &pWS, const V3d &direction)const;
	std::complex<double> evalCoefficient(const P3d &pWS, int coeff_index)const;



	// sampling tools
	// pdf <-
	// sample <- importance sampling


	// io
	void save( const std::string& filename );


	// space transforms
	P3d voxelToLocal(const P3d& pVS)const;
	P3d localToVoxel(const P3d& pLS)const;
	P3d localToWorld( const P3d& pLS )const;
	P3d worldToLocal( const P3d& pWS )const;
	P3d voxelToWorld( const P3d& pVS)const;
	P3d worldToVoxel( const P3d& pWS)const;
	bool contains_voxel( const P3i& voxel )const;

	// info
	int getOrder()const;
	int getNumCoeffs()const;
	V3i getResolution()const;
	P3d getBoundMin()const;
	P3d getBoundMax()const;
	std::complex<double>* data();
private:

	int getIndex( const V3i& voxel_coord ) const;
	int m_order; // SH order
	int m_numCoeffs;
	V3i m_resolution;
	Box3d m_bound;
	V3d m_extend;
	std::vector<std::complex<double>> m_data; // number of voxels*number of coefficients
};

