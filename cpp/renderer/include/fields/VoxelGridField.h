#pragma once

#include <complex>
#include <iostream>
#include <algorithm>

#include <math/common.h>
#include <math/vector.h>
#include <math/ray.h>

#include <util/voxelgrid.h>

#include <Field.h>

#include <Eigen/Dense>


template<typename T>
struct VoxelGridField : public Field<T>
{
	typedef std::shared_ptr<VoxelGridField<T>> Ptr;

	VoxelGridField( const V3i& resolution, T* data ):
		m_resolution(resolution),
		m_voxelgrid()
	{
		m_voxelgrid.resize(m_resolution);
		if(data)
			memcpy( m_voxelgrid.getRawPointer(), data, m_resolution[0]*m_resolution[1]*m_resolution[2]*sizeof(T) );

		/*
		//double max_value = *std::max_element(std::begin(data), std::end(data));
		m_max_value = 0.0;
		for( int i=0;i<domain.numVoxels();++i )
		{
			double value = data[i].real();
			if( value > m_max_value )
				m_max_value = value;
		}
		*/
	}

	//virtual std::pair<T, T> getValueRange()const
	virtual T getMaxValue()const
	{
		throw std::runtime_error("getMaxValue not implemented for VoxelGridField");
		return T(0.0);
		//return std::make_pair(T(0.0), T(0.0));
	}

	virtual T eval( const P3d& pLS )const override
	{
		return m_voxelgrid.evaluate(localToVoxel(pLS));
	}

	/*
	double getMax()const
	{
		return m_max_value;
	}
	*/

	V3i getResolution()const
	{
		return m_resolution;
	}

	VoxelGrid<T>& getVoxelGrid()
	{
		return m_voxelgrid;
	}

	T* getData()
	{
		return m_voxelgrid.getRawPointer();
	}


	P3d localToVoxel(const P3d& pLS)const
	{
		return P3d(pLS[0]*m_resolution[0], pLS[1]*m_resolution[1], pLS[2]*m_resolution[2]);
	}

	P3d voxelToLocal(const P3d& pVS)const
	{
		return P3d(pVS[0]/m_resolution[0], pVS[1]/m_resolution[1], pVS[2]/m_resolution[2]);
	}

private:
	V3i m_resolution;
	VoxelGrid<T> m_voxelgrid;
	//double m_max_value;
};
typedef VoxelGridField<double> VoxelGridFieldd;
typedef VoxelGridField<V3d> VoxelGridField3d;
typedef VoxelGridField<std::complex<double>> VoxelGridFieldcd;


