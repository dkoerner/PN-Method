#pragma once

#include <complex>
#include <iostream>
#include <algorithm>

#include <math/common.h>
#include <math/vector.h>
#include <math/ray.h>

#include <util/voxelgrid.h>

#include <fields/Field.h>

#include <Eigen/Dense>


template<typename T>
struct VoxelGridField : public Field<T>
{
	typedef std::shared_ptr<VoxelGridField<T>> Ptr;

	VoxelGridField( const V3i& resolution, T* data ):
		m_voxelgrid()
	{
		m_voxelgrid.resize(resolution);
		if(data)
			memcpy( m_voxelgrid.getRawPointer(), data, resolution[0]*resolution[1]*resolution[2]*sizeof(T) );

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

	VoxelGridField(const std::string& filename):
		m_voxelgrid(filename)
	{
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
		//std::cout << "VoxelGridField::eval pLS=" << pLS.toString() << std::endl;
		//std::cout << "VoxelGridField::eval pVS=" << localToVoxel(pLS).toString() << std::endl;

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
		return m_voxelgrid.getResolution();
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
		return m_voxelgrid.localToVoxel(pLS);
	}

	P3d voxelToLocal(const P3d& pVS)const
	{
		return m_voxelgrid.voxelToLocal(pLS);
	}

	void save( const std::string& filename )
	{
		m_voxelgrid.save(filename);
	}

private:
	VoxelGrid<T> m_voxelgrid;
	//double m_max_value;
};
typedef VoxelGridField<double> VoxelGridFieldd;
typedef VoxelGridField<V3d> VoxelGridField3d;
typedef VoxelGridField<std::complex<double>> VoxelGridFieldcd;


