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

	// returns a coarser mimap version of this field
	// this is used by the multigrid solver
	virtual Field::Ptr downsample()const
	{
		V3i res_fine = m_voxelgrid.getResolution();
		bool is2D = res_fine[2] == 1;

		// for downsampling, we require the resolution to be even,
		// so that we can do it straight forwardly
		if( (res_fine[0]%2!=0)||(res_fine[1]%2!=0)||(!is2D && (res_fine[2]%2!=0)))
			throw std::runtime_error("VoxelGridField::downsample currently requires even resolution");

		V3i res_coarse( res_fine[0]/2, res_fine[1]/2, is2D ? 1:res_fine[2]/2 );

		Ptr result = std::make_shared<VoxelGridField<T>>( res_coarse, (T*)0 );

		// create coarse from fine voxels
		if(is2D)
		{
			for(int i_coarse=0;i_coarse<res_coarse[0];++i_coarse)
				for(int j_coarse=0;j_coarse<res_coarse[1];++j_coarse)
				{
					auto v = m_voxelgrid.sample(i_coarse*2+0, j_coarse*2+0, 0);
					v+= m_voxelgrid.sample(i_coarse*2+1, j_coarse*2+0, 0);
					v+= m_voxelgrid.sample(i_coarse*2+1, j_coarse*2+1, 0);
					v+= m_voxelgrid.sample(i_coarse*2+0, j_coarse*2+1, 0);
					result->getVoxelGrid().lvalue( i_coarse, j_coarse, 0 ) = v/4.0;
				}
		}else
		{
			for(int i_coarse=0;i_coarse<res_coarse[0];++i_coarse)
				for(int j_coarse=0;j_coarse<res_coarse[1];++j_coarse)
					for(int k_coarse=0;k_coarse<res_coarse[2];++k_coarse)
					{
						auto v = m_voxelgrid.sample(i_coarse*2+0, j_coarse*2+0, k_coarse*2+0);
						v+= m_voxelgrid.sample(i_coarse*2+1, j_coarse*2+0, k_coarse*2+0);
						v+= m_voxelgrid.sample(i_coarse*2+1, j_coarse*2+1, k_coarse*2+0);
						v+= m_voxelgrid.sample(i_coarse*2+0, j_coarse*2+1, k_coarse*2+0);

						v+= m_voxelgrid.sample(i_coarse*2+0, j_coarse*2+0, k_coarse*2+1);
						v+= m_voxelgrid.sample(i_coarse*2+1, j_coarse*2+0, k_coarse*2+1);
						v+= m_voxelgrid.sample(i_coarse*2+1, j_coarse*2+1, k_coarse*2+1);
						v+= m_voxelgrid.sample(i_coarse*2+0, j_coarse*2+1, k_coarse*2+1);
						result->getVoxelGrid().lvalue( i_coarse, j_coarse, k_coarse ) = v/8.0;
					}
		}

		return result;
	}

private:
	VoxelGrid<T> m_voxelgrid;
	//double m_max_value;
};
typedef VoxelGridField<double> VoxelGridFieldd;
typedef VoxelGridField<V3d> VoxelGridField3d;
typedef VoxelGridField<std::complex<double>> VoxelGridFieldcd;


