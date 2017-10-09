#pragma once

#include <complex>
#include <iostream>
#include <algorithm>

#include <math/common.h>
#include <math/vector.h>
#include <math/ray.h>

#include <util/voxelgrid.h>

#include <field/Field.h>




struct VoxelGridField : public Field
{
	typedef std::shared_ptr<VoxelGridField> Ptr;

	VoxelGridField( std::complex<double>* data, const Domain& domain, V3d offset ):
		m_domain(domain),
		m_voxelgrid()
	{
		m_step_x = V3d(m_domain.getVoxelSize()[0]*0.5, 0.0, 0.0);
		m_step_y = V3d(0.0, m_domain.getVoxelSize()[1]*0.5, 0.0);
		m_step_z = V3d(0.0, 0.0, m_domain.getVoxelSize()[2]*0.5);

		V3i res= m_domain.getResolution();
		m_voxelgrid.resize(res);
		memcpy( m_voxelgrid.getRawPointer(), data, domain.numVoxels()*sizeof(std::complex<double>) );
		m_voxelgrid.m_sampleLocation = offset.cast<float>();

		//double max_value = *std::max_element(std::begin(data), std::end(data));
		m_max_value = 0.0;
		for( int i=0;i<domain.numVoxels();++i )
		{
			double value = data[i].real();
			if( value > m_max_value )
				m_max_value = value;
		}
	}

	virtual std::complex<double> eval( const P3d& pWS )const override
	{
		P3d pVS = m_domain.worldToVoxel(pWS);
		return m_voxelgrid.evaluate(pVS);
	}

	double getMax()const
	{
		return m_max_value;
	}


private:
	Domain m_domain;
	V3d m_step_x;
	V3d m_step_y;
	V3d m_step_z;
	VoxelGrid<std::complex<double>> m_voxelgrid;
	double m_max_value;
};


