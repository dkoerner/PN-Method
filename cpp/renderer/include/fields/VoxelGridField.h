#pragma once

#include <complex>
#include <iostream>
#include <algorithm>

#include <math/common.h>
#include <math/vector.h>
#include <math/ray.h>

#include <util/voxelgrid.h>

#include <field/Field.h>

#include <Eigen/Dense>


struct VoxelGridField : public Field
{
	typedef std::shared_ptr<VoxelGridField> Ptr;

	VoxelGridField( std::complex<double>* data, const Domain& domain, V3d offset ):
		m_domain(domain),
		m_voxelgrid()
	{
		V3i res= m_domain.getResolution();
		m_voxelgrid.resize(res);
		if(data)
			memcpy( m_voxelgrid.getRawPointer(), data, domain.numVoxels()*sizeof(std::complex<double>) );

		m_voxelgrid.m_sampleLocation = offset.cast<float>();

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

	virtual std::complex<double> eval( const P3d& pWS )const override
	{
		P3d pVS = m_domain.worldToVoxel(pWS);
		return m_voxelgrid.evaluate(pVS);
	}

	/*
	double getMax()const
	{
		return m_max_value;
	}
	*/

	V3d getOffset()const
	{
		return m_voxelgrid.m_sampleLocation.cast<double>();
	}

	Eigen::MatrixXd getSlice( int k )
	{
		V3i res = m_domain.getResolution();
		Eigen::MatrixXd m = Eigen::MatrixXd::Zero(res[0],res[1]);

		for(int i=0;i<res[0];++i)
			for(int j=0;j<res[1];++j)
			{
				m(i, j) = m_voxelgrid.sample(i, j, k).real();
			}

		return m;
	}

	virtual Field::Ptr createRestricted()const override
	{
		V3i res_fine = m_domain.getResolution();
		bool is2D = res_fine[2] == 1;

		// for restriction, we require the resolution to be even,
		// so that we can do it straight forwardly
		if( (res_fine[0]%2!=0)||(res_fine[1]%2!=0)||(!is2D && (res_fine[2]%2!=0)))
			throw std::runtime_error("VoxelGridField::createRestricted restriction currently requires even resolution");

		V3i res_coarse( res_fine[0]/2, res_fine[1]/2, is2D ? 1:res_fine[2]/2 );

		Domain domain_coarse( m_domain.getBound().getExtents(), res_coarse, m_domain.getBound().min );

		Ptr result = std::make_shared<VoxelGridField>( (std::complex<double>*) 0, domain_coarse, getOffset() );

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
					result->m_voxelgrid.lvalue( i_coarse, j_coarse, 0 ) = v/4.0;
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
						result->m_voxelgrid.lvalue( i_coarse, j_coarse, k_coarse ) = v/8.0;
					}
		}


		return result;
	}

private:
	Domain m_domain;
	VoxelGrid<std::complex<double>> m_voxelgrid;
	//double m_max_value;
};


