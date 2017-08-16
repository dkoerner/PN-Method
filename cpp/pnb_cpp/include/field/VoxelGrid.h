#pragma once

#include <complex>
#include <iostream>

#include <math/common.h>
#include <math/vector.h>
#include <math/ray.h>

#include <field/Field.h>





struct VoxelGrid : public Field
{
	typedef std::shared_ptr<VoxelGrid> Ptr;

	VoxelGrid( std::complex<double>* data, const Domain& domain, V2d offset ):
		m_data(data),
		m_domain(domain),
		m_offset(offset)
	{
		m_data2 = std::vector<std::complex<double>>(data, data+domain.numVoxels());
		m_step_x = V2d(m_domain.voxelSize()[0]*0.5, 0.0);
		m_step_y = V2d(0.0, m_domain.voxelSize()[1]*0.5);
	}

	std::complex<double> sample( const V2i& voxel )const
	{
		int index = voxel[0]*m_domain.resolution()[1] + voxel[1];
		return m_data2[index];
	}

	virtual std::complex<double> eval( const P2d& pWS )const override
	{
		P2d pVS = m_domain.worldToVoxel(pWS)-m_offset;
		P2i v0;
		v0[0] = (int)floor(pVS[0]);
		v0[1] = (int)floor(pVS[1]);
		v0[2] = (int)floor(pVS[2]);
		P2i v1(v0[0]+1, v0[1]+1);

		double tx = pVS.x() - v0[0];
		double ty = pVS.y() - v0[1];

		// clamp to domain boundaries
		V2i res = m_domain.resolution();
		if( v0[0] < 0 )
			v0[0] = 0;
		else
		if( v0[0] >= res[0] )
			v0[0] = res[0]-1;

		if( v0[1] < 0 )
			v0[1] = 0;
		else
		if( v0[1] >= res[1] )
			v0[1] = res[1]-1;

		if( v1[0] < 0 )
			v1[0] = 0;
		else
		if( v1[0] >= res[0] )
			v1[0] = res[0]-1;

		if( v1[1] < 0 )
			v1[1] = 0;
		else
		if( v1[1] >= res[1] )
			v1[1] = res[1]-1;

		//return self.sample(v0[0], v0[1])

		std::complex<double> result = 0.0;
		result += sample(V2i(v0[0], v0[1]))*(1.0-tx)*(1.0-ty);
		result += sample(V2i(v1[0], v0[1]))*tx*(1.0-ty);
		result += sample(V2i(v0[0], v1[1]))*(1.0-tx)*ty;
		result += sample(V2i(v1[0], v1[1]))*tx*ty;
		return result;
	}
	virtual std::complex<double> dx(const P2d& pWS)const override
	{
		std::complex<double> a = eval(pWS-m_step_x);
		std::complex<double> b = eval(pWS+m_step_x);
		return (b-a)/m_domain.voxelSize()[0];
	}
	virtual std::complex<double> dxdx(const P2d& pWS)const override
	{
		std::complex<double> a = dx(pWS-m_step_x);
		std::complex<double> b = dx(pWS+m_step_x);
		return (b-a)/m_domain.voxelSize()[0];
	}
	virtual std::complex<double> dxdy(const P2d& pWS)const override
	{
		std::complex<double> a = dx(pWS-m_step_y);
		std::complex<double> b = dx(pWS+m_step_y);
		return (b-a)/m_domain.voxelSize()[1];
	}
	virtual std::complex<double> dy(const P2d& pWS)const override
	{
		std::complex<double> a = eval(pWS-m_step_y);
		std::complex<double> b = eval(pWS+m_step_y);
		return (b-a)/m_domain.voxelSize()[1];
	}
	virtual std::complex<double> dydy(const P2d& pWS)const override
	{
		std::complex<double> a = dy(pWS-m_step_y);
		std::complex<double> b = dy(pWS+m_step_y);
		return (b-a)/m_domain.voxelSize()[1];
	}
	virtual std::complex<double> dydx(const P2d& pWS)const override
	{
		return dxdy(pWS);
	}
	virtual std::complex<double> dz(const P2d& pWS)const override
	{
		return 0.0;
	}

	void test()const
	{
		for(int s=0;s<10;++s)
			for(int t=0;t<10;++t)
			{
				std::cout << sample(V2i(s,t)) << std::endl;
			}

	}

private:
	Domain m_domain;
	V2d m_offset;
	V2d m_step_x;
	V2d m_step_y;
	std::complex<double>* m_data;
	std::vector<std::complex<double>> m_data2;
};


