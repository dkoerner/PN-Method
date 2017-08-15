#pragma once

#include <complex>
#include <iostream>

#include <math/common.h>
#include <math/vector.h>
#include <math/ray.h>

#include <common/Domain.h>
#include <common/GridLocation.h>






struct CoefficientGrid
{
	CoefficientGrid( std::complex<double>* data, const Domain& domain, V2i offset ):
		m_data(data),
		m_domain(domain),
		m_offset(offset)
	{
		m_data2 = std::vector<std::complex<double>>(data, data+domain.numVoxels());
	}

	std::complex<double> sample( const V2i& voxel )const
	{
		int index = voxel[0]*m_domain.resolution()[1] + voxel[1];
		return m_data2[index];
	}

	std::complex<double> eval( GridLocation& loc )const
	{
		V2i location_offset = loc.getOffset();
		V2i voxel = loc.getVoxel();

		if( (location_offset[0] == m_offset[0])&&
			(location_offset[1] == m_offset[1]))
		{
			return sample(voxel);
		}else
		if(location_offset[0] == m_offset[0])
		{
			std::complex<double> l = sample(loc.getShiftedLocation(V2i(0, 1)).getVoxel());
			std::complex<double> r = sample(loc.getShiftedLocation(V2i(0, -1)).getVoxel());
			return 0.5*(l+r);
		}else
		if(location_offset[1] == m_offset[1])
		{
			std::complex<double> l = sample(loc.getShiftedLocation(V2i(1, 0)).getVoxel());
			std::complex<double> r = sample(loc.getShiftedLocation(V2i(-1,0)).getVoxel());
			return 0.5*(l+r);
		}else
		{
			const std::vector<V2i> offsets = {V2i(-1,-1), V2i(-1,1), V2i(1,-1), V2i(1,1)};
			std::complex<double> result = 0.0;
			for( auto offset:offsets )
			{
				result += sample(loc.getShiftedLocation(offset).getVoxel());
			}
			return result/double(offsets.size());
		}

		return std::complex<double>(0.0,0.0);
	}
	std::complex<double> dx(const GridLocation& loc)
	{
		V2i step(1,0);
		std::complex<double> a = eval(loc.getShiftedLocation(-step));
		std::complex<double> b = eval(loc.getShiftedLocation(step));
		return (b-a)/m_domain.voxelSize()[0];
	}
	std::complex<double> dxdx(const GridLocation& loc)
	{
		V2i step(1,0);
		std::complex<double> a = dx(loc.getShiftedLocation(-step));
		std::complex<double> b = dx(loc.getShiftedLocation(step));
		return (b-a)/m_domain.voxelSize()[0];
	}
	std::complex<double> dxdy(const GridLocation& loc)
	{
		V2i step(0,1);
		std::complex<double> a = dx(loc.getShiftedLocation(-step));
		std::complex<double> b = dx(loc.getShiftedLocation(step));
		return (b-a)/m_domain.voxelSize()[1];
	}
	std::complex<double> dy(const GridLocation& loc)
	{
		V2i step(0,1);
		std::complex<double> a = eval(loc.getShiftedLocation(-step));
		std::complex<double> b = eval(loc.getShiftedLocation(step));
		return (b-a)/m_domain.voxelSize()[1];
	}
	std::complex<double> dydy(const GridLocation& loc)
	{
		V2i step(0,1);
		std::complex<double> a = dy(loc.getShiftedLocation(-step));
		std::complex<double> b = dy(loc.getShiftedLocation(step));
		return (b-a)/m_domain.voxelSize()[1];
	}
	std::complex<double> dydx(const GridLocation& loc)
	{
		return dxdy(loc);
	}
	std::complex<double> dz(const GridLocation& loc)
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
	V2i m_offset;
	std::complex<double>* m_data;
	std::vector<std::complex<double>> m_data2;
};


