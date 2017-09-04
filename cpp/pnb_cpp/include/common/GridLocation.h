#pragma once




#include <math/common.h>
#include <math/vector.h>
#include <math/ray.h>



#include <common/Domain.h>



inline div_t py_divmod( int a, int b )
{
	div_t result;
	result.quot = (a - (((a % b) + b) % b)) / b;
	result.rem = ((a % b) + b) % b;
	return result;
}


struct GridLocation
{
	GridLocation( const Domain& domain, const V2i& voxel, const V2i& offset ):
	    m_domain(domain)
	{
		auto d0 = py_divmod(offset[0], 2);
		auto d1 = py_divmod(offset[1], 2);

		m_voxel = V2i(voxel[0] + d0.quot,voxel[1] + d1.quot );
		m_offset = V2i( d0.rem, d1.rem );
		m_pWS = domain.voxelToWorld( P2d(m_voxel[0] + m_offset[0]*0.5, m_voxel[1] + m_offset[1]*0.5) );
	}

	V2i getOffset()const
	{
		return m_offset;
	}
	V2i getVoxel()const
	{
		return m_voxel;
	}
	V2d getPWS()const
	{
		return m_pWS;
	}
	GridLocation getShiftedLocation(const V2i& offset)const
	{
		return GridLocation(m_domain, m_voxel, m_offset+offset);
	}

private:
	Domain m_domain;
	V2i m_offset;
	V2i m_voxel;
	V2d m_pWS;
};
