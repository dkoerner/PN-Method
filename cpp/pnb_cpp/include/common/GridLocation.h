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

		//std::cout << "GridLocation: div(" << offset[0]  << ",2)=" << d0.quot << " " << d0.rem << std::endl;

		//std::cout << offset[0] << " " << offset[1] << std::endl;
		//std::cout << "quotient: " <<d0.quot << " " << d1.quot << std::endl;
		//std::cout << "remainder: "<<d0.rem << " " << d1.rem << std::endl;

		m_voxel = V2i(voxel[0] + d0.quot,voxel[1] + d1.quot );
		m_offset = V2i( d0.rem, d1.rem );

		// TODO: m_pWS
		//#self.pWS = self.domain.bound_min + np.multiply(np.array([self.voxel_i+self.offset[0]*0.5, self.voxel_j+self.offset[1]*0.5]), self.domain.voxelsize)
		//#self.pWS = self.domain.bound_min + np.multiply(self.voxel+self.offset*0.5, self.domain.voxelsize)
		//x = (self.voxel[0]+self.offset[0]*0.5)*self.domain.voxelsize[0]
		//y = (self.voxel[1]+self.offset[1]*0.5)*self.domain.voxelsize[1]
		//self.pWS = np.array([x, y])

	}

	V2i getOffset()const
	{
		return m_offset;
	}
	V2i getVoxel()const
	{
		return m_voxel;
	}
	GridLocation getShiftedLocation(const V2i& offset)const
	{
		return GridLocation(m_domain, m_voxel, m_offset+offset);
	}

private:
	Domain m_domain;
	V2i m_offset;
	V2i m_voxel;
	//V2d m_pWS;
};
