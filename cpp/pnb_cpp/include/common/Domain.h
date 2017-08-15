#pragma once

#include <math/common.h>
#include <math/bbox.h>
#include <math/vector.h>
#include <math/ray.h>







struct Domain
{
	Domain( const V2d& size, V2i resolution, const P2d& offset ):
	    m_resolution(resolution),
	    m_bound(offset, offset+size),
	    m_voxelsize( size[0]/resolution[0], size[1]/resolution[1] )
	{
		m_extend = m_bound.getExtents();
	}

	P2i resolution()const
	{
		return m_resolution;
	}
	P2d voxelSize()const
	{
		return m_voxelsize;
	}

	int numVoxels()const
	{
		return m_resolution[0]*m_resolution[1];
	}

	P2d voxelToLocal(const P2d& pVS)const
	{
		return P2d(pVS[0]/m_resolution[0], pVS[1]/m_resolution[1]);
	}

	P2d localToWorld( const P2d& pLS )const
	{
		return P2d(pLS[0]*m_extend[0] + m_bound.min[0], pLS[1]*m_extend[1] + m_bound.min[1]);
	}

	P2d voxelToWorld( const P2d& pVS)const
	{
		return localToWorld(voxelToLocal(pVS));
	}

private:
	V2i m_resolution;
	V2d m_voxelsize;
	Box2d m_bound;
	V2d m_extend;
};
