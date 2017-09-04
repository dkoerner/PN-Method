#pragma once

#include <math/common.h>
#include <math/bbox.h>
#include <math/vector.h>
#include <math/ray.h>






// The Domain datastructure is used to define a voxelgrid over a rectangular (axis aligned)
// section of the worldspace. It provides mappings between world and voxel space.
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

	void setResolution(const V2i& resolution)
	{
		m_resolution = resolution;
		V2d size = m_bound.getExtents();
		m_voxelsize( size[0]/resolution[0], size[1]/resolution[1] );
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
	P2d localToVoxel(const P2d& pLS)const
	{
		return P2d(pLS[0]*m_resolution[0], pLS[1]*m_resolution[1]);
	}

	P2d localToWorld( const P2d& pLS )const
	{
		return P2d(pLS[0]*m_extend[0] + m_bound.min[0], pLS[1]*m_extend[1] + m_bound.min[1]);
	}

	P2d worldToLocal( const P2d& pWS )const
	{
		return P2d((pWS[0]-m_bound.min[0])/m_extend[0] , (pWS[1]-m_bound.min[1])/m_extend[1] );
	}

	P2d voxelToWorld( const P2d& pVS)const
	{
		return localToWorld(voxelToLocal(pVS));
	}


	P2d worldToVoxel( const P2d& pWS)const
	{
		return localToVoxel(worldToLocal(pWS));
	}

	bool contains_voxel( const P2i& voxel )const
	{
		if( (voxel[0] < 0)||(voxel[0] >= m_resolution[0])||
			(voxel[1] < 0)||(voxel[1] >= m_resolution[1]))
			return false;
		return true;
	}

private:
	V2i m_resolution;
	V2d m_voxelsize;
	Box2d m_bound;
	V2d m_extend;
};
