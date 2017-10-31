#pragma once

#include <math/common.h>
#include <math/bbox.h>
#include <math/vector.h>
#include <math/ray.h>






// The Domain datastructure is used to define a voxelgrid over a rectangular (axis aligned)
// section of the worldspace. It provides mappings between world and voxel space.
struct Domain
{
	Domain( const V3d& size, V3i resolution, const P3d& offset ):
	    m_resolution(resolution),
	    m_bound(offset, offset+size),
		m_voxelsize( size[0]/resolution[0], size[1]/resolution[1], size[2]/resolution[2] )
	{
		m_extend = m_bound.getExtents();
	}

	P3i getResolution()const
	{
		return m_resolution;
	}

	void setResolution(const V3i& resolution)
	{
		m_resolution = resolution;
		V3d size = m_bound.getExtents();
		m_voxelsize = V3d( size[0]/resolution[0], size[1]/resolution[1], size[2]/resolution[2] );
	}

	P3d getVoxelSize()const
	{
		return m_voxelsize;
	}

	int numVoxels()const
	{
		return m_resolution[0]*m_resolution[1]*m_resolution[2];
	}

	P3d voxelToLocal(const P3d& pVS)const
	{
		return P3d(pVS[0]/m_resolution[0], pVS[1]/m_resolution[1], pVS[2]/m_resolution[2]);
	}
	P3d localToVoxel(const P3d& pLS)const
	{
		return P3d(pLS[0]*m_resolution[0], pLS[1]*m_resolution[1], pLS[2]*m_resolution[2]);
	}
	P3d localToWorld( const P3d& pLS )const
	{
		return P3d(pLS[0]*m_extend[0] + m_bound.min[0], pLS[1]*m_extend[1] + m_bound.min[1], pLS[2]*m_extend[2] + m_bound.min[2]);
	}

	P3d worldToLocal( const P3d& pWS )const
	{
		return P3d((pWS[0]-m_bound.min[0])/m_extend[0], (pWS[1]-m_bound.min[1])/m_extend[1], (pWS[2]-m_bound.min[2])/m_extend[2] );
	}

	P3d voxelToWorld( const P3d& pVS)const
	{
		return localToWorld(voxelToLocal(pVS));
	}


	P3d worldToVoxel( const P3d& pWS)const
	{
		return localToVoxel(worldToLocal(pWS));
	}

	bool contains_voxel( const P3i& voxel )const
	{
		if( (voxel[0] < 0)||(voxel[0] >= m_resolution[0])||
			(voxel[1] < 0)||(voxel[1] >= m_resolution[1])||
			(voxel[2] < 0)||(voxel[2] >= m_resolution[2]))
			return false;
		return true;
	}

	Box3d getBound()const
	{
		return m_bound;
	}

	P3d getBoundMax()const
	{
		return m_bound.max;
	}

	P3d getBoundMin()const
	{
		return m_bound.min;
	}

	V3d getSize()const
	{
		return m_bound.getExtents();
	}

	Domain downsample()const
	{
		V3i res_fine = getResolution();

		bool is2D = res_fine[2] == 1;

		// for multigrid, we require the resolution to be even,
		// so that we can do restriction and interpolation straigh forwardly
		if( (res_fine[0]%2!=0)||(res_fine[1]%2!=0)||(!is2D && (res_fine[2]%2!=0)))
			throw std::runtime_error("Domain::downsample currently requires even resolution");

		V3i res_coarse( res_fine[0]/2, res_fine[1]/2, is2D ? 1:res_fine[2]/2 );

		return Domain(  getBound().getExtents(),
						res_coarse,
						getBound().min );
	}

private:
	V3i m_resolution;
	V3d m_voxelsize;
	Box3d m_bound;
	V3d m_extend;
};
