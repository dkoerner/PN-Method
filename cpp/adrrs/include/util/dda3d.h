#pragma once


#include <math/common.h>
#include <math/vector.h>
#include <math/transform.h>
#include <math/bbox.h>







struct DDA3D
{
	struct VoxelStep
	{
		double min, max;
		V3i coord;
	};


	double diff_distance( double x, double ds )
	{
		double s;

		if(x<0.0)
			s = x-int(-1+x);
		else
			s = x-int(x);

		if(ds > 0.0)
			return (1.0-s)/ds;
		else
			return -s/ds;
		return 0.0;
	}


	bool initialize( const V3i& resolution, const Transformd& localToWorld, const Ray3d& rayWS, bool debug = false);

	bool step( VoxelStep& info, bool debug = false );
	V3i res;




	double m_tstart; // worldspace position along rayWS of intersection with the volume
	double m_tend; // worldspace position along rayWS of end of intersection with volume
	double m_lastmax; // t for maximum t at last voxel step
	double ratio;


	V3i vox; // current voxel
	V3i end; // first invalid voxel coordinates according to the current ray
	V3i m_step; // integers of -1 or 1, indicating wether voxel indices are incremented or decremented (determined by sign of ray direction)
	V3d tmax; // distance along the ray until we hit the next voxel boundary in each dimension
	V3d delta; // how far along a ray with must move (in units of t) for the horizontal component of such a movement to equal the width of a voxel
	int axis; // axis along which we did the last step
};
