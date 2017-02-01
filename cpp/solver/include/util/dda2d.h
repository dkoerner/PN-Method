#pragma once


#include <math/common.h>
#include <math/vector.h>
#include <math/transform.h>
#include <math/bbox.h>







struct DDA2D
{
	struct PixelStep
	{
		double min, max;
		V2i coord;
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


	bool initialize( const V2i& resolution, const Transformd& localToWorld, const Ray2d& rayWS, bool debug = false);
	bool initialize( const V2i& resolution, const P2d& p0, const P2d& p1, bool debug = false);

	bool step( PixelStep& info, bool debug = false );
	V2i res;




	double m_tstart; // worldspace position along rayWS of intersection with the volume
	double m_tend; // worldspace position along rayWS of end of intersection with volume
	double m_lastmax; // t for maximum t at last voxel step
	double ratio;


	V2i vox; // current voxel
	V2i end; // first invalid voxel coordinates according to the current ray
	V2i m_step; // integers of -1 or 1, indicating wether voxel indices are incremented or decremented (determined by sign of ray direction)
	V2d tmax; // distance along the ray until we hit the next voxel boundary in each dimension
	V2d delta; // how far along a ray with must move (in units of t) for the horizontal component of such a movement to equal the width of a voxel
	int axis; // axis along which we did the last step
};

// pDC is where line is drawn. Can be memory device context.
// (X0,Y0) is start point of line.
// (X1, Y1) is end point of line.
// BaseColor is intensity of line. Pass 0 for black line.
// NumLevels is number of gray scale levels. Pass 256.
// IntensityBits denotes bits used to represent color component. Pass 8.
typedef std::function<void(short, short, short)> DrawPixelFunction;
void draw_wu_line( short X0, short Y0, short X1, short Y1, short BaseColor, short NumLevels, unsigned short IntensityBits, DrawPixelFunction DrawPixel);
