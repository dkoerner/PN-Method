#include <util/dda2d.h>

bool DDA2D::initialize( const V2i& resolution, const Transformd& localToWorld, const Ray2d& rayWS, bool debug)
{
	res = resolution;

	// intersect with bounding box
	//Ray2d rayLS = localToWorld.inverse() * rayWS;
	Ray2d rayLS = rayWS;

	BoundingBox2d bboxLS;
	bboxLS.reset();
	bboxLS.min = P2d(0.0f,0.0f);
	bboxLS.max = P2d(1.0f,1.0f);

	double tnearLS, tfarLS;
	double maxt;
	P2d minP, maxP, minP_ls, maxP_ls;
	ratio = 1.0;
	if( bboxLS.rayIntersect(rayLS, tnearLS, tfarLS) )
	{
		// here we add Epsilon in order to make sure that the intersection point is within its first voxel
		minP_ls = rayLS(tnearLS+Epsilon);
		maxP_ls = rayLS(tfarLS);
		//minP = localToWorld*minP_ls;
		minP = minP_ls;
		//maxP = localToWorld*maxP_ls;
		maxP = maxP_ls;
		m_tstart = (minP - rayWS.o).norm()*sign(tnearLS);
		maxt = (maxP - rayWS.o).norm()*sign(tfarLS);

		m_tstart = std::max( m_tstart, rayWS.mint );
		maxt = std::min( maxt, rayWS.maxt );
		if(maxt-m_tstart<=0)
		{
			if(debug)
				std::cout << "DDA2D::initialize maxt-m_tstart<=0\n";

			// no intersection
			return false;
		}

		ratio = (maxt-m_tstart)/(tfarLS-tnearLS+Epsilon);
		m_tend = maxt;
	}else
	{
		if(debug)
			std::cout << "DDA3D::initialize ray does not intersect grid bound\n";
		return false;
	}


	// convert localspace position of entry point to voxelspace
	vox.x() = std::min(res.x(), std::max(int(minP_ls.x()*res.x()), 0));
	vox.y() = std::min(res.y(), std::max(int(minP_ls.y()*res.y()), 0));

	for( int i=0;i<2;++i )
	{
		if( rayLS.d[i] > 0.0 )
		{
			m_step[i] = 1;
			end[i] = res[i];
		}else
		{
			m_step[i] = -1;
			end[i] = -1;
		}
	}

	// delta
	delta.x() = m_step.x()/rayLS.d.x();
	delta.y() = m_step.y()/rayLS.d.y();

	// tmax
	tmax.x()=diff_distance( minP_ls.x()*res.x(), rayLS.d.x() );
	tmax.y()=diff_distance( minP_ls.y()*res.y(), rayLS.d.y() );

	if(std::isinf(delta.x()))
		tmax.x() = std::numeric_limits<double>::infinity();
	if(std::isinf(delta.y()))
		tmax.y() = std::numeric_limits<double>::infinity();
	if( debug )
	{
		std::cout << "DDA2D::initialize delta=" << delta.toString() << std::endl;
		std::cout << "DDA2D::initialize tmax=" << tmax.toString() << std::endl;
		std::cout << "DDA2D::initialize rayLS.d=" << rayLS.d.toString() << std::endl;
	}

	m_lastmax = m_tstart;

	if( tmax.x() < tmax.y() )
	{
		axis = 0;
	}else
	{
		axis = 1;
	}

	return true;
}

bool DDA2D::initialize( const V2i& resolution, const P2d& p0, const P2d& p1, bool debug)
{
	res = resolution;

	//double maxt;
	//V2d d = normalized(p1-p0, maxt);
	V2d d = p1-p0;
	double maxt = d.norm();
	if(maxt > 0.0)
		d /= maxt;
	else
		return false;


	m_tstart = 0.0;
	m_tend = maxt;
	ratio = 1.0;

	// convert localspace position of entry point to voxelspace
	vox.x() = std::min(res.x(), std::max(int(p0.x()), 0));
	vox.y() = std::min(res.y(), std::max(int(p0.y()), 0));

	// convert localspace position of end point to voxelspace
	end.x() = std::min(res.x(), std::max(int(p1.x()), 0));
	end.y() = std::min(res.y(), std::max(int(p1.y()), 0));
	//std::cout << "end=" << end.toString() << std::endl;
	//return false;

	for( int i=0;i<2;++i )
	{
		if( d[i] > 0.0 )
		{
			m_step[i] = 1;
			//end[i] = res[i];
		}else
		{
			m_step[i] = -1;
			//end[i] = -1;
		}
	}

	// delta
	delta.x() = m_step.x()/d.x();
	delta.y() = m_step.y()/d.y();

	// tmax
	tmax.x()=diff_distance( p0.x(), d.x() );
	tmax.y()=diff_distance( p0.y(), d.y() );

	if(std::isinf(delta.x()))
		tmax.x() = std::numeric_limits<double>::infinity();
	if(std::isinf(delta.y()))
		tmax.y() = std::numeric_limits<double>::infinity();


	m_lastmax = m_tstart;

	if( tmax.x() < tmax.y() )
	{
		axis = 0;
	}else
	{
		axis = 1;
	}

	return true;
}

bool DDA2D::step( PixelStep& info, bool debug )
{
	// check if we have left the grid
	if(vox[axis] == end[axis])
		return false;

	// do next step
	if( tmax.x() < tmax.y() )
	{
		axis = 0;
	}else
	{
		axis = 1;
	}

	info.min = m_lastmax;
	info.max = m_tstart + tmax[axis]*ratio/res[axis];
	m_lastmax = info.max;
	if(debug)
	{
		std::cout << "DDA2D info.min=" << info.min << std::endl;
		std::cout << "DDA2D m_lastmax=" << m_lastmax << std::endl;
		std::cout << "DDA2D m_tstart=" << m_tstart << std::endl;
		std::cout << "DDA2D axis=" << axis << std::endl;
		std::cout << "DDA2D tmax[axis]=" << tmax[axis] << std::endl;
		std::cout << "DDA2D ratio=" << ratio << std::endl;
		std::cout << "DDA2D res[axis]=" << res[axis] << std::endl;
	}
	info.coord = vox;

	// update
	vox[axis] += m_step[axis];
	tmax[axis] += delta[axis];
	return true;
}





void draw_wu_line( short X0, short Y0, short X1, short Y1, short BaseColor, short NumLevels, unsigned short IntensityBits, DrawPixelFunction DrawPixel)
{
	unsigned short IntensityShift, ErrorAdj, ErrorAcc;
	unsigned short ErrorAccTemp, Weighting, WeightingComplementMask;
	short DeltaX, DeltaY, Temp, XDir;

	// Make sure the line runs top to bottom
	if (Y0 > Y1) {
		Temp = Y0; Y0 = Y1; Y1 = Temp;
		Temp = X0; X0 = X1; X1 = Temp;
	}
	// Draw the initial pixel, which is always exactly intersected by the line and so needs no weighting
	DrawPixel(X0, Y0, BaseColor);

	if ((DeltaX = X1 - X0) >= 0)
	{
		XDir = 1;
	} else {
		XDir = -1;
		DeltaX = -DeltaX; // make DeltaX positive
	}
	// Special-case horizontal, vertical, and diagonal lines, which
	// require no weighting because they go right through the center of
	// every pixel
	if ((DeltaY = Y1 - Y0) == 0)
	{
		// Horizontal line
		while (DeltaX-- != 0)
		{
			X0 += XDir;
			DrawPixel(X0, Y0, BaseColor);
		}
		return;
	}
	if (DeltaX == 0)
	{
		// Vertical line
		do {
			Y0++;
			DrawPixel(X0, Y0, BaseColor);
		} while (--DeltaY != 0);
		return;
	}
	if (DeltaX == DeltaY)
	{
		// Diagonal line
		do {
			X0 += XDir;
			Y0++;
			DrawPixel(X0, Y0, BaseColor);
		} while (--DeltaY != 0);
		return;
	}
	// Line is not horizontal, diagonal, or vertical
	ErrorAcc = 0;  // initialize the line error accumulator to 0
				   // # of bits by which to shift ErrorAcc to get intensity level
	IntensityShift = 16 - IntensityBits;
	// Mask used to flip all bits in an intensity weighting, producing the
	// result (1 - intensity weighting) */
	WeightingComplementMask = NumLevels - 1;
	// Is this an X-major or Y-major line?
	if (DeltaY > DeltaX) {
		// Y-major line; calculate 16-bit fixed-point fractional part of a
		// pixel that X advances each time Y advances 1 pixel, truncating the
		// result so that we won't overrun the endpoint along the X axis
		ErrorAdj = ((unsigned long) DeltaX << 16) / (unsigned long) DeltaY;
		// Draw all pixels other than the first and last
		while (--DeltaY)
		{
			ErrorAccTemp = ErrorAcc;   // remember currrent accumulated error
			ErrorAcc += ErrorAdj;      // calculate error for next pixel
			if (ErrorAcc <= ErrorAccTemp) {
				// The error accumulator turned over, so advance the X coord
				X0 += XDir;
			}
			Y0++; // Y-major, so always advance Y
				  // The IntensityBits most significant bits of ErrorAcc give us the
				  // intensity weighting for this pixel, and the complement of the
				  // weighting for the paired pixel
			Weighting = ErrorAcc >> IntensityShift;
			DrawPixel(X0, Y0, BaseColor + Weighting);
			DrawPixel(X0 + XDir, Y0, BaseColor + (Weighting ^ WeightingComplementMask));
		}
		// Draw the final pixel, which is 
		// always exactly intersected by the line
		// and so needs no weighting
		DrawPixel(X1, Y1, BaseColor);
		return;
	}
	// It's an X-major line; calculate 16-bit fixed-point fractional part of a
	// pixel that Y advances each time X advances 1 pixel, truncating the
	// result to avoid overrunning the endpoint along the X axis */
	ErrorAdj = ((unsigned long) DeltaY << 16) / (unsigned long) DeltaX;
	// Draw all pixels other than the first and last
	while (--DeltaX)
	{
		ErrorAccTemp = ErrorAcc;   // remember currrent accumulated error
		ErrorAcc += ErrorAdj;      // calculate error for next pixel
		if (ErrorAcc <= ErrorAccTemp)
		{
			// The error accumulator turned over, so advance the Y coord
			Y0++;
		}
		X0 += XDir; // X-major, so always advance X
					// The IntensityBits most significant bits of ErrorAcc give us the
					// intensity weighting for this pixel, and the complement of the
					// weighting for the paired pixel
		Weighting = ErrorAcc >> IntensityShift;
		DrawPixel(X0, Y0, BaseColor + Weighting);
		DrawPixel(X0, Y0 + 1, BaseColor + (Weighting ^ WeightingComplementMask));
	}
	// Draw the final pixel, which is always exactly intersected by the line
	// and so needs no weighting
	DrawPixel(X1, Y1, BaseColor);
}