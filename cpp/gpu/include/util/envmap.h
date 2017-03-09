#pragma once
#include <util/bitmap.h>

#include <vector>
#include <string>
#include <iostream>
#include <memory>
#include <math/vector.h>
#include <math/rng.h>
#include <math/color.h>
#include <math/transform.h>

struct EnvMap
{
	typedef std::shared_ptr<EnvMap> Ptr;

	EnvMap( const std::string& filename )
	{
		m_bitmap = Bitmap(filename);
		m_transform = Transformd();
	}

	EnvMap( const V2i& res = V2i(512, 256) )
	{
		m_bitmap = Bitmap(res);
		m_transform = Transformd();
	}

	// evaluate environment map
	Color3f eval( double theta, double phi )const
	{
		V3d d = sphericalDirection<double>(theta, phi);
		P2d uv = directionToUV(d);
		return m_bitmap.eval(uv);
	}
	P2d directionToUV( const V3d& d )const
	{
		// using formulas given in http://gl.ict.usc.edu/Data/HighResProbes/
		// with the difference that u=[0,1] (instead of [0,2]) and we negate z
		P2d uv( (1+std::atan2(d.x(), d.z())/M_PI)/2,
				 safe_acos(d.y())/M_PI );
		return uv;
	}
	V3d uvToDirection( const P2d& uv )const
	{
		// using formulas given in http://gl.ict.usc.edu/Data/HighResProbes/
		// with the difference that u=[0,1] (instead of [0,2]) and we negate z
		// azimuthal angle
		double phi = M_PI*(uv.x()*2.0-1.0);
		// elevation angle
		double theta = M_PI*uv.y();
		// TODO: I am confused...
		return V3d( std::sin(theta)*std::sin(phi), std::cos(theta), std::sin(theta)*cos(phi) );
		//return sphericalDirection( theta, phi );
	}
	P2d uvToXY(const P2d& uv)const
	{
		P2d xy(
			(uv.x()*(m_bitmap.cols()-1)),
			(uv.y()*(m_bitmap.rows()-1))
			);
		return xy;
	}

	P2d xyToUV(const P2d& xy)const
	{
		return P2d(
			(xy.x())/double(m_bitmap.cols()-1),
			(xy.y())/double(m_bitmap.rows()-1)
			);
	}

	V3d xyToDirection( const P2d& xy )const
	{
		return uvToDirection( xyToUV(xy) );
	}
	P2d directionToXY( const V3d& d )const
	{
		return uvToXY(directionToUV(d));
	}

	Bitmap& bitmap()
	{
		return m_bitmap;
	}

	void clear()
	{
		m_bitmap.fill(Color3f(0.0f, 0.0f, 0.0f));
	}

	void saveGeo( const std::string& filename, double exposure = 0.0);
//private:
	Transformd m_transform;
	Bitmap m_bitmap;
};


void rasterizeSphericalFunctionSphere(const std::string& filename, std::function<Color3f (double, double)> func, double exposure =0.0);
