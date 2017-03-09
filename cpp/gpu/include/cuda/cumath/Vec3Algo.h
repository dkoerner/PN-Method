#pragma once
#include "Vec3.h"


namespace cumath
{



	template<typename T>
	HOST_DEVICE inline Vec3<T> operator-( const Vec3<T> &rhs )
	{
		return Vec3<T>( -rhs.x, -rhs.y, -rhs.z );
	}

	template<typename T>
	HOST_DEVICE inline Vec3<T> operator+( const Vec3<T> &lhs, const Vec3<T> &rhs )
	{
		return Vec3<T>( lhs.x+rhs.x, lhs.y+rhs.y, lhs.z+rhs.z );
	}

	template<typename T>
	HOST_DEVICE inline Vec3<T> operator-( const Vec3<T> &lhs, const Vec3<T> &rhs )
	{
		return Vec3<T>( lhs.x-rhs.x, lhs.y-rhs.y, lhs.z-rhs.z );
	}

	template<typename T>
	HOST_DEVICE inline Vec3<T> operator/( const Vec3<T> &lhs, const Vec3<T> &rhs )
	{
		return Vec3<T>( lhs.x/rhs.x, lhs.y/rhs.y, lhs.z/rhs.z );
	}		

	template<typename T>
	HOST_DEVICE inline Vec3<T> operator+( const Vec3<T> &lhs, const T &rhs )
	{
		return Vec3<T>( lhs.x+rhs, lhs.y+rhs, lhs.z+rhs );
	}

	template<typename T>
	HOST_DEVICE inline Vec3<T> operator-( const Vec3<T> &lhs, const T &rhs )
	{
		return Vec3<T>( lhs.x-rhs, lhs.y-rhs, lhs.z-rhs );
	}

	template<typename T>
	HOST_DEVICE inline Vec3<T> operator*( const Vec3<T> &lhs, const T &rhs )
	{
		return Vec3<T>( lhs.x*rhs, lhs.y*rhs, lhs.z*rhs );
	}

	template<typename T>
	HOST_DEVICE inline Vec3<T> operator/( const Vec3<T> &lhs, const T &rhs )
	{
		return Vec3<T>( lhs.x/rhs, lhs.y/rhs, lhs.z/rhs );
	}

	template<typename T>
	HOST_DEVICE inline Vec3<T> operator+( const T &lhs, const Vec3<T> &rhs )
	{
		return (rhs+lhs);
	}

	template<typename T>
	HOST_DEVICE inline Vec3<T> operator-( const T &lhs, const Vec3<T> &rhs )
	{
		return Vec3<T>( lhs-rhs.x, lhs-rhs.y, lhs-rhs.z );
	}

	template<typename T>
	HOST_DEVICE inline Vec3<T> operator*( const T &lhs, const Vec3<T> &rhs )
	{
		return (rhs*lhs);
	}

	template<typename T>
	HOST_DEVICE inline Vec3<T> operator*( const Vec3<T> &lhs, const Vec3<T> &rhs )
	{
		return Vec3<T>(rhs.x*lhs.x, rhs.y*lhs.y, rhs.z*lhs.z);
	}

	template<typename T>
	HOST_DEVICE inline Vec3<T> operator/( const T &lhs, const Vec3<T> &rhs )
	{
		return (rhs/lhs);
	}

	template<typename T>
	HOST_DEVICE inline Vec3<T> normalize( const Vec3<T> &vector )
	{
		Vec3<T> result = vector;
		result.normalize();
		return result;
	}
		/*
	template<typename T>
	inline T dotProduct( const Vec3<T> &lhs, const Vec3<T> &rhs )
	{
		return (lhs.x*rhs.x + lhs.y*rhs.y + lhs.z*rhs.z);
	}
	*/

	template<typename T>
	HOST_DEVICE inline T dot( const Vec3<T> &a, const Vec3<T> &b )
	{ 
		return a.x * b.x + a.y * b.y + a.z * b.z;
	}


	/*
	template<typename T>
	inline void dotProduct( T &result, const Vec3<T> &lhs, const Vec3<T> &rhs )
	{
		result = (lhs.x*rhs.x + lhs.y*rhs.y + lhs.z*rhs.z);
	}


	template<typename T>
	inline Vec3<T> crossProduct( const Vec3<T> &lhs, const Vec3<T> &rhs )
	{
		return Vec3<T>( lhs.y*rhs.z - lhs.z*rhs.y,
					  lhs.z*rhs.x - lhs.x*rhs.z,
					  lhs.x*rhs.y - lhs.y*rhs.x );
	}

	template<typename T>
	inline void crossProduct( Vec3<T> &result, const Vec3<T> &lhs, const Vec3<T> &rhs )
	{
		result.x = lhs.y*rhs.z - lhs.z*rhs.y;
		result.y = lhs.z*rhs.x - lhs.x*rhs.z;
		result.z = lhs.x*rhs.y - lhs.y*rhs.x;
	}


	template<typename T>
	inline Vec3<T> cross( const Vec3<T> &lhs, const Vec3<T> &rhs )
	{
		return Vec3<T>( lhs.y*rhs.z - lhs.z*rhs.y,
					  lhs.z*rhs.x - lhs.x*rhs.z,
					  lhs.x*rhs.y - lhs.y*rhs.x );
	}

	template<typename T>
	inline void cross( Vec3<T> &result, const Vec3<T> &lhs, const Vec3<T> &rhs )
	{
		result.x = lhs.y*rhs.z - lhs.z*rhs.y;
		result.y = lhs.z*rhs.x - lhs.x*rhs.z;
		result.z = lhs.x*rhs.y - lhs.y*rhs.x;
	}

	// For a given incident vector I and surface normal N reflect returns the reflection direction calculated as I - 2.0 * dot(N, I) * N. N should be normalized in order to achieve the desired result.
	template<typename T>
	inline Vec3<T> reflect( const math::Vec3<T> &i, const math::Vec3<T> &n )
	{
		return -i + 2.0f*dotProduct(n,i)*n;
	}


	// returns index of the component with the smallest value
	template<typename T>
	int nondominantAxis( const math::Vec3<T> &v )
	{
		if( abs(v.x) <= abs(v.y) )
		{
			if( abs(v.z) <= abs(v.x) )
				return 2;
			else
				return 0;
		}else
		{
			if( abs(v.z) <= abs(v.y) )
				return 2;
			else
				return 1;
		}
	}

	template<typename T>
	Vec3<T> baseVec3( int i )
	{
		switch(i)
		{
		case 0:return Vec3<T>((T)(1.0), (T)(0.0), (T)(0.0));
		case 1:return Vec3<T>((T)(0.0), (T)(1.0), (T)(0.0));
		default:
		case 2:return Vec3<T>((T)(0.0), (T)(0.0), (T)(1.0));
		};
	}

	template<typename T>
	inline T length( Vec3<T> &v )
	{
		return v.getLength();
	}


	*/

	template<typename T>
	HOST_DEVICE inline Vec3<T> min( Vec3<T> &v0, Vec3<T> &v1 )
	{
		return Vec3<T>( v0.x < v1.x ? v0.x : v1.x, v0.y < v1.y ? v0.y : v1.y, v0.z < v1.z ? v0.z : v1.z );
	}

	template<typename T>
	HOST_DEVICE inline Vec3<T> max( Vec3<T> &v0, Vec3<T> &v1 )
	{
		return Vec3<T>( v0.x > v1.x ? v0.x : v1.x, v0.y > v1.y ? v0.y : v1.y, v0.z > v1.z ? v0.z : v1.z );
	}

}