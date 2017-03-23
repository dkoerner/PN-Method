#pragma once
#include "Vec2.h"




namespace cumath
{
	template<typename T>
	HOST_DEVICE inline Vec2<T> operator+( const Vec2<T> &lhs, const Vec2<T> &rhs )
	{
		return Vec2<T>( lhs.x+rhs.x, lhs.y+rhs.y );
	}

	template<typename T>
	HOST_DEVICE inline Vec2<T> operator-( const Vec2<T> &lhs, const Vec2<T> &rhs )
	{
		return Vec2<T>( lhs.x-rhs.x, lhs.y-rhs.y );
	}


	template<typename T>
	HOST_DEVICE inline Vec2<T> operator/( const Vec2<T> &lhs, const Vec2<T> &rhs )
	{
		return Vec2<T>( lhs.x/rhs.x, lhs.y/rhs.y );
	}

	template<typename T>
	HOST_DEVICE inline Vec2<T> operator+( const Vec2<T> &lhs, const T &rhs )
	{
		return Vec2<T>( lhs.x+rhs, lhs.y+rhs );
	}

	template<typename T>
	HOST_DEVICE inline Vec2<T> operator-( const Vec2<T> &lhs, const T &rhs )
	{
		return Vec2<T>( lhs.x-rhs, lhs.y-rhs);
	}

	template<typename T>
	HOST_DEVICE inline Vec2<T> operator*( const Vec2<T> &lhs, const T &rhs )
	{
		return Vec2<T>( lhs.x*rhs, lhs.y*rhs);
	}


	template<typename T>
	HOST_DEVICE inline Vec2<T> normalize( const Vec2<T> &vector )
	{
		Vec2<T> result = vector;
		result.normalize();
		return result;
	}
} // namespace cumath
