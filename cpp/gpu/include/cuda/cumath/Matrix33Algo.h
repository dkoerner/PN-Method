#pragma once


#include "Matrix33.h"


namespace cumath
{

// same as transform
template<typename T>
HOST_DEVICE inline Vec2<T> operator*( const Vec2<T> &lhs, const Matrix33<T> &rhs )
{
	return Vec2<T>(lhs.x*rhs._11 + lhs.y*rhs._21 + rhs._31,
				 lhs.x*rhs._12 + lhs.y*rhs._22 + rhs._32 );
}



} // namespace cumath
