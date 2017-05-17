#pragma once


#include "Vec2.h"
#include "Vec2Algo.h"
#include "Vec3.h"
#include "Vec3Algo.h"
#include "Vec4.h"
#include "Matrix22.h"
#include "Matrix33.h"
#include "Matrix33Algo.h"
#include "Matrix44.h"
#include "Matrix44Algo.h"
#include "Ray2.h"
#include "Ray3.h"
#include "BoundingBox3.h"
#include "RNG.h"

#include <math_constants.h>

#define MATH_PIf 3.14159265f
#define MATH_PI 3.14159265
#define MATH_NEG_INFTY_F 0xff800000
#define MATH_INV_FOURPIf   0.07957747154594766788f


namespace cumath
{

	HOST_DEVICE inline unsigned int randhash(unsigned int a)
	{
		a = (a+0x7ed55d16) + (a<<12);
		a = (a^0xc761c23c) ^ (a>>19);
		a = (a+0x165667b1) + (a<<5);
		a = (a+0xd3a2646c) ^ (a<<9);
		a = (a+0xfd7046c5) + (a<<3);
		a = (a^0xb55a4f09) ^ (a>>16);
		return a;
	}


	template<typename T, typename R>
	HOST_DEVICE inline T lerp( T x0, T x1, R t )
	{
		//return x0*((R)(1.0)-t) + x1*t;
		return x0 + t*(x1-x0);
	}

	template<typename T>
	HOST_DEVICE inline T max( T x, T y )
	{
		return x > y ? x : y;
	}

	template<typename T>
	HOST_DEVICE inline T min( T x, T y )
	{
		return x < y ? x : y;
	}

	template<typename T>
	HOST_DEVICE inline void swap( T &x, T &y )
	{
		T t = x;
		x = y;
		y = t;
	}

	template<typename T>
	HOST_DEVICE inline double safe_acos( const T& in )
	{
		return acos( max( T(-1.0), min(T(1.0), in)) );
	}


	template<typename T>
	HOST_DEVICE Vec3<T> sampleSphere( const RNG &rng, T *pdf = 0 )
	{
		if( pdf )
			*pdf = (T)((T)1.0/((T)4.0*MATH_PI));

		T u1 = rng.randomFloat();
		T u2 = rng.randomFloat();
		T z = (T)(1.0 - 2.0 * u1);
		T s = (T)(sqrt(max(0.0, 1.0 - z*z)));
		T phi =(T)(2.0 * MATH_PI * u2);
		return Vec3<T>(s*cos(phi), s*sin(phi), z);
	}

	template<typename T>
	HOST_DEVICE Vec3<T> sphericalDirection(T theta, T phi)
	{
		T sinTheta=sin(theta), cosTheta=cos(theta), sinPhi=sin(phi), cosPhi=cos(phi);

		//sincos(theta, &sinTheta, &cosTheta);
		//sincos(phi, &sinPhi, &cosPhi);

		return Vec3<T>(
		   sinTheta * cosPhi,
		   sinTheta * sinPhi,
		   cosTheta
		);
	}


}
