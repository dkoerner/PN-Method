#pragma once
#include "Vec2.h"

//#include "BoundingBox3.h"
#include <limits>




namespace cumath
{

	/// \brief a class which holds an origin and a target and is used by intersection routines
	template<typename T>
	class Ray2
	{
	public:

		HOST_DEVICE Ray2();                                                                                   // constructor
		HOST_DEVICE Ray2( const cumath::Vec2<T> &origin, const cumath::Vec2<T> &direction, const T &_tmin, const T &_tmax );       // constructor


		HOST_DEVICE Vec2<T>                                              getPosition( T t )const; // returns origin+direction*t

		Vec2<T>                                                                                o; // point in space where the ray originates from
		Vec2<T>                                                                                d; // normalized direction of the ray
		T                                                                             tmin, tmax; // valid ray segment
	};



	// constructor
	template<typename T>
	HOST_DEVICE Ray2<T>::Ray2() : tmin(0.0f), tmax(FLT_MAX)
	{
	}

	// constructor
	template<typename T>
	HOST_DEVICE Ray2<T>::Ray2( const cumath::Vec2<T> &origin, const cumath::Vec2<T> &direction, const T &_tmin, const T &_tmax ) : o(origin), d(direction), tmin(_tmin), tmax(_tmax)
	{
	}

	// returns origin+direction*t
	template<typename T>
	HOST_DEVICE Vec2<T> Ray2<T>::getPosition( T t )const
	{
		return o + d*t;
	}



	typedef Ray2<float> Ray2f;
	typedef Ray2<double> Ray2d;

} // namespace cumath


