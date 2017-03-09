// http://stackoverflow.com/questions/12340384/are-there-advantages-to-using-the-cuda-vector-types
#pragma once

#include "common.h"

namespace cumath
{
	///
	/// \brief simple vector class
	///
	template<typename T>
	struct Vec4
	{

		HOST_DEVICE Vec4( );
		HOST_DEVICE Vec4( const T &x, const T &y, const T &z, const T &w );
		HOST_DEVICE Vec4( const T &xyzw );
		template <typename S> HOST_DEVICE Vec4(const Vec4<S> &v);
		HOST_DEVICE ~Vec4( );



		HOST_DEVICE bool           operator+=( const Vec4<T> &rhs );
		HOST_DEVICE const T&                          operator[]( int i ) const;
		HOST_DEVICE T&                                      operator[]( int i );

		union
		{
			struct
			{
				T x, y, z, w;
			};
			T v[4];
		};
	};



	template<typename T>
	HOST_DEVICE Vec4<T>::Vec4()
	{
		x=(T)0.0; y=(T)0.0; z=(T)0.0; w=(T)0.0;
    }

	template<typename T>
	HOST_DEVICE Vec4<T>::Vec4( const T &x, const T &y, const T &z, const T &w )
    {
		this->x=x; this->y=y; this->z=z; this->w=w;
    }

	template<typename T>
	HOST_DEVICE Vec4<T>::Vec4( const T &xyzw ) : x(xyzw), y(xyzw), z(xyzw), w(xyzw)
	{
	}

	template <typename T>
	template <typename S>
	inline
	HOST_DEVICE Vec4<T>::Vec4 (const Vec4<S> &v)
	{
		x = T (v.x);
		y = T (v.y);
		z = T (v.z);
		w = T (v.w);
	}

	template<typename T>
	HOST_DEVICE Vec4<T>::~Vec4( )
	{
	}



	template<typename T>
	bool Vec4<T>::operator+=( const Vec4<T> &rhs )
	{
		x+=rhs.x;
		y+=rhs.y;
		z+=rhs.z;
		w+=rhs.w;

		return true;
	}

	template<typename T>
	HOST_DEVICE const T& Vec4<T>::operator[]( int i ) const
	{
		return v[i];
	}

	template<typename T>
	HOST_DEVICE T& Vec4<T>::operator[]( int i )
	{
		return v[i];
	}


	// shortcuts
	typedef Vec4<float> Vec4f;
	typedef Vec4<double> Vec4d;
	typedef Vec4<int> Vec4i;
	typedef Vec4<float> V4f;
	typedef Vec4<double> V4d;
	typedef Vec4<int> V4i;

}
