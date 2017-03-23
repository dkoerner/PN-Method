// http://stackoverflow.com/questions/12340384/are-there-advantages-to-using-the-cuda-vector-types
#pragma once

#include "common.h"

namespace cumath
{
	///
	/// \brief simple vector class
	///
	template<typename T>
	struct Vec2
	{

		HOST_DEVICE Vec2( );
		HOST_DEVICE Vec2( const T &x, const T &y );
		HOST_DEVICE Vec2( const T &xy );
		template <typename S> HOST_DEVICE Vec2(const Vec2<S> &v);
		HOST_DEVICE ~Vec2( );


		/*
        void              set( const T &x, const T &y, const T &z );
		*/

		HOST_DEVICE T                                    getLength( void )const; ///< returns the cartesian length of the vector
		/*
		T                             getSquaredLength( void )const; ///< returns the un-square-rooted length of the vector 
		void                          setLength( const T &fLength ); ///< scales the vector to the specified length
		*/
		HOST_DEVICE void           normalize( void ); ///< normalizes the vector
		HOST_DEVICE Vec2<T>        normalized( void ) const; ///< returnes normalized version if this the vector
		HOST_DEVICE Vec2<T>        normalized( float& length ) const; ///< returnes normalized version and the length
		/*
		void                                         negate( void ); ///< negates the vector

		void                       reflect( const Vec3<T> &normal ); ///< reflects the vector at the given normal



		bool                       operator==( const Vec3<T> &rhs );
		bool                       operator!=( const Vec3<T> &rhs );
		*/
		HOST_DEVICE bool           operator+=( const Vec2 &rhs );
		/*
		bool                       operator-=( const Vec3<T> &rhs );

		bool                             operator+=( const T &rhs );
		bool                             operator-=( const T &rhs );
		bool                             operator*=( const T &rhs );
		bool                             operator/=( const T &rhs );
		*/
		HOST_DEVICE const T&                          operator[]( int i ) const;
		HOST_DEVICE T&                                      operator[]( int i );

		union
		{
			struct
			{
				T x, y;
			};
			T v[2];
		};
	};



	template<typename T>
	HOST_DEVICE Vec2<T>::Vec2()
	{
		x=(T)0.0; y=(T)0.0;
    }

	template<typename T>
	HOST_DEVICE Vec2<T>::Vec2( const T &x, const T &y )
    {
		this->x=x; this->y=y;
    }

	template<typename T>
	HOST_DEVICE Vec2<T>::Vec2( const T &xy ) : x(xy), y(xy)
	{
	}

	template <typename T>
	template <typename S>
	inline
	HOST_DEVICE Vec2<T>::Vec2 (const Vec2<S> &v)
	{
		x = T (v.x);
		y = T (v.y);
	}

	template<typename T>
	HOST_DEVICE Vec2<T>::~Vec2( )
	{
	}
		/*
	template<typename T>
	void Vec2<T>::set( const T &x, const T &y, const T &z )
	{
		this->x=x; this->y=y; this->z=z;
	}
	*/
	template<typename T>
	HOST_DEVICE T Vec2<T>::getLength( void ) const
	{
		return sqrt( x*x + y*y);
	}

	/*

	template<typename T>
	T Vec2<T>::getSquaredLength( void ) const
	{
		return x*x + y*y + z*z;
	}

	template<typename T>
	void Vec2<T>::setLength( const T &length )
	{
		normalize();

		x*=length;
		y*=length;
		z*=length;
	}
	*/
	template<typename T>
	HOST_DEVICE void Vec2<T>::normalize( void )
	{
		T length = getLength();

		if( length != (T)0.0 )
		{
			x /= length;
			y /= length;
		}
	}


	///< returnes normalized version if this the vector
	template<typename T>
	HOST_DEVICE Vec2<T> Vec2<T>::normalized( void ) const
	{
	}

	template<typename T>
	HOST_DEVICE Vec2<T> Vec2<T>::normalized( float& length ) const
	{
		length = getLength();

		if( length != (T)0.0 )
			return Vec2<T>(x/length, y/length);
		else
			return Vec2<T>((T)0.0, (T)0.0);
	}

	/*
	template<typename T>
	void Vec2<T>::negate( void )
	{
		x=-x; y=-y; z=-z;
	}

	// reflects the vector at the given normal
	template<typename T>
	void Vec2<T>::reflect( const Vec2<T> &normal )
	{
		normalize();

		Vec2<T>	temp( x, y, z );
		float	value;
		
		value = dot( normal, *this );
		value *= (T)(2.0);

		x = normal.x * value;
		y = normal.y * value;
		z = normal.z * value;

		x -= temp.x;
		y -= temp.y;
		z -= temp.z;
	}

	template<typename T>
	bool Vec2<T>::operator==( const Vec2<T> &rhs )
	{
		// TODO: fix this
		if( (abs(x - rhs.x) < (T)0.00001) && (abs(y - rhs.y) < (T)0.00001) && (abs(z - rhs.z) < (T)0.00001) )
			return true;
		else
			return false;
	}

	template<typename T>
	bool Vec2<T>::operator!=( const Vec2<T> &rhs )
	{
		return !((*this)==rhs);
	}
	*/
	template<typename T>
	bool Vec2<T>::operator+=( const Vec2<T> &rhs )
	{
		x+=rhs.x;
		y+=rhs.y;

		return true;
	}
	/*
	template<typename T>
	bool Vec2<T>::operator-=( const Vec2<T> &rhs )
	{
		x-=rhs.x;
		y-=rhs.y;
		z-=rhs.z;

		return true;
	}

	template<typename T>
	bool Vec2<T>::operator+=( const T &rhs )
	{
		x+=rhs;
		y+=rhs;
		z+=rhs;

		return true;
	}

	template<typename T>
	bool Vec2<T>::operator-=( const T &rhs )
	{
		x-=rhs;
		y-=rhs;
		z-=rhs;

		return true; 
	}

	template<typename T>
	bool Vec2<T>::operator*=( const T &rhs )
	{
		x*=rhs;
		y*=rhs;
		z*=rhs;

		return true;
	}

	template<typename T>
	bool Vec2<T>::operator/=( const T &rhs )
	{
		x/=rhs;
		y/=rhs;
		z/=rhs;

		return true; 
	}
	*/

	template<typename T>
	HOST_DEVICE const T& Vec2<T>::operator[]( int i ) const
	{
		return v[i];
	}

	template<typename T>
	HOST_DEVICE T& Vec2<T>::operator[]( int i )
	{
		return v[i];
	}

	/*

	// stream output
	template <class T>
	std::ostream &
	operator << (std::ostream &s, const Vec2<T> &v)
	{
		return s << '(' << v.x << ' ' << v.y << ' ' << v.z <<')';
	}

		*/

	// shortcuts
	typedef Vec2<float> Vec2f;
	typedef Vec2<double> Vec2d;
	typedef Vec2<int> Vec2i;
	typedef Vec2<float> V2f;
	typedef Vec2<double> V2d;
	typedef Vec2<int> V2i;

}
