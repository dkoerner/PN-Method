#pragma once
#include "Vec2.h"

namespace cumath
{

	template<typename T>
	struct Matrix33
	{
		union
		{
			struct
			{
				T _11, _12, _13;
				T _21, _22, _23;
				T _31, _32, _33;
			};
			T m[3][3];
			T ma[9];
		};
	};



	typedef Matrix33<float> Matrix33f;
	typedef Matrix33<double> Matrix33d;
	typedef Matrix33<float> M33f;
	typedef Matrix33<double> M33d;


} // namespace cumath
