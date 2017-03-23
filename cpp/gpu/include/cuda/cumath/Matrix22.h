#pragma once
#include "Vec2.h"

namespace cumath
{

	template<typename T>
	struct Matrix22
	{
		union
		{
			struct
			{
				T _11, _12;
				T _21, _22;
			};
			T m[2][2];
			T ma[4];
		};
	};



	typedef Matrix22<float> Matrix22f;
	typedef Matrix22<double> Matrix22d;
	typedef Matrix22<float> M22f;
	typedef Matrix22<double> M22d;


} // namespace cumath
