#pragma once
#include "common.h"
#include <thrust/random.h>


namespace cumath
{
	struct RNG
	{
		HOST_DEVICE RNG()
		{
		}

		HOST_DEVICE RNG( int seed )
		{
			rng.seed(seed);
		}


		HOST_DEVICE float randomFloat()const
		{
			return float(rng()-thrust::default_random_engine::min)/float(thrust::default_random_engine::max - thrust::default_random_engine::min);
		}

		mutable thrust::default_random_engine rng;
	};
}