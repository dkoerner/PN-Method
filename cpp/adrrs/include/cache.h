#pragma once

#include <memory>
#include <math/vector.h>
#include <math/color.h>

struct Cache
{
	typedef std::shared_ptr<Cache> Ptr;

	virtual Color3f eval( const P3d& pWS )const=0; // returns fluence (in units of power) at position pWS
	virtual Color3f eval( const P3d& pWS, const V3d& d )const=0;
};



