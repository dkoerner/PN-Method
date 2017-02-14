#pragma once
#include <memory>
#include <math/vector.h>
#include <math/color.h>






struct LightSample
{
	P3d refP; // reference point towards which we want the lightsource to be sampled
	P3d p;
	V3d d; // unit length direction vector pointing from p towards refP (p->refP)
	double distance; // distance between refP and p
};

struct Light
{
	typedef std::shared_ptr<Light> Ptr;


	// samples a position on the light source and returns radiance/pdf
	// the geometry term between the light source position and the reference location will be accounted for
	virtual Color3f sample( LightSample& ls )const =0;
};


namespace lights
{
	Light::Ptr point( const P3d& position = P3d(0.0, 0.0, 0.0) );
	Light::Ptr directional( const P3d& direction, const Box3d& scene_bound );
}


