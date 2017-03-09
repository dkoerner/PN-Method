#pragma once

#include <cuda/cumath/cumath.h>


struct LightSample
{
	cumath::V3f refP; // reference point towards which we want the lightsource to be sampled
	cumath::V3f p;
	cumath::V3f d; // unit length direction vector pointing from p towards refP (p->refP)
	float distance; // distance between refP and p
};

// NB: hardcoded to a directional light
struct CudaLight
{
	cumath::V3f m_direction;
	float m_distance;
	cumath::V3f m_radiance;

	DEVICE cumath::V3f sample( LightSample& ls, cumath::RNG& rng )const
	{
		ls.p = ls.refP - m_direction*m_distance;
		ls.d = m_direction;
#ifdef __CUDACC__
		ls.distance = CUDART_INF_F;
#else
		ls.distance = std::numeric_limits<float>::infinity();
#endif

		//NB: there is no inverse squared falloff due to the directional lightsource
		// being a delta distribution in solid angle domain
		//TODO: take surface normal at reference location into account once we do surfaces too
		return m_radiance;
	}
};

/*
struct PhaseFunction
{
	typedef std::shared_ptr<PhaseFunction> Ptr;

	// wi is the incident light direction (pointing from the scattering point towards the light)
	// wo is the outgoing light direction (pointing from the scattering point towards the camera)
	virtual Color3f eval( const V3d& wi, const V3d& wo )const=0;

	// wi point outwards, away from the scattering event
	// wo points outwards, away from the scattering event
	// this is inline with the convention for bsdfs and mitsuba
	// pdf is given in solid angle measure
	virtual Color3f sample( const V3d& wi, V3d& wo, double& pdf, RNGd& rng )const=0;
};
*/
