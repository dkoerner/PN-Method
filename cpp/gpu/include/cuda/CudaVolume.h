#pragma once

#include <cuda/CudaVoxelGrid.h>



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


namespace cumath
{
	struct Transform
	{
		cumath::M44f matrix;
		cumath::M44f inverse;
	};
}

struct CudaVolume
{
	DEVICE cumath::V3f evalExtinction( const cumath::V3f& pWS, bool debug = false )const
	{
		cumath::V3f pVS = m_sigma_t->localToVoxel(pWS*m_localToWorld.inverse);
		float sigma_t = m_sigma_t->eval(pVS);
		return cumath::V3f(sigma_t);
	}

	DEVICE cumath::V3f evalAlbedo( const cumath::V3f& pWS, bool debug = false )const
	{
		return cumath::V3f(0.8f);
	}

	//const PhaseFunction* getPhaseFunction()const=0;

	DEVICE cumath::V3f getMaxExtinction()const
	{
		return m_sigma_t_max;
	}

	DEVICE const cumath::Box3f& getBound()const
	{
		return m_boundWS;
	}

	DEVICE cumath::M44f getLocalToWorld()const
	{
		return m_localToWorld.matrix;
	}

	// the worldspace bounding box returned by getBound may be too conservative if the volume is rotated
	// intersectBound can transform the ray into the volumes localspace and intersect the bounding box there...
	DEVICE bool intersectBound(const cumath::Ray3f& rayWS, float& mint, float& maxt, bool debug = false)const
	{
		return cumath::intersectionRayBox( rayWS, m_boundWS, mint, maxt, debug );
	}

	// wi is the incident light direction (pointing from the scattering point towards the light)
	// wo is the outgoing light direction (pointing from the scattering point towards the camera)
	DEVICE cumath::V3f evalPhaseFunction( const cumath::V3f& wi, const cumath::V3f& wo )const
	{
		return cumath::V3f(MATH_INV_FOURPIf, MATH_INV_FOURPIf, MATH_INV_FOURPIf);
	}

	DEVICE cumath::V3f samplePhaseFunction( const cumath::V3f& wi, cumath::V3f& wo, float& pdf, cumath::RNG& rng ) const
	{
		wo = cumath::sampleSphere<double>(rng);
		pdf = MATH_INV_FOURPIf;
		return cumath::V3f(1.0); // eval/pdf = 1.0
	}

	cumath::Transform m_localToWorld;
	cumath::Box3f m_boundWS;
	CudaVoxelGrid<float>* m_sigma_t;
	cumath::V3f m_sigma_t_max;
};
