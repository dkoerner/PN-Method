#pragma once

#include <cuda/CudaVoxelGrid.h>
#include <cuda/CudaPixelBuffer.h>



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
		M44f matrix;
		M44f inverse;
	};
	struct Transform2D
	{
		M33f matrix;
		M33f inverse;
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


struct CudaSGGX2D
{
	HOST CudaSGGX2D( cumath::V2f frame_s = cumath::V2f(1.0f, 0.0f), cumath::V2f projectedDistances = cumath::V2f(1.0f, 1.0f) )
	{
		float S_11 = projectedDistances.x*projectedDistances.x;
		float S_22 = projectedDistances.y*projectedDistances.y;

		cumath::V2f frame_t(-frame_s.y, frame_s.x);

		float S_xx = S_11*frame_s.x*frame_s.x + S_22*frame_t.x*frame_t.x;
		float S_xy = S_11*frame_s.x*frame_s.y + S_22*frame_t.x*frame_t.y;
		float S_yx = S_xy;
		float S_yy = S_11*frame_s.y*frame_s.y + S_22*frame_t.y*frame_t.y;

		S._11 = S_xx;
		S._12 = S_xy;
		S._21 = S_yx;
		S._22 = S_yy;
	}



	HOST_DEVICE float projectedDistance( const cumath::V2f& d )const
	{
		float sigma_squared = S._11*d.x*d.x + S._22*d.y*d.y + 2.0*S._12*d.x*d.y;
		if(sigma_squared>0.0)
			return sqrt(sigma_squared);
		// conditional to avoid numerical errors
		return 0.0;
	}

	cumath::M22f S;
};

struct CudaVolume2D
{

	DEVICE cumath::V3f evalExtinction( const cumath::V2f& pWS, const cumath::V2f& d, bool debug = false )const
	{
		return evalDensity(pWS)*evalSGGX(pWS).projectedDistance(d);
	}

	DEVICE cumath::V3f evalDensity( const cumath::V2f& pWS)const
	{
		//cumath::V2f pVS = localToVoxel(pWS*m_localToWorld.inverse);
		//return m_density->eval(pVS);
		return m_density;
	}

	DEVICE const CudaSGGX2D& evalSGGX( const cumath::V2f& pWS )const
	{
		// homogeneous microflake distribution for now...
		return m_sggx;
	}

	/*
	DEVICE cumath::V3f evalAlbedo( const cumath::V3f& pWS, bool debug = false )const
	{
		return cumath::V3f(0.8f);
	}
*/

	// returns maximum extinction for all directions
	DEVICE cumath::V3f getMaxExtinction()const
	{
		return m_sigma_t_max;
	}

	/*
	DEVICE const cumath::Box3f& getBound()const
	{
		return m_boundWS;
	}
	*/

	DEVICE cumath::M33f getLocalToWorld()const
	{
		return m_localToWorld.matrix;
	}

	/*
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
	*/

	/*
	DEVICE cumath::V2f localToVoxel( const cumath::V2f& pLS )const
	{
		return cumath::V2f( pLS.x*m_density->m_width,
							pLS.y*m_density->m_height);
	}

	DEVICE cumath::V2f voxelToLocal( const cumath::V2f& pVS )const
	{
		return cumath::V2f( pVS.x/m_density->m_width,
							pVS.y/m_density->m_height);
	}
	*/

	cumath::Transform2D m_localToWorld;
	//cumath::Box3f m_boundWS;
	//CudaPixelBuffer* m_density;
	cumath::V3f m_density;
	cumath::V3f m_sigma_t_max;
	CudaSGGX2D m_sggx;
};
