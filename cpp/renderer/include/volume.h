#pragma once

#include <memory>
#include <fields/Field.h>
#include <math/frame.h>
#include <math/transform.h>
#include <math/bbox.h>
#include <math/color.h>
#include <math/rng.h>






struct PhaseFunction
{
	typedef std::shared_ptr<PhaseFunction> Ptr;

	// wi is the incident light direction (pointing from the scattering point towards the light)
	// wo is the outgoing light direction (pointing from the scattering point towards the camera)
	virtual double eval( const V3d& wi, const V3d& wo )const=0;

	// wi point outwards, away from the scattering event
	// wo points outwards, away from the scattering event
	// this is inline with the convention for bsdfs and mitsuba
	// pdf is given in solid angle measure
	virtual double sample( const V3d& wi, V3d& wo, double& pdf, RNGd& rng )const=0;
};

struct IsotropicPhase : public PhaseFunction
{
	typedef std::shared_ptr<IsotropicPhase> Ptr;

	virtual double sample( const V3d& wi, V3d& wo, double& pdf, RNGd& rng ) const override
	{
		wo = sampleSphere<double>(rng);
		pdf = INV_FOURPI;
		return 1.0; // eval/pdf = 1.0
	}

	virtual double eval( const V3d& wi, const V3d& wo )const override
	{
		return INV_FOURPI;
	}
};

struct HGPhase : public PhaseFunction
{
	typedef std::shared_ptr<HGPhase> Ptr;

	virtual double sample( const V3d& wi, V3d& wo, double& pdf, RNGd& rng ) const override
	{
		/*
		wo = sampleSphere<double>(rng);
		pdf = INV_FOURPI;
		return 1.0; // eval/pdf = 1.0
		*/
		return 0.0;
	}

	virtual double eval( const V3d& wi, const V3d& wo )const override
	{
		double temp = 1.0 + m_g*m_g + 2.0 * m_g * dot(wi, wo);
		return INV_FOURPI * (1 - m_g*m_g) / (temp * std::sqrt(temp));
	}


private:
	double m_g;
};

/*
struct PhaseFunctionField
{
	typedef std::shared_ptr<PhaseFunctionField> Ptr;

	// pWS is the world position at which the phase function is to be evaluated
	// wi is the incident light direction (pointing from the scattering point towards the light)
	// wo is the outgoing light direction (pointing from the scattering point towards the camera)
	virtual double eval( const P3d& pWS, const V3d& wi, const V3d& wo )const=0;


	// pWS is the world position at which the phase function is to be sampled
	// wi point outwards, away from the scattering event
	// wo points outwards, away from the scattering event
	// this is inline with the convention for bsdfs and mitsuba
	// pdf is given in solid angle measure
	virtual double sample( const P3d& pWS, const V3d& wi, V3d& wo, double& pdf, RNGd& rng )const=0;

	// transforms the phase function field from local space to some worldspace
	virtual void setLocalToWorld( const Transformd& localToWorld )=0;
};
*/

struct Volume
{
	typedef std::shared_ptr<Volume> Ptr;

	Volume();

	void setExtinctionAlbedo( Field3d::Ptr extinction, Field3d::Ptr albedo );
	//void setAbsorptionScattering( Field3d::Ptr absorption, Field3d::Ptr scattering );

	V3d evalExtinction( const P3d& pWS, bool debug = false )const;
	V3d evalAlbedo( const P3d& pWS, bool debug = false )const;

	double evalPhase( const P3d& pWS, const V3d& wi, const V3d& wo )const;
	double samplePhase( const P3d& pWS, const V3d& wi, V3d& wo, double& pdf, RNGd& rng )const;

	//virtual Transformd getLocalToWorld()const{return Transformd();}

	P3d localToWorld(const P3d& pLS)const;
	P3d worldToLocal(const P3d& pWS)const;

	//Box3d getBound()const=0; // aa worldspace bounding box
	void setLocalToWorld( const Transformd& localToWorld );
	void setBound( Box3d boundWS );
	const Box3d& getBound()const;

	// the worldspace bounding box returned by getBound may be too conservative if the volume is rotated
	// intersectBound can transform the ray into the volumes localspace and intersect the bounding box there...
	bool intersectBound(const Ray3d& rayWS, double& mint, double& maxt, bool debug = false)const;

	// this is required for delta tracking
	V3d getMaxExtinction()const;

	std::string toString()const;

	Eigen::MatrixXd getSlice( double offset )const;

private:
	// rgb dependent albedo and extinction values (evaluated in local space)
	Field3d::Ptr       m_field_extinction;
	Field3d::Ptr       m_field_albedo;
	PhaseFunction::Ptr m_phaseFunction;
	//PhaseFunctionField::Ptr phase_field;

	//V3d          m_extinction_min; // used for residual ratio tracking
	V3d           m_extinction_max; // used for delta sampling
	Transformd    m_localToWorld; // transform
	Transformd    m_worldToLocal; // transform
	Box3d         m_bboxLS; // bounding box in local space (0,0,0->1,1,1)
	Box3d         m_bboxWS; // bounding box in world space
};

/*



//VoxelGridField::Ptr loadBGEO(const std::string &path);


namespace volumes
{
	Volume::Ptr C60();
	Volume::Ptr nebulae();
	Volume::Ptr homogeneous( const V3d& extinction, const V3d& albedo = V3d(0.8), Transformd localToWorld = Transformd() );
	Volume::Ptr scarf();

	PhaseFunction::Ptr phase_isotropic();
	PhaseFunction::Ptr phase_sggx_specular( Framed frame = Framed(V3d(0.0, 0.0, 1.0)), V3d projectedAreas=V3d(1.0, 1.0, 1.0) );


	PhaseFunctionField::Ptr phasefield_sggx();
}
*/
