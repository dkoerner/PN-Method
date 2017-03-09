#pragma once

#include <util/field.h>
#include <util/data.h>
#include <math/transform.h>
#include <math/bbox.h>
#include <math/rng.h>






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


struct Volume
{
	typedef std::shared_ptr<Volume> Ptr;


	virtual Color3f evalExtinction( const P3d& pWS, bool debug = false )const=0;
	virtual Color3f evalAlbedo( const P3d& pWS, bool debug = false )const=0;
	virtual const PhaseFunction* getPhaseFunction()const=0;
	virtual V3d getMaxExtinction()const=0;

	// returns bounding box in worldspace
	virtual Box3d getBound()const=0;

	virtual Transformd getLocalToWorld()const{return Transformd();}

	// the worldspace bounding box returned by getBound may be too conservative if the volume is rotated
	// intersectBound can transform the ray into the volumes localspace and intersect the bounding box there...
	virtual bool intersectBound(const Ray3d& rayWS, double& mint, double& maxt, bool debug = false)const=0;
};





//VoxelGridField::Ptr loadBGEO(const std::string &path);


namespace volumes
{
	Volume::Ptr C60();
	Volume::Ptr nebulae();


	/*
	Scene C60Detail();
	Scene rbf();
	Scene engine();
	Scene chameleon();
	*/
}

