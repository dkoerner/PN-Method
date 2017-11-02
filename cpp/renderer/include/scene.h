#pragma once
#include <string>
#include <math/common.h>
#include <math/vector.h>
#include <math/bbox.h>
#include <math/rng.h>

//#include <light.h>

struct Intersection
{
	void reset()
	{
		t = std::numeric_limits<double>::infinity();
	}

	P3d p;
	double t;
};

struct Light;
struct LightSample;
struct Camera;
struct Volume;
struct Integrator;
struct Scene
{
	// scene info
	//std::string id;
	//Box3d bound;   // bounding box of the scene (used for view finding)

	// camera
	const Camera* camera;

	// light
	const Light* light;

	// volume
	const Volume* volume;

	// integrator
	const Integrator* integrator;


	bool intersect( const Ray3d& ray, Intersection& its )const;
	V3d sample_attenuated_directlight( LightSample& ls, RNGd& rng, bool debug = false )const;
};


