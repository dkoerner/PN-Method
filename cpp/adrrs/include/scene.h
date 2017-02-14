#pragma once
#include <string>
#include <math/common.h>
#include <math/vector.h>
#include <math/bbox.h>



struct Camera;
struct Light;
struct Volume;
struct Integrator;
struct Scene
{
	// scene info
	std::string id;
	Box3d bound;   // bounding box of the scene (used for view finding)

	// camera
	const Camera* camera;

	// light
	const Light* light;

	// volume
	const Volume* volume;

	// integrator
	const Integrator* integrator;
};






