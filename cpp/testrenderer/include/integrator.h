#pragma once
#include <math/color.h>


#include <volume.h>
#include <light.h>
#include <camera.h>
#include <util/bitmap.h>

#include <math/common.h>
#include <math/rng.h>

#include <cache.h>



struct Scene;
struct RenderTaskInfo;


struct RadianceQuery
{
	RadianceQuery():
		debug(false),
		pixel(-1, -1)
	{
	}

	Ray3d ray;
	V2i pixel;
	bool debug;

	Color3f transmittance;
};


struct Integrator
{
	typedef std::shared_ptr<Integrator> Ptr;

	virtual Color3f Li( const Scene* scene, RadianceQuery& rq, RNGd& rng )const=0;
	virtual Color3f sample_transmittance( const Scene* scene, const Ray3d& ray, double maxt, RNGd& rng )const=0;
	virtual std::string getId()const=0;
};


struct RenderTaskInfoGlobal
{
	RenderTaskInfoGlobal():
		debug_pixel(-1, -1)
	{

	}

	// always used ---
	const Scene* scene;

	// used for image rendering ---
	Bitmap* image;
	Bitmap* image_transmittance;
	P2i debug_pixel;
	Box2i crop_window;

	// used for fluence rendering ---
	field::VoxelGridField<double>::Ptr fluence_grid;
	Transformd localToWorld;
};


struct RenderTaskInfo
{
	RenderTaskInfo(int taskid, int numTasks):
		taskid(taskid),
		numTasks(numTasks),
		samples(0),
		rng(123+taskid)
	{
	}

	static RenderTaskInfo create( int taskid, int numTasks )
	{
		return RenderTaskInfo(taskid, numTasks);
	}

	// per task info ---
	int taskid;
	int numTasks;
	int samples;
	RNGd rng;

	static RenderTaskInfoGlobal g; // global data across all tasks
};


namespace integrators
{
	Integrator::Ptr dummy();
	Integrator::Ptr raymarcher(double stepsize = 0.1);
	Integrator::Ptr volume_path_tracer_cached(Cache::Ptr cache );
	Integrator::Ptr volume_path_tracer(int maxDepth = std::numeric_limits<int>::max(), bool doSingleScattering = true);
	Integrator::Ptr adrrs_volume_path_tracer( Bitmap::Ptr pixel_estimates, Field3d::Ptr fluence );
}

