


#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>

#include <volume.h>
#include <light.h>
#include <camera.h>
#include <integrator.h>
#include <Image.h>

#include <lights/pointlight.h>
#include <lights/directionallight.h>
#include <integrators/simplept.h>
#include <integrators/pnispt.h>
#include <integrators/jispt.h>
#include <integrators/directpn.h>

#include <fields/VoxelGridField.h>

#include <util/threadpool.h>
#include <houio/HouGeoIO.h>

#include <math/sph.h>

#include <PNSolution.h>

namespace py = pybind11;


struct GlobalRenderInfo
{
	GlobalRenderInfo():
		debug_pixel(-1, -1)
	{

	}

	// always used ---
	const Scene* scene;
	Image* image;
	Box2i crop_window;

	// used for debugging
	P2i debug_pixel;
	std::vector<RNGd> *debug_rngs;
};

struct GlobalComputeFluenceInfo
{
	GlobalComputeFluenceInfo():
		do_direct_light(true),
		do_indirect_light(true)
	{
	}

	bool do_direct_light;
	bool do_indirect_light;
	const Scene* scene;

	std::vector<P3d>* sensor_points;
	std::vector<V3d>* fluence;
};

struct GlobalIntegratePNSolutionInfo
{
	GlobalIntegratePNSolutionInfo()
	{
	}

	const PNSolution* pns;

	std::vector<P3d>* sensor_points;
	std::vector<V3d>* fluence;
};

struct GlobalComputeUncollidedLightInfo
{
	GlobalComputeUncollidedLightInfo()
	{
	}

	const Scene* scene;
	VoxelGrid3d* grid;
};


double luminance( const V3d& c )
{
	return c[0] * 0.212671 + c[1] * 0.715160 + c[2] * 0.072169;
}


void compute_fluence_thread( MonteCarloTaskInfo& ti, const GlobalComputeFluenceInfo* gi )
{
	int numPoints = gi->sensor_points->size();

	// note that y goes from bottom=0 to top=max (rasterspace)
	for (int point_index = ti.taskid; point_index < numPoints; point_index += ti.numTasks)
	{
		V3d f(0.0, 0.0, 0.0);

		bool debug = false;

		P3d pWS = (*gi->sensor_points)[point_index];

		// do raycast ---
		try
		{
			// direct light ---
			if( gi->do_direct_light )
			{
				// lightsample is used to sample direct light
				LightSample ls;
				ls.refP = pWS;
				f += gi->scene->sample_attenuated_directlight( ls, ti.rng, debug );
			}

			// indirect light ---
			if( gi->do_indirect_light )
			{
				V3d direction = sampleSphere(ti.rng);
				double direction_pdf = sampleSpherePDF();

				Ray3d rayWS( pWS, direction );

				RadianceQuery rq;
				rq.ray = rayWS;
				rq.pixel = V2i(-1, -1);
				rq.debug = debug;
				if(gi->scene->volume->getBound().contains(pWS))
					rq.volume = gi->scene->volume;
				f += gi->scene->integrator->Li(gi->scene, rq, ti.rng)/direction_pdf;
			}
		}
		catch (std::exception& e)
		{
			std::cout << "render_volume: caught exception at task=" << ti.taskid << " point_index=" <<  point_index << " sample=" << ti.samples << std::endl;
			std::cout << e.what() << '\n';
			std::flush(std::cout);
			throw e;
		}


		if( std::isnan(luminance(f)) )
			std::cout << "PathTracingTask::run: got NaN value @ index=" << point_index << " sample=" << ti.samples << std::endl;

		// update pixel color
		V3d& c = (*gi->fluence)[point_index];
		c += (f - c)/double(ti.samples+1);

		if(debug)
		{
			c = V3d(1.0, 0.0, 0.0);
			//c_transmittance = Color3f(1.0f, 0.0f, 0.0f);
		}
	} // point


	++ti.samples;
}



Eigen::VectorXd compute_fluence( Volume::Ptr volume,
								 Light::Ptr light,
								 Integrator::Ptr integrator,
								 const Eigen::MatrixXd& points,
								 int seed,
								 bool do_direct_light,
								 bool do_indirect_light,
								 int numSamples)
{
	std::cout << "computing fluence:\n";
	std::cout << volume->toString();
	std::cout << light->toString();
	std::cout << integrator->toString();

	int numPoints = points.rows();

	std::vector<P3d> sensor_points;
	std::vector<V3d> fluence(numPoints, V3d(0.0, 0.0, 0.0));

	for( int i=0;i<numPoints;++i )
	{
		P3d p(points.coeffRef(i, 0),
			  points.coeffRef(i, 1),
			  points.coeffRef(i, 2));
		sensor_points.push_back( p );
	}

	Scene scene;
	scene.camera = 0;
	scene.integrator = integrator.get();
	scene.light = light.get();
	scene.volume = volume.get();

	GlobalComputeFluenceInfo gi;
	gi.scene = &scene;
	gi.fluence = &fluence;
	gi.sensor_points = &sensor_points;
	gi.do_direct_light = do_direct_light;
	gi.do_indirect_light = do_indirect_light;

	MonteCarloTaskInfo::g_seed = seed;

	Terminator terminator(numSamples);
	std::cout << "computing fluence..."<< std::endl;std::flush(std::cout);
	runGenericTasks<MonteCarloTaskInfo, GlobalComputeFluenceInfo>( compute_fluence_thread, &gi, terminator, ThreadPool::getNumSystemCores() );

	Eigen::VectorXd result( fluence.size() );
	for( int i=0;i<numPoints;++i )
		result.coeffRef(i) = luminance(fluence[i]);
	return result;
}



void render_thread( MonteCarloTaskInfo& ti, const GlobalRenderInfo* gi )
{
	//RNGd& rng = ti.rng;
	int width = gi->scene->camera->getResolutionX();
	int height = gi->scene->camera->getResolutionY();

	// note that y goes from bottom=0 to top=max (rasterspace)
	for (int y = ti.taskid; y < gi->crop_window.max.y(); y += ti.numTasks)
	{
		if( y < gi->crop_window.min.y() )
			continue;

		for (int x = gi->crop_window.min.x(); x < gi->crop_window.max.x(); ++x)
		{
			int index = y*width + x;
			//RNGd& rng = (*gi->debug_rngs)[ index ];
			RNGd& rng = ti.rng;
			V3d f(0.0, 0.0, 0.0);
			//Color3f T(1.0f);

			bool debug = false;

			if( gi->debug_pixel.x()>0 )
			{
				debug = (x == gi->debug_pixel.x()) &&
						(y == gi->debug_pixel.y());

				if(!debug)
					continue;
			}

			Ray3d rayWS;
			gi->scene->camera->sampleRay( P2d(x+0.5, y+0.5), rayWS, debug );
			//ti.g.scene->camera->sampleRay( P2d(x+rng.next1D(), y+rng.next1D()), rayWS );

			if(debug)
			{
				std::cout << "rayWS=" << rayWS.toString() << std::endl;
				std::cout << "cam2world=\n" << gi->scene->camera->m_cameraToWorld.getMatrix() << std::endl;
			}

			// do raycast ---
			try
			{
				RadianceQuery rq;
				rq.ray = rayWS;
				rq.pixel = V2i(x, y);
				rq.debug = debug;
				f = gi->scene->integrator->Li(gi->scene, rq, rng);
				//T = rq.transmittance;

				if( debug )
					std::cout << "f=" << f << std::endl;
			}
			catch (std::exception& e)
			{
				std::cout << "render_volume: caught exception at task=" << ti.taskid << " x=" << x << " y=" << y << " " << " index=" <<  index << " sample=" << ti.samples << std::endl;
				std::cout << e.what() << '\n';
				std::flush(std::cout);
				throw e;
			}


			if( std::isnan(luminance(f)) )
				std::cout << "PathTracingTask::run: got NaN value @ index=" << index << " sample=" << ti.samples << std::endl;

			// update pixel color
			Eigen::Vector3d& c = gi->image->pixel(x, y);
			c += (f - c)/double(ti.samples+1);

			if(debug)
			{
				//c = V3d(1.0, 0.0, 0.0);
			}
		} // pixel x
	} // pixel y



	++ti.samples;
}


Image::Ptr render( Volume::Ptr volume, Light::Ptr light, Camera::Ptr camera, Integrator::Ptr integrator, int numSamples )
{
	std::cout << "rendering:\n";

	std::cout << volume->toString();
	std::cout << light->toString();
	std::cout << camera->toString();
	std::cout << integrator->toString();

	Image::Ptr result = std::make_shared<Image>(camera->getResolution());

	Scene scene;
	scene.camera = camera.get();
	scene.integrator = integrator.get();
	scene.light = light.get();
	scene.volume = volume.get();

	GlobalRenderInfo gi;
	gi.scene = &scene;
	gi.image = result.get();
	gi.crop_window = Box2i( V2i(0,0), camera->getResolution() );
	//gi.debug_pixel = V2i(378, 110);
	//gi.debug_pixel = V2i(256, 256);

	/*
	// for debugging, intitialize rngs
	std::vector<RNGd> debug_rngs;
	int index = 0;
	for( int i=0;i<camera->getResolution()[0];++i )
		for( int j=0;j<camera->getResolution()[1];++j, ++index )
			debug_rngs.push_back( RNGd(123+index) );
	gi.debug_rngs = &debug_rngs;
	*/

	//int numSamples = 1;

	/*
	{
		P3d pWS(0.229053, 0.408338, 0.433446);
		volume->evalExtinction(pWS, true);
	}
	*/



	///*
	Terminator terminator(numSamples);
	//Terminator terminator(90.0);
	std::cout << "rendering image..."<< std::endl;std::flush(std::cout);
	runGenericTasks<MonteCarloTaskInfo, GlobalRenderInfo>( render_thread, &gi, terminator, ThreadPool::getNumSystemCores() );
	//*/


	return result;
}



void compute_fluence_pnsolution_thread( MonteCarloTaskInfo& ti, const GlobalIntegratePNSolutionInfo* gi )
{
	int numPoints = gi->sensor_points->size();

	// note that y goes from bottom=0 to top=max (rasterspace)
	for (int point_index = ti.taskid; point_index < numPoints; point_index += ti.numTasks)
	{
		V3d f(0.0, 0.0, 0.0);

		bool debug = false;

		P3d pWS = (*gi->sensor_points)[point_index];
		V3d direction = sampleSphere(ti.rng);
		double direction_pdf = sampleSpherePDF();

		Ray3d rayWS( pWS, direction );

		f = gi->pns->eval(pWS, direction)/direction_pdf;

		if( std::isnan(luminance(f)) )
			std::cout << "PathTracingTask::run: got NaN value @ index=" << point_index << " sample=" << ti.samples << std::endl;

		// update pixel color
		V3d& c = (*gi->fluence)[point_index];
		c += (f - c)/double(ti.samples+1);

		if(debug)
		{
			c = V3d(1.0, 0.0, 0.0);
			//c_transmittance = Color3f(1.0f, 0.0f, 0.0f);
		}
	} // point


	++ti.samples;
}

Eigen::VectorXd compute_fluence_pnsolution( PNSolution::Ptr pns, const Eigen::MatrixXd& points )
{
	std::cout << "computing fluence from pnsolution:\n";

	int numPoints = points.rows();

	std::vector<P3d> sensor_points;
	std::vector<V3d> fluence(numPoints, V3d(0.0, 0.0, 0.0));

	for( int i=0;i<numPoints;++i )
	{
		P3d p(points.coeffRef(i, 0),
			  points.coeffRef(i, 1),
			  points.coeffRef(i, 2));
		sensor_points.push_back( p );
	}

	GlobalIntegratePNSolutionInfo gi;
	gi.pns = pns.get();
	gi.fluence = &fluence;
	gi.sensor_points = &sensor_points;

	int numSamples = 10000;

	Terminator terminator(numSamples);
	std::cout << "computing fluence from pnsolution..."<< std::endl;std::flush(std::cout);
	runGenericTasks<MonteCarloTaskInfo, GlobalIntegratePNSolutionInfo>( compute_fluence_pnsolution_thread, &gi, terminator, ThreadPool::getNumSystemCores() );

	Eigen::VectorXd result( fluence.size() );
	for( int i=0;i<numPoints;++i )
		result.coeffRef(i) = luminance(fluence[i]);
	return result;
}




void compute_unscattered_fluence_thread( MonteCarloTaskInfo& taskinfo, const GlobalComputeUncollidedLightInfo* gi )
{
	V3i res = gi->grid->getResolution();
	int numPoints = res[0]*res[1]*res[2];

	// note that y goes from bottom=0 to top=max (rasterspace)
	for (int point_index = taskinfo.taskid; point_index < numPoints; point_index += taskinfo.numTasks)
	{
		bool debug = false;

		// retrieve worldspace position of current voxel
		int i,j,k;
		gi->grid->getCoordFromIndex( point_index, i, j, k );
		P3d pLS = gi->grid->voxelToLocal( V3d(i+0.5, j+0.5, k+0.5) );
		P3d pWS = gi->scene->volume->localToWorld(pLS);

		// sample light source and transmittance ---
		LightSample ls;
		ls.refP = pWS;
		V3d f = gi->scene->sample_attenuated_directlight( ls, taskinfo.rng );

		// nee gives the uncollided light multiplied by phase function and sigma_s
		//V3d f = nee(ti, taskinfo.rng);

		if( std::isnan(luminance(f)) )
			std::cout << "compute_unscattered_light_thread: got NaN value @ index=" << point_index << " sample=" << taskinfo.samples << std::endl;

		// update pixel color
		V3d& c = gi->grid->lvalue(i, j, k);
		c += (f - c)/double(taskinfo.samples+1);
		//c = V3d(double(i)/double(res[0]));
	} // point


	++taskinfo.samples;
}




//
// computes the emission field from unscattered light
// the unscattered light is estimated and multiplied by sigma_s*phase
// TODO: generalize this to arbitrary phase functions
//
VoxelGridField3d::Ptr compute_unscattered_fluence(  Volume::Ptr volume,
													Light::Ptr light,
													int numSamples,
													Eigen::Vector3i resolution)
{
	std::cout << "compute_unscattered_light:\n";
	std::cout << volume->toString();
	std::cout << light->toString();
	//std::cout << integrator->toString();

	// we need an integrator only for computing the volume attenuation
	Integrator::Ptr integrator = std::make_shared<SimplePT>(-1, true);

	VoxelGridField3d::Ptr result = std::make_shared<VoxelGridField3d>( V3i(resolution.x(), resolution.y(), resolution.z()), (V3d*)0);

	Scene scene;
	scene.integrator = integrator.get();
	scene.light = light.get();
	scene.volume = volume.get();

	GlobalComputeUncollidedLightInfo gi;
	gi.scene = &scene;
	gi.grid = &result->getVoxelGrid();

	Terminator terminator(numSamples);
	std::cout << "computing uncollided light..."<< std::endl;std::flush(std::cout);
	runGenericTasks<MonteCarloTaskInfo, GlobalComputeUncollidedLightInfo>( compute_unscattered_fluence_thread, &gi, terminator, ThreadPool::getNumSystemCores() );

	return result;
}
struct GlobalBakeUnscatteredDirectionalLightInfo
{
	Volume* volume;
	VoxelGrid3d* grid; // result where all threads write to
	V3d lightDirection;
	double lightRadiance;
	double maxDistance;
};


void bake_unscattered_directional_light_thread( MonteCarloTaskInfo& taskinfo, const GlobalBakeUnscatteredDirectionalLightInfo* gi )
{
	RNGd& rng = taskinfo.rng;
	V3i res = gi->grid->getResolution();
	int numPoints = res[0]*res[1]*res[2];

	// note that y goes from bottom=0 to top=max (rasterspace)
	for (int point_index = taskinfo.taskid; point_index < numPoints; point_index += taskinfo.numTasks)
	{
		bool debug = false;

		// retrieve worldspace position of current voxel
		int i,j,k;
		gi->grid->getCoordFromIndex( point_index, i, j, k );
		P3d pLS = gi->grid->voxelToLocal( V3d(i+0.5, j+0.5, k+0.5) );
		P3d pWS = gi->volume->localToWorld(pLS);

		// compute directional light contribution to current voxel  -----
		V3d f(0.0);
		if( gi->volume->getBound().contains(pWS) )
		{
			// find the maximum distance for attenuation computation ---
			double maxt = gi->maxDistance;
			{
				double mint;
				bool result = gi->volume->intersectBound(Ray3d( pWS, -gi->lightDirection ), mint, maxt);
				if( !result )
					// thats strange
					throw std::runtime_error("expected intersection");
			}

			//NB: the geometry term between light source position and current vertex is already been dealt with during light sampling

			// sample transmittance ---
			int component = 0;
			V3d sigma_t;
			double distance = delta_tracking( gi->volume, Ray3d( pWS, -gi->lightDirection ), maxt, component, rng, sigma_t, debug );

			if( distance < maxt )
				// since we set the contribution to 0, we dont need to produce the pdf in the denominator
				f = V3d(0.0);
			else
				// transmittance sampling pdf cancels out with the transmittance term
				// TODO: we need to divide by extinction, dont we?
				// because we do not apply scattering and therefore we dont have
				// that sigma_s/sigma_t = albedo thing
				f = gi->lightRadiance;
		}

		// nee gives the uncollided light multiplied by phase function and sigma_s
		//V3d f = nee(ti, taskinfo.rng);

		if( std::isnan(luminance(f)) )
			std::cout << "compute_unscattered_light_thread: got NaN value @ index=" << point_index << " sample=" << taskinfo.samples << std::endl;

		// update pixel color
		V3d& c = gi->grid->lvalue(i, j, k);
		c += (f - c)/double(taskinfo.samples+1);
		//c = V3d(double(i)/double(res[0]));
	} // point


	++taskinfo.samples;
}


VoxelGridField3d::Ptr bake_unscattered_directional_light(  Volume::Ptr volume,
														   Eigen::Vector3d direction, // light direction
														   int numSamples,
														   Eigen::Vector3i resolution)
{
	std::cout << "bake_unscattered_directional_light:\n";
	std::cout << volume->toString();

	// we need an integrator only for computing the volume attenuation
	//Integrator::Ptr integrator = std::make_shared<SimplePT>(-1, true);

	VoxelGridField3d::Ptr result = std::make_shared<VoxelGridField3d>( V3i(resolution.x(), resolution.y(), resolution.z()), (V3d*)0);

	//Scene scene;
	//scene.integrator = integrator.get();
	//scene.light = light.get();
	//scene.volume = volume.get();

	GlobalBakeUnscatteredDirectionalLightInfo gi;
	gi.volume = volume.get();
	gi.grid = &result->getVoxelGrid();
	gi.lightDirection = direction;
	gi.lightRadiance = 1.0;
	gi.maxDistance = volume->getBound().getExtents().norm()*1.1;

	Terminator terminator(numSamples);
	std::cout << "baking unscattered directional light..."<< std::endl;std::flush(std::cout);
	runGenericTasks<MonteCarloTaskInfo, GlobalBakeUnscatteredDirectionalLightInfo>( bake_unscattered_directional_light_thread, &gi, terminator, ThreadPool::getNumSystemCores() );

	return result;
}


VoxelGridField3d::Ptr resample( VoxelGridField3d::Ptr field, const Eigen::Vector3i& res )
{
	V3i in_res = field->getResolution();
	VoxelGridField3d::Ptr result = std::make_shared<VoxelGridField3d>( V3i(res.x(), res.y(), res.z()), (V3d*)0);

	// iterate all new voxels
	for( int i=0;i<res[0];++i )
		for( int j=0;j<res[1];++j )
			for( int k=0;k<res[2];++k )
			{
				// compute bounding box of new voxel in field voxelspace
				V3d minVS = field->localToVoxel(result->voxelToLocal( V3d(double(i), double(j), double(k)) ) );
				V3d maxVS = field->localToVoxel(result->voxelToLocal( V3d(double(i+1), double(j+1), double(k+1)) ) );

				// get min and max voxels of the new voxel in field
				V3i vmin( std::max(int(floor(minVS.x())),0), std::max(int(floor(minVS.y())),0), std::max(int(floor(minVS.z())),0) );
				V3i vmax( std::min(int(ceil(maxVS.x())),in_res.x()-1), std::min(int(ceil(maxVS.y())),in_res.y()-1), std::min(int(ceil(maxVS.z())),in_res.z()-1) );


				// iterate each field voxel we touch
				V3d value(0.0, 0.0, 0.0);
				double w_sum = 0.0; // sum of all voxel volumes which intersect

				//printf( "%d\n", vmin.x );
				//printf( "%d\n", vmax.x );
				//printf( "%d    %d    %d\n", vmax.x-vmin.x, vmax.y-vmin.y, vmax.z-vmin.z );
				//printf( "%f   %f   %f\n", minVS.x, minVS.y, minVS.z );

				for( int k_src=vmin.z(); k_src<=vmax.z();++k_src )
				{
					for( int j_src=vmin.y(); j_src<=vmax.y();++j_src )
					{
						for( int i_src=vmin.x(); i_src<=vmax.x();++i_src )
						{
							// compute intersection volume (http://stackoverflow.com/a/5556796)
							double w =   std::max(std::min(double(i_src+1),maxVS.x())-std::max(double(i_src),minVS.x()),0.0)
										*std::max(std::min(double(j_src+1),maxVS.y())-std::max(double(j_src),minVS.y()),0.0)
										*std::max(std::min(double(k_src+1),maxVS.z())-std::max(double(k_src),minVS.z()),0.0);
							V3d v = field->getVoxelGrid().sample(i_src, j_src, k_src);
							value += v*w;
							w_sum += w;
						}
					}
				}
				if(w_sum > 0.0)
					result->getVoxelGrid().lvalue(i, j, k) = value/w_sum;
			}


	return result;
}


Light::Ptr create_point_light( const Eigen::Vector3d& position, const Eigen::Vector3d& power )
{
	return std::make_shared<PointLight>(position, power);
}

Light::Ptr create_directional_light(const Eigen::Vector3d& direction, const Eigen::Vector3d& radiance)
{
	return std::make_shared<DirectionalLight>(direction, radiance);
}


Camera::Ptr create_perspective_camera( int xres, int yres, double hfov )
{
	return std::make_shared<PerspectiveCamera>(xres, yres, hfov);
}


Camera::Ptr create_sphere_camera( const Eigen::Vector3d& pWS, int xres, int yres )
{
	return std::make_shared<SphereCamera>( pWS, xres, yres );
}

Integrator::Ptr create_directpn_integrator( PNSolution::Ptr pns, bool doSingleScattering, bool doMultipleScattering )
{
	return std::make_shared<DirectPN>(pns.get(), doSingleScattering, doMultipleScattering);
}


Integrator::Ptr create_simplept_integrator( bool doSingleScattering, int maxDepth )
{
	int md = std::numeric_limits<int>::max();
	if(maxDepth >= 0)
		md = maxDepth;
	return std::make_shared<SimplePT>(md, doSingleScattering);
}

Integrator::Ptr create_pnispt_integrator( PNSolution::Ptr pns, bool doSingleScattering, int maxDepth )
{
	int md = std::numeric_limits<int>::max();
	if(maxDepth >= 0)
		md = maxDepth;
	return std::make_shared<PNISPT>(pns, doSingleScattering, md);
}

Integrator::Ptr create_jispt_integrator()
{
	int maxDepth = std::numeric_limits<int>::max();
	bool doSingleScattering = true;
	return std::make_shared<JISPT>(maxDepth, doSingleScattering);
}


Volume::Ptr create_volume()
{
	return std::make_shared<Volume>();
}

Field3d::Ptr create_constant_field3d( const Eigen::Vector3d& value )
{
	return std::make_shared<ConstantField3d>(value);
}


std::tuple<Field3d::Ptr, Eigen::Vector3d, Eigen::Vector3d> load_bgeo( const std::string& filename )
{
	houio::ScalarField::Ptr sigma_t = houio::HouGeoIO::importVolume(filename);

	if(!sigma_t)
		throw std::runtime_error("load_bgeo");

	houio::math::V3i res = sigma_t->getResolution();
	int numVoxels = res.x*res.y*res.z;

	std::vector<V3d> data(numVoxels);
	for( int i=0;i<res.x;++i )
		for( int j=0;j<res.y;++j )
			for( int k=0;k<res.z;++k )
			{
				float value = sigma_t->sample(i, j, k);
				// numpy indexing (k ist fastest)
				// TODO: figure out, why we have to flip i and k indices...
				// ANSWER: something has changed within houdini in regards to the transformation of volumes
				// the old nebulae dataset had its x and z axis flipped. This was fixed manually on the dataset.
				//int index_dst = (res.x-i-1)*res.y*res.z + j*res.z + (res.z-k-1);

				//
				int index_dst = i*res.y*res.z + j*res.z + k;
				/*
				if(value>0.0)
				{
					std::cout << "i=" << i << " j=" << j << " k=" << k << " value=" << value << std::endl;
					std::cout << "index_src=" << (k*res.x*res.y + j*res.x + i) << std::endl;
					std::cout << "index_dst=" << index_dst << std::endl;
				}
				*/


				data[index_dst] = V3d(value, value, value);
			}

	/*
	for( int i=0;i<numVoxels;++i )
	{
		if( data[i].x() > 0.0 )
		{
			std::cout << "data[" << i*3+0 << "]=" << data[i].x() << std::endl;
			std::cout << "data[" << i*3+1 << "]=" << data[i].y() << std::endl;
			std::cout << "data[" << i*3+2 << "]=" << data[i].z() << std::endl;
		}
	}
	*/


	VoxelGridField3d::Ptr result = std::make_shared<VoxelGridField3d>(  V3i(res.x, res.y, res.z),
																		data.data());
	Eigen::Vector3d bound_min(sigma_t->bound().minPoint.x, sigma_t->bound().minPoint.y, sigma_t->bound().minPoint.z);
	Eigen::Vector3d bound_max(sigma_t->bound().maxPoint.x, sigma_t->bound().maxPoint.y, sigma_t->bound().maxPoint.z);
	return std::make_tuple( result, bound_min, bound_max );
}

void save_bgeo( const std::string& filename, VoxelGridField3d::Ptr field, Eigen::Vector3d bound_min, Eigen::Vector3d bound_max )
{
	V3i res = field->getResolution();
	houio::ScalarField::Ptr sigma_t = houio::ScalarField::create( houio::math::V3i(res[0], res[1], res[2]),
																  houio::math::Box3f( houio::math::V3f(bound_min[0], bound_min[1], bound_min[2]),
																					  houio::math::V3f(bound_max[0], bound_max[1], bound_max[2])));

	int component = 0;
	for( int i=0;i<res.x();++i )
		for( int j=0;j<res.y();++j )
			for( int k=0;k<res.z();++k )
			{
				// TODO: figure out, why we have to flip i and k indices...
				//sigma_t->lvalue(i, j, k) = field->getVoxelGrid().lvalue(res.x()-i-1, j, res.z()-k-1)[component];
				sigma_t->lvalue(i, j, k) = field->getVoxelGrid().lvalue(i, j, k)[component];
			}

	houio::HouGeoIO::xport( filename, sigma_t );
}

PNSolution::Ptr load_pnsolution( const std::string& filename )
{
	return std::make_shared<PNSolution>(filename);
}


Image::Ptr load_image( const std::string& filename )
{
	return std::make_shared<Image>(filename);
}

void save_exr( const std::string& filename, const Eigen::MatrixXd& values )
{
	Image::Ptr img = std::make_shared<Image>( V2i(values.rows(), values.cols()) );
	for( int i=0;i<values.rows();++i )
		for( int j=0;j<values.cols();++j )
		{
			double value = values.coeffRef(i, j);
			img->pixel(j, i) = V3d(value, value, value);
		}
	img->save(filename);
}

ImageSampler::Ptr create_image_sampler( Image::Ptr image )
{
	return std::make_shared<ImageSampler>(image);
}

#include <math/shvector.h>

#include <mitsuba/core/shvector.h>

struct SHTest
{
	typedef std::shared_ptr<SHTest> Ptr;

	SHTest( int order, int depth ):
		m_order(order),
		shsampler(new SHSampler(order+1, depth)),
		shsampler2(new PNSolution::SHSampler(order+1, depth)),
		shv(order+1),
		mitsuba_shvector(order+1),
		mitsuba_shsampler(new mitsuba::SHSampler(order+1, depth)),
		rng(123)
	{
		SHVector::staticInitialization();
		mitsuba::SHVector::staticInitialization();
		sph::staticInit();

		/*
		shv(0,0) = 1.23;
		shv(1,-1) = 1.0;
		shv(1,0) = .42;
		shv(1,1) = 0.1;
		*/
		for( int l=0;l<=order;++l )
			for( int m=-l;m<=l;++m )
			{
				double value = rng.next1D()*2.0-1.0;
				//std::cout << "SHTest::SHTest f(" << l << ", " << m << ") = " << value << std::endl;
				shv(l,m) = value;
				mitsuba_shvector(l,m) = value;
			}


		//shv(1,-1) = 1.0;
		//shv(1,0) = 3.42;
		//shv(1,1) = -0.86;
	}

	void setCoeffs( const Eigen::VectorXd& coeffs )
	{
		//if( coeffs.rows() != shv. )
		int index = 0;
		for( int l=0;l<=m_order;++l )
			for( int m=-l;m<=l;++m,++index )
			{
				double value = coeffs.coeffRef(index);
				shv(l,m) = value;
				mitsuba_shvector(l,m) = value;
			}
	}

	double eval( double theta, double phi )
	{
		//V3d t = sphericalDirection(theta, phi);
		//return shv.eval(t);
		return sph::eval( theta, phi, shv.m_coeffs.data(), m_order );
		//return mitsuba_shvector.eval(theta, phi);
	}

	//Eigen::Vector3d sample()
	std::tuple<double, double, double, double, double> sample()
	{
		double r1 = rng.next1D();
		double r2 = rng.next1D();

		P2d sample(r1, r2);
		//double pdf = shsampler->warp(shv, sample);
		double pdf = shsampler2->warp(shv.m_coeffs.data(), sample);

		//mitsuba::Point2 sample2(r1, r2);
		//double pdf2 = mitsuba_shsampler->warp(mitsuba_shvector, sample2);
		/*
		//V3d test = sphericalDirection(sample[0], sample[1]);
		//return Eigen::Vector3d(test[0], test[1], test[2]);
		return std::make_pair(sample[0], sample[1]);
		*/
		return std::make_tuple(sample[0], sample[1], pdf, r1, r2);
		//return std::make_tuple(sample2[0], sample2[1], pdf2, r1, r2);
	}

	//Eigen::Vector3d sample()
	std::tuple<double, double, double, double, double> sample2( double r1, double r2)
	{
		//double r1 = rng.next1D();
		//double r2 = rng.next1D();

		P2d sample(r1, r2);
		double pdf = shsampler->warp(shv, sample);

		return std::make_tuple(sample[0], sample[1], pdf, r1, r2);
	}

	double pdf(double r1, double r2)
	{
		//P2d sample(r1, r2);
		//double pdf = shsampler->warp(shv, sample);

		mitsuba::Point2 sample2(r1, r2);
		double pdf2 = mitsuba_shsampler->warp(mitsuba_shvector, sample2);
		/*
		//V3d test = sphericalDirection(sample[0], sample[1]);
		//return Eigen::Vector3d(test[0], test[1], test[2]);
		return std::make_pair(sample[0], sample[1]);
		*/
		//return std::make_tuple(sample[0], sample[1], pdf);
		return pdf2;
	}

	Eigen::MatrixXd get_blocks( int depth, bool use_org )
	{
		return shsampler->getBlocks(depth, shv, use_org);

	}

	SHVector shv;
	SHSampler* shsampler;
	PNSolution::SHSampler* shsampler2;
	RNGd rng;

	mitsuba::SHVector mitsuba_shvector;
	mitsuba::SHSampler* mitsuba_shsampler;

	int m_order;

};



SHTest::Ptr create_shtest( int order, int depth )
{
	return std::make_shared<SHTest>( order, depth );
}


VoxelGridField3d::Ptr load_voxelgridfield3d(const std::string& filename)
{
	return std::make_shared<VoxelGridField3d>(filename);
}


template<typename T>
using SphericalFunction = std::function<T (double, double)>;
// resolution.x -> resolution along theta angle
// resolution.y -> resolution along phi angle
template<typename T>
T sample_spherical_function( SphericalFunction<T> f, V2i resolution, std::vector<houio::math::V3f>& points )
{
	T result = 0.0;

	double min_theta = 0;
	double max_theta = M_PI;
	double min_phi = 0;
	double max_phi = 2.0*M_PI;

	int resolution_theta=resolution.x(); // resolution along theta angle
	int resolution_phi=resolution.y(); // resolution along phi angle
	double dtheta = (max_theta-min_theta)/double(resolution_theta);
	double dphi = (max_phi-min_phi)/double(resolution_phi);

	double pixel_area = dtheta*dphi;

	for( int t=0;t<resolution_theta;++t )
		for( int p=0;p<resolution_phi;++p )
		{
			double phi = dphi*p;
			double theta = dtheta*t;
			double f_value = f(theta, phi);
			result+=f_value*pixel_area*std::sin(theta);

			V3d pos = sphericalDirection(theta, phi)*f_value;
			points.push_back( houio::math::V3f(pos.x(), pos.y(), pos.z()) );
		}

	return result;
}

void test()
{
	/*
	sph::staticInit();
	RNGd rng(123);
	double result = 0.0;
	int numSamples = 100;
	for( int i=0;i<numSamples;++i )
	{
		V3d d = sampleSphere(rng);
		double pdf = sampleSpherePDF();
		double theta, phi;
		sphericalCoordinates(d, theta, phi);

		double coeff = 1.0;
		double f = sph::eval(theta, phi, &coeff, 0)/pdf;
		std::cout << sph::eval(theta, phi, &coeff, 0) << std::endl;

		result += (f-result)/double(i+1);
	}

	std::cout << "test result=" << result << std::endl;
	*/

	/*

	// add constant factor

	int res = 16;
	double zStep = -2 / (double) res;
	double phiStep = 2 * (double) M_PI / (double) res;

	//double phi_integral = phiStep;

	double scale = 0.001;
	double sum = 0.0;
	for( int phiBlock=0;phiBlock<res;++phiBlock )
	{
		double phi_integral = phiStep;
		for( int zBlock=0;zBlock<res;++zBlock )
		{
			double theta_max = 1.0+zBlock*zStep;
			double theta_min = 1.0+(zBlock+1)*zStep;
			double theta_integral = -theta_min + theta_max;
			sum += theta_integral*phi_integral*scale;
		}
	}
	std::cout << "sum=" << sum << std::endl;
	*/


	V3d wi(0.0, 0.0, 1.0);
	V2i sampling_res(128, 256);
	double g = 0.6;
	int order = 3;

	HGPhase hg(g);


	// hg groundtruth
	{
		std::function<double(double, double)> hg_fun = [&]( double theta, double phi )
		{
			V3d wo = sphericalDirection(theta, phi);
			return hg.eval(wi, wo);
		};
		std::vector<houio::math::V3f>  hg_gt_points;

		sample_spherical_function( hg_fun, sampling_res, hg_gt_points );
		houio::HouGeoIO::xport( "hg_gt.bgeo", hg_gt_points );
	}

	// hg sampling evaluation
	{
		//std::vector<houio::math::V3f>  hg_samples_scaled;
		std::vector<houio::math::V3f>  hg_samples;
		std::vector<houio::math::V3f>  hg_samples_values;
		RNGd rng;
		int numSamples = 100000;
		for( int i=0;i<numSamples;++i )
		{
			V3d wo;
			double pdf;
			double throughput_over_pdf = hg.sample(wi, wo, pdf, rng);
			double throughput = hg.eval(wi, wo);

			//hg_samples_scaled.push_back( houio::math::V3d(wo.x()*throughput, wo.y()*scale, wo.z()*scale) );
			hg_samples.push_back( houio::math::V3f(wo.x(), wo.y(), wo.z()) );
			hg_samples_values.push_back(houio::math::V3f(throughput, pdf, throughput_over_pdf));
		}
		std::map<std::string,std::vector<houio::math::V3f>> pattrs;
		pattrs["P"] = hg_samples;
		pattrs["eval"] = hg_samples_values;
		houio::HouGeoIO::xport( "hg_samples.bgeo", pattrs );
		//houio::HouGeoIO::xport( "hg_samples.bgeo", hg_samples );
		//houio::HouGeoIO::xport( "hg_samples_scaled.bgeo", hg_samples_scaled );
	}

	// evaluation of sh reconstruction for HG
	{
		int order = 10;
		std::function<double(double, double)> hg_fun = [&]( double theta, double phi )
		{
			//V3d wo = sphericalDirection(theta, phi);
			//return hg.eval(wi, wo);
			double sum = 0.0;
			for( int l=0;l<=order;++l )
				for( int m=-l;m<=l;++m )
					sum += hg.shCoeff(l,m)*sph::basis(l,m, theta, phi);
			return sum;
		};
		std::vector<houio::math::V3f>  hg_sh_points;

		sample_spherical_function( hg_fun, sampling_res, hg_sh_points );
		houio::HouGeoIO::xport( "hg_sh.bgeo", hg_sh_points );
	}


}

/*
std::complex<double> sh_basis( int l, int m, double theta, double phi )
{
	//double theta, phi;
	//sphericalCoordinates(V3d(direction), theta, phi);
	return sph::complex_basis(l,m, theta, phi);
}

std::complex<double> sh_basis_conj( int l, int m, double theta, double phi)
{
	//double theta, phi;
	//sphericalCoordinates(V3d(direction), theta, phi);
	return std::conj(sph::complex_basis(l,m, theta, phi));
}
*/

double sh_basis( int l, int m, double theta, double phi )
{
	return sph::basis(l,m, theta, phi);
}

double sh_basis_conj( int l, int m, double theta, double phi)
{
	return std::conj(sph::basis(l,m, theta, phi));
}

PYBIND11_MODULE(renderer, m)
{
	sph::staticInit();

	// Light ============================================================
	py::class_<Light, Light::Ptr> class_light(m, "Light");
	class_light
	;

	// Camera ============================================================
	py::class_<Camera, Camera::Ptr> class_camera(m, "Camera");
	class_camera
		.def( "setCameraToWorldTransform",
		[](Camera &m, const Eigen::Matrix4Xd& cameraToWorld )
		{
			m.setCameraToWorld(Transformd(cameraToWorld));
		})
	;

	// Integrators ============================================================
	py::class_<Integrator, Integrator::Ptr> class_integrator(m, "Integrator");
	class_integrator
	;

	// ...
	py::class_<PNISPT, PNISPT::Ptr> class_pnispt(m, "PNISPT", class_integrator);
	class_pnispt
	.def("dbgGet", &PNISPT::dbgGet)
	;



	// Volume ============================================================
	py::class_<Volume, Volume::Ptr> class_volume(m, "Volume");
	class_volume
		.def( "setBound",
		[](Volume &m, const Eigen::Vector3d& min, const Eigen::Vector3d& max )
		{
			m.setBound(Box3d(min, max));
		})
		.def( "setExtinctionAlbedo", &Volume::setExtinctionAlbedo )
		.def( "setPhaseFunction", &Volume::setPhaseFunction )
		.def("worldToLocal",
		[]( Volume &m, const Eigen::Vector3d& pWS )
		{
			Eigen::Vector3d pLS = m.worldToLocal(pWS);
			return pLS;
		})
		.def("localToWorld",
		[]( Volume &m, const Eigen::Vector3d& pLS )
		{
			Eigen::Vector3d pWS = m.localToWorld(pLS);
			return pWS;
		})
		.def( "getBoundMin", []( Volume &m )
		{
			Eigen::Vector3d pWS = m.getBound().min;
			return pWS;
		})
		.def( "getBoundMax", []( Volume &m )
		{
			Eigen::Vector3d pWS = m.getBound().max;
			return pWS;
		})
		.def("getSlice", &Volume::getSlice )
	;

	// Image ============================================================
	py::class_<Image, Image::Ptr> class_image(m, "Image");
	class_image
		.def( "save", &Image::save )
		.def( "asMatrix",
		[]( Image& m )
		{
			int channel = 0;
			V2i res = m.getResolution();
			Eigen::MatrixXd result( res[1], res[0] );

			for( int i=0;i<res[0];++i )
				for( int j=0;j<res[1];++j )
					result.coeffRef(j, i) = m.pixel(i, j)[channel];

			return result;
		})
		.def( "uvToRaster",
		[]( Image& m, const Eigen::Vector2d& uv )
		{
			Eigen::Vector2d pRaster = m.uvToRaster(uv);
			return pRaster;
		})
		.def( "rasterToUV",
		[]( Image& m, const Eigen::Vector2d& pRaster )
		{
			Eigen::Vector2d uv = m.rasterToUV(pRaster);
			return uv;
		})
		.def( "eval",
		[]( Image& m, const Eigen::Vector2d& pRaster )
		{
			Eigen::Vector3d result = m.eval(pRaster);
			return result;
		})

	;

	// ImageSampler ============================================================
	py::class_<ImageSampler, ImageSampler::Ptr> class_imagesampler(m, "ImageSampler");
	class_imagesampler
		.def( "sample",
		[]( ImageSampler& m, double r1, double r2 )
		{
			double pdf=0.0;
			Eigen::Vector2d pRaster = m.sample(pdf, P2d(r1, r2));
			return pRaster;
		})
		/*
		.def( "asMatrix",
		[]( ImageSampler& m )
		{
			Image map( V2i(xres, yres) );
			Eigen::MatrixXd result( res[1], res[0] );
			for( int i=0;i<xres;++i )
				for( int j=0;j<yres;++j )
				{
					double pdf = m_colDpdf[i]*m_rowDpdf[i][j];
					double c = pdf;
					map.pixel(j, i) = V3d(c, c, c);
				}
			map.save( "c:/projects/epfl/epfl17/temp/test_is.exr" );


			int channel = 0;
			V2i res = m.getResolution();
			Eigen::MatrixXd result( res[1], res[0] );

			for( int i=0;i<res[0];++i )
				for( int j=0;j<res[1];++j )
					result.coeffRef(j, i) = m.pixel(i, j)[channel];

			return result;
		})
		*/
	;

	// Field3d ============================================================
	py::class_<Field3d, Field3d::Ptr> class_field3d(m, "Field3d");
	class_field3d
	;
	py::class_<Fieldcd, Fieldcd::Ptr> class_fieldcd(m, "Fieldcd");
	class_fieldcd
	.def("eval",
	 [](Fieldcd &m,
		const Eigen::Matrix<double, 3, 1>& pWS )
	 {
		return m.eval(pWS);
	 })
	;

	// VoxelGrid (complex valued) ============================================================
	py::class_<VoxelGridFieldcd, VoxelGridFieldcd::Ptr> class_VoxelGridFieldcd(m, "VoxelGridFieldcd", class_fieldcd);
	class_VoxelGridFieldcd
	.def("__init__",
	[](VoxelGridFieldcd &m, py::array b)
	{
		typedef Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic> Strides;

		// Request a buffer descriptor from Python
		py::buffer_info info = b.request();

		// Some sanity checks ...
		if (info.format != py::format_descriptor<std::complex<double>>::format())
			throw std::runtime_error("Incompatible format: expected a complex array!");

		if (info.ndim != 3)
			throw std::runtime_error("Incompatible buffer dimension!");

		V3i resolution( int(info.shape[0]),
						int(info.shape[1]),
						int(info.shape[2]) );

		auto data = static_cast<std::complex<double> *>(info.ptr);
		new (&m) VoxelGridFieldcd( resolution, data );
	})
	//.def("test", &VoxelGridField::test)
	//.def("getSlice", &VoxelGridField::getSlice)
	;

	// VoxelGrid (real valued) ============================================================
	py::class_<VoxelGridField3d, VoxelGridField3d::Ptr> class_VoxelGridField3d(m, "VoxelGridField3d", class_field3d);
	class_VoxelGridField3d
		.def("__init__",
		[](VoxelGridField3d &m, py::array b)
		{
			typedef Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic> Strides;

			// Request a buffer descriptor from Python
			py::buffer_info info = b.request();

			// Some sanity checks ...
			if (info.format != py::format_descriptor<double>::format())
				throw std::runtime_error("Incompatible format: expected a double array!");

			if (info.ndim != 4)
				throw std::runtime_error("Incompatible buffer dimension!");

			if (info.shape[3] != 3)
				throw std::runtime_error("Incompatible buffer dimension!");

			V3i resolution( int(info.shape[0]),
							int(info.shape[1]),
							int(info.shape[2]) );

			auto data = static_cast<double*>(info.ptr);
			//for( int i=0;i<resolution.x()*resolution.y()*resolution.z()*3;++i )
			//	std::cout << "data[" << i << "]=" << data[i] << std::endl;

			new (&m) VoxelGridField3d(resolution, (V3d*)data);
		})
		.def("save", &VoxelGridField3d::save)
		.def("asArray",
			[]( VoxelGridField3d &m )
			{
				V3i res = m.getResolution();
				py::array b( py::dtype::of<double>(),
							{res[0], res[1], res[2], 3},
							{int(sizeof(double))*res[2]*res[1]*3, int(sizeof(double))*res[2]*3, int(sizeof(double))*3, int(sizeof(double))},
							m.getData());
				return b;
			})
	;


	// PNSolution ==============================
	py::class_<PNSolution, PNSolution::Ptr> class_pnsolution(m, "PNSolution");
	class_pnsolution
	.def("__init__",
	///*
	[](PNSolution &m,
		int order,
		const Eigen::Matrix<int, 3, 1>& resolution,
		const Eigen::Matrix<double, 3, 1>& bound_min,
		const Eigen::Matrix<double, 3, 1>& bound_max,
		const Eigen::VectorXd& data
		)
	//*/
	/*
	[](PNSolution &m,
		const std::string& filename
	)
	*/
	{
		//int order, const V3i& resolution, const Box3d& bound, const std::complex<double> *data
		new (&m) PNSolution( order, resolution, Box3d(bound_min, bound_max), data.data() );
		//new (&m) PNSolution( filename );
	})
	.def("eval",
	 [](PNSolution &m,
		const Eigen::Matrix<double, 3, 1>& pWS,
		const Eigen::Matrix<double, 3, 1>& direction )
	 {
		return m.eval(pWS, direction);
	 })
	.def("sample",
	 [](PNSolution &m,
		const Eigen::Matrix<double, 3, 1>& pWS,
		const Eigen::Matrix<double, 2, 1>& sample)
	 {
		double pdf=0.0;
		V3d dir = m.sample(pWS, pdf, sample);
		std::pair<Eigen::Vector3d, double> result = std::make_pair(dir, pdf);
		return result;
	 })
	.def("pdf",
	 [](PNSolution &m,
		const Eigen::Matrix<double, 3, 1>& pWS,
		const Eigen::Matrix<double, 3, 1>& direction)
	 {
		return m.pdf(pWS, direction);
	 })
	.def("getBlocks",
	 [](PNSolution &m,
		const Eigen::Matrix<double, 3, 1>& pWS,
		 int depth)
	 {
		return m.getBlocks(pWS, depth);
	 })
	.def("evalCoefficient",
	 [](PNSolution &m,
		const Eigen::Matrix<double, 3, 1>& pWS,
		int coeff_index)
	 {
		return m.evalCoefficient(pWS, coeff_index);
	 })
	.def("evalCoefficients",
	 [](PNSolution &m,
		const Eigen::Matrix<double, 3, 1>& pWS)
	 {
		return m.evalCoefficients(pWS);
	 })
	.def("save", &PNSolution::save)
	.def("getResolution", &PNSolution::getResolution )
	.def("getNumCoeffs", &PNSolution::getNumCoeffs )
	.def("getBoundMin", &PNSolution::getBoundMin )
	.def("getBoundMax", &PNSolution::getBoundMax )
	/*
	.def("voxelToLocal", &PNSolution::voxelToLocal )
	.def("localToVoxel", &PNSolution::localToVoxel )
	*/
	.def("localToWorld",
	[]( PNSolution &m,
		const Eigen::Matrix<double, 3, 1>& pLS)
	{
		return m.localToWorld(pLS);
	})
	/*
	.def("worldToLocal", &PNSolution::worldToLocal )
	.def("voxelToWorld", &PNSolution::voxelToWorld )
	.def("worldToVoxel", &PNSolution::worldToVoxel )
	*/
	//.def("test", &PNSolution::test )
	.def("setFilter", &PNSolution::setFilter )
	.def("getCoefficientField",
		[]( PNSolution &m, int coeff_index )
		{
			int numCoeffs = m.getNumCoeffs();
			V3i res = m.getResolution();
			std::vector<double> data(res[0]*res[1]*res[2]);

			for( int i=0;i<res[0];++i )
				for( int j=0;j<res[1];++j )
					for( int k=0;k<res[2];++k )
					{
						int index_dst = i*res[2]*res[1] + j*res[2] + k;
						int index_src = m.getIndex(V3i(i, j, k))+coeff_index;
						data[index_dst] = m.data()[index_src];
					}
			py::array b( //py::dtype(std::string("float64")),
						py::dtype::of<double>(),
					   {res[0], res[1], res[2]},
					   {int(sizeof(double))*res[2]*res[1], int(sizeof(double))*res[2], int(sizeof(double))},
					   data.data());
			return b;
		})
	;


	// PhaseFunction ============================================================
	py::class_<PhaseFunction, PhaseFunction::Ptr> class_phasefunction(m, "PhaseFunction");
	class_phasefunction
	;

	// HGPhase ============================================================
	py::class_<HGPhase, HGPhase::Ptr> class_hgphase(m, "HGPhase", class_phasefunction);
	class_hgphase
		.def("__init__",
		[](HGPhase &m, double g)
		{
			new (&m) HGPhase(g);
		})
	;

	// SHTest ============================================================
	py::class_<SHTest, SHTest::Ptr> class_shtest(m, "SHTest");
	class_shtest
		.def( "eval", &SHTest::eval )
		.def( "sample", &SHTest::sample )
		.def( "sample2", &SHTest::sample2 )
		.def( "pdf", &SHTest::pdf )
		.def( "setCoeffs", &SHTest::setCoeffs )
		.def( "get_blocks", &SHTest::get_blocks )
	;


	// main render function
	m.def( "render", &render );
	m.def( "compute_fluence", &compute_fluence );
	m.def( "create_shtest", &create_shtest );
	m.def( "compute_unscattered_fluence", &compute_unscattered_fluence );
	m.def( "bake_unscattered_directional_light", &bake_unscattered_directional_light );
	// lights
	m.def( "create_point_light", &create_point_light );
	m.def( "create_directional_light", &create_directional_light );
	// cameras
	m.def( "create_perspective_camera", &create_perspective_camera );
	m.def( "create_sphere_camera", &create_sphere_camera );
	// integrators
	m.def( "create_simplept_integrator", &create_simplept_integrator );
	m.def( "create_pnispt_integrator", &create_pnispt_integrator );
	m.def( "create_jispt_integrator", &create_jispt_integrator );
	m.def( "create_directpn_integrator", &create_directpn_integrator );
	// volumes
	m.def( "create_volume", &create_volume );
	// fields
	m.def( "create_constant_field3d", &create_constant_field3d );
	m.def( "load_bgeo", &load_bgeo );
	m.def( "save_bgeo", &save_bgeo );
	m.def( "load_voxelgridfield3d", &load_voxelgridfield3d );
	m.def( "resample", &resample );
	// pnsolution
	m.def( "load_pnsolution", &load_pnsolution );
	m.def( "compute_fluence_pnsolution", &compute_fluence_pnsolution );
	// image
	m.def( "load_image", &load_image );
	m.def( "save_exr", &save_exr );
	m.def( "create_image_sampler", &create_image_sampler );
	m.def( "blur_image", &blur_image );
	m.def( "ldr_image", &ldr_image );
	// misc
	m.def( "test", &test );
	m.def( "sh_basis", &sh_basis );
	m.def( "sh_basis_conj", &sh_basis_conj );



}
