#include <iostream>
#include <complex>

#include <scene.h>
#include <integrator.h>
#include <util/threadpool.h>
#include <util/field.h>
#include <util/wedge.h>
#include <util/cas.h>
#include <util/moexp.h>



struct GlobalRenderInfo
{
	GlobalRenderInfo():
		debug_pixel(-1, -1)
	{

	}

	// always used ---
	const RenderScene* scene;

	// used for image rendering ---
	Bitmap* image;
	Bitmap* image_transmittance;
	P2i debug_pixel;
	Box2i crop_window;

	// used for fluence rendering ---
	field::VoxelGridField<double>::Ptr fluence_grid;
	Transformd localToWorld;
};



void render_volume( MonteCarloTaskInfo& ti, const GlobalRenderInfo* gi )
{
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
			Color3f f(0.0f);
			Color3f T(1.0f);

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
			//ti.g.scene->camera->sampleRay( P2d(x+ti.rng.next1D(), y+ti.rng.next1D()), rayWS );

			// do raycast ---
			try
			{
				RadianceQuery rq;
				rq.ray = rayWS;
				rq.pixel = V2i(x, y);
				rq.debug = debug;
				f = gi->scene->integrator->Li(gi->scene, rq, ti.rng);

				T = rq.transmittance;
			}
			catch (std::exception& e)
			{
				std::cout << "render_volume: caught exception at task=" << ti.taskid << " x=" << x << " y=" << y << " " << " index=" <<  index << " sample=" << ti.samples << std::endl;
				std::cout << e.what() << '\n';
				std::flush(std::cout);
				throw e;
			}


			if( std::isnan(f.getLuminance()) )
				std::cout << "PathTracingTask::run: got NaN value @ index=" << index << " sample=" << ti.samples << std::endl;

			// update pixel color
			Color3f& c = gi->image->coeffRef(index);
			c += (f - c)/float(ti.samples+1);

			// update transmittance
			Color3f& c_transmittance = gi->image_transmittance->coeffRef(index);
			c_transmittance += (T - c_transmittance)/float(ti.samples+1);

			if(debug)
			{
				c = Color3f(1.0f, 0.0f, 0.0f);
				c_transmittance = Color3f(1.0f, 0.0f, 0.0f);
			}
		} // pixel x
	} // pixel y


	++ti.samples;
}

// globals
VoxelGridf::Ptr g_fluencegrid;
Transformd      g_localToWorld;
RenderScene*          g_scene;

void render_fluence_slice( MonteCarloTaskInfo& ti )
{
	V3i resolution = g_fluencegrid->getResolution();
	int numVoxels = resolution.x()*resolution.y()*resolution.z();


	// note that y goes from bottom=0 to top=max

	for (int voxel = ti.taskid; voxel< numVoxels; voxel += ti.numTasks)
	{
		int i, j, k;
		g_fluencegrid->getCoordFromIndex(voxel, i, j, k);


		// work out voxel center
		P3d pWS = g_localToWorld*g_fluencegrid->voxelToLocal(P3d(i+0.5, j+0.5, k+0.5));
		V3d d = sampleSphere<double>(ti.rng);
		double d_pdf = sampleSpherePDF();

		Ray3d rayWS( pWS, d );

		Color3f f(0.0f);

		// do raycast ---
		try
		{
			RadianceQuery rq;
			rq.ray = rayWS;
			//rq.debug = debug;
			f = g_scene->integrator->Li(g_scene, rq, ti.rng)/d_pdf;
		}
		catch (std::exception& e)
		{
			std::cout << "render_fluence_slice: caught exception at task=" << ti.taskid << " voxel=" << voxel << " sample=" << ti.samples << std::endl;
			std::cout << e.what() << '\n';
			std::flush(std::cout);
			throw e;
		}


		if( std::isnan(f.getLuminance()) )
			std::cout << "PathTracingTask::run: got NaN value @ voxel=" << voxel << " sample=" << ti.samples << std::endl;

		// update voxel value
		float& c = g_fluencegrid->lvalue(i, j, k);
		c += (f.r() - c)/float(ti.samples+1);
	} // voxel


	++ti.samples;
}




namespace houio
{
	namespace util
	{
		// http://www.sidefx.com/docs/houdini/ref/cameralenses
		double fovy( int res_x, int res_y, double aperture_x, double pixel_aspect, double focal_length )
		{
			double aperture_y = (res_y*aperture_x) / (res_x*pixel_aspect);
			double fovy = 2.0*std::atan( (aperture_y/2.0) / focal_length );
			return radToDeg(fovy);
		}
	}

	// NB: houdini uses a y-up, right-handed coordinate system
	struct HScene
	{
		houio::json::ObjectPtr getCamera(const std::string& name)
		{
			return m_cameras[name];
		}


		std::map<std::string, houio::json::ObjectPtr> m_cameras;
	};


	HScene loadScene( const std::string & filename )
	{
		houio::json::JSONReader reader;
		houio::json::Parser p;

		std::ifstream in( filename );

		if(!p.parse( &in, &reader ))
			throw std::runtime_error("houio::loadScene parse failed");

		houio::json::ObjectPtr root = reader.getRoot().asObject();

		if( !root )
			throw std::runtime_error("houio::loadScene parse failed (no root)");

		// ======================================================
		HScene scene;


		// get cameras ----------------------
		houio::json::ObjectPtr cameras = root->getObject("cameras");
		if(cameras)
		{
			std::vector<std::string> cam_names;
			cameras->getKeys(cam_names);
			for( auto cam_name:cam_names )
				scene.m_cameras[cam_name] = cameras->getObject(cam_name);
		}

		return scene;
	}
}

// right-handed: rotation is CW
// left-handed: rotation is CCW

// NB: we use a y-up, left-handed coordinate system
struct Scene
{
	typedef std::shared_ptr<Scene> Ptr;

	std::map<std::string, Camera::Ptr> m_cameras;
};


Scene::Ptr loadScene( const std::string& filename )
{
	houio::HScene hscene = houio::loadScene( filename );

	Scene::Ptr scene = std::make_shared<Scene>();

	// cameras ---------
	for( auto it:hscene.m_cameras )
	{
		std::string cam_name = it.first;
		houio::json::ObjectPtr cam = it.second;


		V2i res( cam->get<int>("resx"), cam->get<int>("resy") );

		double hfov_deg = houio::util::fovy( res.x(), res.y(),
											 cam->get<float>("camera.horizontalFilmAperture", 41.4214f),
											 1.0,
											 cam->get<float>("camera.fl", 50.0f));

		// extract camera parameters...
		double nearClip = 1.0e-4;
		double farClip = 1000.0;
		V3d translation( cam->get<float>("transform.tx"),
						 cam->get<float>("transform.ty"),
						 cam->get<float>("transform.tz")
						 );
		V3d rotatin_euler_angles( cam->get<float>("transform.rx"),
								  cam->get<float>("transform.ry"),
								  cam->get<float>("transform.rz")
								  );

		Eigen::Affine3d cameraToWorld;

		cameraToWorld = Eigen::Translation3d(V3d(translation.x(), translation.y(), translation.z()))*
						Eigen::AngleAxis<double>(degToRad(rotatin_euler_angles.z()), TVector<double, 3>::UnitZ())*
						Eigen::AngleAxis<double>(degToRad(rotatin_euler_angles.y()), TVector<double, 3>::UnitY())*
						Eigen::AngleAxis<double>(degToRad(rotatin_euler_angles.x()), TVector<double, 3>::UnitX());

		//Camera::Ptr camera = std::make_shared<OrthographicCamera>(res.x(), res.y(), 2.0, 2.0, 1e-4f, 5000.0);
		PerspectiveCamera::Ptr camera = std::make_shared<PerspectiveCamera>(res.x(), res.y(), hfov_deg, nearClip, farClip );
		//camera->setCameraToWorld( Eigen::Affine3d(lookAt<double>(targetWS-viewDir*b, targetWS)) );
		//Camera::Ptr camera = createView( V3d(0.0, 0.0, 1.0), scene.bound, res );
		camera->setCameraToWorld(Transformd(cameraToWorld));

		scene->m_cameras[cam_name] = camera;

		std::cout << camera->getResolutionX() << std::endl;
		std::cout << camera->getResolutionY() << std::endl;

		// make some tests
		/*
		{
			RNGd rng(123);
			int numPoints = 30;

			//houio::Geometry::Ptr geo_rays = houio::Geometry::createLineGeometry();
			std::vector<houio::math::V3f> points_0_world;
			std::vector<houio::math::V3f> points_1_perspective;
			std::vector<houio::math::V3f> points_2_perspective_shift;
			std::vector<houio::math::V3f> points_3_perspective_shift_scaled;
			std::vector<houio::math::V3f> points_rotated;

			std::vector<houio::math::V3f> points_0_clip;
			std::vector<houio::math::V3f> points_1_camera;
			std::vector<houio::math::V3f> points_2_raster;

			for( int i=0;i<numPoints;++i )
			{
				double aspect = 1.0;
				// create point in clip space
				P3d p_0(rng.next1D()*2.0-1.0, rng.next1D()*2.0-1.0, rng.next1D()*2.0-1.0);

				P3d p_1 =  (camera->m_clipToCamera * p_0.colwise().homogeneous()).colwise().hnormalized();
				P3d p_2 = camera->m_cameraToRaster*p_1;


				points_0_clip.push_back(houio::math::V3f(float(p_0.x()), float(p_0.y()), float(p_0.z())));
				points_1_camera.push_back(houio::math::V3f(float(p_1.x()), float(p_1.y()), float(p_1.z())));
				points_2_raster.push_back(houio::math::V3f(float(p_2.x()), float(p_2.y()), float(p_2.z())));

			}

			houio::HouGeoIO::xport("cam_0.bgeo", points_0_clip);
			houio::HouGeoIO::xport("cam_1.bgeo", points_1_camera);
			houio::HouGeoIO::xport("cam_2.bgeo", points_2_raster);

			std::vector<houio::Geometry::Ptr> lines;
			Ray3d ray;
			camera->sampleRay( P2d(1.0, 1.0), ray, true );
			lines.push_back( houio::Geometry::createLine( houio::math::V3f(ray.o.x(), ray.o.y(), ray.o.z()), houio::math::V3f(ray.o.x()+ray.d.x(), ray.o.y()+ray.d.y(), ray.o.z()+ray.d.z()) ) );

			houio::HouGeoIO::xport( "cam_rays.bgeo", houio::Geometry::merge(lines) );
		}
		*/
	}


	return scene;
}



int main()
{


	// lets test loading houdini scene information
	Scene::Ptr scene = loadScene("c:/projects/epfl/data/scarf.scn");
	//return 0;



	//std::string basePath = "c:/projects/visus/data";
	std::string basePath = ".";
	//std::string outpath = basePath + "/noisereduction/fitting4";
	std::string outpath = basePath + "";




	// output image ---
	V2i res = V2i(512, 512);
	//V2i res = V2i(2048, 2048);
	Bitmap image_color(res);
	Bitmap image_transmittance(res);


	// scene elements ---
	//Volume::Ptr volume = volumes::C60();
	//Volume::Ptr volume = volumes::nebulae();
	//Volume::Ptr volume = volumes::homogeneous(V3d(5.0), V3d(0.8));
	Volume::Ptr volume = volumes::scarf();
	//return 0;
	//Light::Ptr light = lights::point( volume->getBound().getCenter() + P3d(0.0, volume->getBound().getExtents().y()*0.6, 0.0) );
	Light::Ptr light = lights::directional( V3d(0.0, -1.0, 0.0), volume->getBound() );

	RenderScene renderscene;
	renderscene.bound = volume->getBound();
	renderscene.id = "nebulae";
	renderscene.volume = volume.get();
	renderscene.light = light.get();
	renderscene.camera = scene->m_cameras["cam1"].get();


	// RENDERING -----------------------------------
	//for( int i=0;i<=30;++i )
	{
		image_color.fill(Color3f(0.0f, 0.0f, 0.0f));
		image_transmittance.fill(Color3f(0.0f, 0.0f, 0.0f));
		//shcache->m_evaluationOrder = i;

		//Integrator::Ptr integrator = integrators::raymarcher(0.005);
		Integrator::Ptr integrator = integrators::volume_path_tracer(100);
		//Integrator::Ptr integrator = integrators::volume_path_tracer_cached(shcache);
		//Integrator::Ptr integrator = integrators::cache_raymarcher(0.05, fluencecache_gpu);
		//Bitmap::Ptr pixel_estimates = readImage(outpath + "/nebulae_pixel_estimate.exr");
		//Integrator::Ptr integrator = integrators::adrrs_volume_path_tracer(pixel_estimates, fluence_field);
		renderscene.integrator = integrator.get();

		int numSamples = 10;

		// prepare render info ---
		GlobalRenderInfo gi;
		gi.scene = &renderscene;
		gi.image = &image_color;
		gi.image_transmittance = &image_transmittance;
		gi.crop_window = Box2i( V2i(0,0), res );
		//RenderTaskInfo::g.debug_pixel = P2i(318, 209);
		//gi.debug_pixel = P2i(256, 256);
		//gi.debug_pixel = P2i(128, 109);
		//RenderTaskInfo::g.debug_pixel = P2i(340, 340);
		//RenderTaskInfo::g.debug_pixel = P2i(195, 512-25);


		// execute multithreaded render ---
		Terminator terminator(numSamples);
		std::cout << "rendering image..."<< std::endl;std::flush(std::cout);
		runGenericTasks<MonteCarloTaskInfo, GlobalRenderInfo>( render_volume, &gi, terminator, ThreadPool::getNumSystemCores() );


		// save results ---
		image_color.saveEXR( outpath + "/" + renderscene.id + "_" + integrator->getId() + "_color.exr" );
	}


	// compute fluence grid ----------------
	/*
	{
		int numSamples = 100;
		V3i resolution(64, 64, 64);

		Integrator::Ptr integrator = integrators::volume_path_tracer(100);
		scene.integrator = integrator.get();

		g_fluencegrid = VoxelGrid<float>::create(resolution);
		g_localToWorld = volume->getLocalToWorld();
		g_scene = &scene;

		// execute multithreaded render ---
		Terminator terminator(numSamples);
		std::cout << "rendering fluence grid..."<< std::endl;std::flush(std::cout);
		runGenericTasks<MonteCarloTaskInfo>( render_fluence_slice, terminator, ThreadPool::getNumSystemCores() );

		//
		Fieldf::Ptr fluencefield = std::make_shared<field::VoxelGridFieldf>(resolution, g_fluencegrid->getRawPointer());
		field::write( "nebulae_fluence_cpu.bgeo", field::xform<float>(fluencefield, g_localToWorld), resolution, g_localToWorld );
	}
	*/




	return 0;
}


