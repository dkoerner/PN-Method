#include <iostream>
#include <sys/stat.h>

#include <scene.h>
#include <integrator.h>
#include <util/threadpool.h>
#include <util/field.h>
#include <util/dda3d.h>

#include <houio/Geometry.h>


inline bool file_exists(const std::string& name)
{
	struct stat buffer;
	return (stat (name.c_str(), &buffer) == 0);
}



void render_volume( RenderTaskInfo& ti )
{
	int width = ti.g.scene->camera->getResolutionX();
	int height = ti.g.scene->camera->getResolutionY();


	// note that y goes from bottom=0 to top=max

	for (int scanline = ti.taskid; scanline < ti.g.crop_window.max.y(); scanline += ti.numTasks)
	{
		if( scanline < ti.g.crop_window.min.y() )
			continue;

		for (int x = ti.g.crop_window.min.x(); x < ti.g.crop_window.max.x(); ++x)
		{
			int y = height-1-scanline;
			int index = y*width + x;
			Color3f f(0.0f);
			Color3f T(1.0f);


			bool debug = false;

			if( ti.g.debug_pixel.x()>0 )
			{
				debug = (x == ti.g.debug_pixel.x()) &&
						(y == ti.g.debug_pixel.y());

				if(!debug)
					continue;
			}

			Ray3d rayWS;
			ti.g.scene->camera->sampleRay( P2d(x+0.5, y+0.5), rayWS, debug );
			//ti.g.scene->camera->sampleRay( P2d(x+ti.rng.next1D(), y+ti.rng.next1D()), rayWS );

			// do raycast ---
			try
			{
				RadianceQuery rq;
				rq.ray = rayWS;
				rq.pixel = V2i(x, y);
				rq.debug = debug;
				f = ti.g.scene->integrator->Li(ti.g.scene, rq, ti.rng);
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
			Color3f& c = ti.g.image->coeffRef(index);
			c += (f - c)/float(ti.samples+1);

			// update transmittance
			Color3f& c_transmittance = ti.g.image_transmittance->coeffRef(index);
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

void render_fluence( RenderTaskInfo& ti )
{
	V3i resolution = ti.g.fluence_grid->grid.getResolution();
	int numVoxels = resolution.x()*resolution.y()*resolution.z();
	for (int voxel = ti.taskid; voxel< numVoxels; voxel += ti.numTasks)
	{
		// work out voxel coordinate
		int i,j, k;
		ti.g.fluence_grid->grid.getCoordFromIndex(voxel, i, j, k);

		// sample position within voxel
		P3d pVS(i+ti.rng.next1D(), j+ti.rng.next1D(), k+ti.rng.next1D());
		P3d pLS = ti.g.fluence_grid->grid.voxelToLocal(pVS);
		P3d pWS = ti.g.localToWorld*pLS;

		// sample direction
		V3d d = sampleSphere<double>(ti.rng);
		double d_pdf = sampleSpherePDF();

		// integrate
		RadianceQuery rq;
		rq.ray = Ray3d(pWS, d);
		//rq.debug = true;
		double fluence_sample = ti.g.scene->integrator->Li( ti.g.scene, rq, ti.rng )[0]*INV_FOURPI;

		// accumulate
		double& fluence = ti.g.fluence_grid->grid.lvalue(i, j, k);
		fluence += (fluence_sample/d_pdf - fluence)/float(ti.samples+1);
	}

	++ti.samples;
}

int main()
{
	//std::string basePath = "c:/projects/visus/data";
	std::string basePath = ".";
	//std::string outpath = basePath + "/noisereduction/fitting4";
	std::string outpath = basePath + "";



	int numSamples = 10000;

	// output image ---
	V2i res = V2i(512, 512);
	//V2i res = V2i(2048, 2048);
	Bitmap image_color(res);
	Bitmap image_transmittance(res);


	// scene elements ---
	//Volume::Ptr volume = volumes::C60();
	Volume::Ptr volume = volumes::nebulae();
	//Light::Ptr light = lights::point( volume->getBound().getCenter() + P3d(0.0, volume->getBound().getExtents().y()*0.6, 0.0) );
	Light::Ptr light = lights::directional( V3d(0.0, -1.0, 0.0), volume->getBound() );
	//Integrator::Ptr integrator = integrators::raymarcher(0.005);
	Integrator::Ptr integrator = integrators::volume_path_tracer();

	Scene scene;
	scene.bound = volume->getBound();
	scene.id = "nebulae";
	scene.volume = volume.get();
	scene.light = light.get();
	scene.integrator = integrator.get();

	//Camera::Ptr camera = createView( V3d(0.0, 0.0, 1.0), scene.bound, res );


	Eigen::Affine3d cameraToWorld;
	cameraToWorld = Eigen::Translation3d(V3d(0.0, 0.0, -2.5));

	//Camera::Ptr camera = std::make_shared<OrthographicCamera>(res.x(), res.y(), 2.0, 2.0, 1e-4f, 5000.0);
	Camera::Ptr camera = std::make_shared<PerspectiveCamera>(res.x(), res.y() );
	//camera->setCameraToWorld( Eigen::Affine3d(lookAt<double>(targetWS-viewDir*b, targetWS)) );
	camera->setCameraToWorld(Transformd(cameraToWorld));

	scene.camera = camera.get();

	// compute fluence grid -----------------------
	{
		int numSamples = 500;
		int res = 128;
		V3i resolution(res, res, res);
		Transformd localToWorld = scene.volume->getLocalToWorld();

		// rasterize local space
		// parallel version ---
		RenderTaskInfo::g.scene = &scene;
		RenderTaskInfo::g.fluence_grid = std::make_shared<field::VoxelGridField<double>>(resolution);
		RenderTaskInfo::g.localToWorld = localToWorld;

		Terminator terminator(numSamples);
		std::cout << "rendering fluence..."<< std::endl;std::flush(std::cout);
		runGenericTasks<RenderTaskInfo>( render_fluence, terminator, ThreadPool::getNumSystemCores() );

		field::write(outpath + "/nebulae_fluence.bgeo", field::xform<double>(RenderTaskInfo::g.fluence_grid, RenderTaskInfo::g.localToWorld), V3i(res, res, res), RenderTaskInfo::g.localToWorld );
	}

	/*
	// RENDERING -----------------------------------
	{
		// prepare render info ---
		RenderTaskInfo::g.scene = &scene;
		RenderTaskInfo::g.image = &image_color;
		RenderTaskInfo::g.image_transmittance = &image_transmittance;
		RenderTaskInfo::g.crop_window = Box2i( V2i(0,0), res );
		//RenderTaskInfo::g.sample_function = integrator.sample_function;
		//RenderTaskInfo::g.debug_pixel = P2i(318, 209);
		//RenderTaskInfo::g.debug_pixel = P2i(256, 256);
		//RenderTaskInfo::g.debug_pixel = P2i(340, 340);


		// execute multithreaded render ---
		Terminator terminator(numSamples);
		std::cout << "rendering image..."<< std::endl;std::flush(std::cout);
		runGenericTasks<RenderTaskInfo>( render_volume, terminator, ThreadPool::getNumSystemCores() );


		// save results ---
		flip(image_color);
		flip(image_transmittance);
		image.saveEXR( outpath + "/" + scene.id + "_color.exr" );
		image_color.saveEXR(color_exr);
		std::string transmittance_exr = outpath + "/test_transmittance.exr";
		image_transmittance.saveEXR(transmittance_exr);
	}
	*/




	return 0;
}

