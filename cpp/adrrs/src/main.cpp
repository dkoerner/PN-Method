#include <iostream>
#include <complex>

#include <scene.h>
#include <integrator.h>
#include <util/threadpool.h>
#include <util/field.h>
#include <util/wedge.h>
#include <util/cas.h>
#include <util/moexp.h>



#include <shcache.h>





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

// globals
VoxelGridf::Ptr g_fluencegrid;
Transformd      g_localToWorld;
Scene*          g_scene;

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


int main()
{
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
	Volume::Ptr volume = volumes::nebulae();
	//Light::Ptr light = lights::point( volume->getBound().getCenter() + P3d(0.0, volume->getBound().getExtents().y()*0.6, 0.0) );
	Light::Ptr light = lights::directional( V3d(0.0, -1.0, 0.0), volume->getBound() );

	Scene scene;
	scene.bound = volume->getBound();
	scene.id = "nebulae";
	scene.volume = volume.get();
	scene.light = light.get();


	//Camera::Ptr camera = createView( V3d(0.0, 0.0, 1.0), scene.bound, res );



	// load sh cache -----------------------
	SHCache::Ptr shcache = std::make_shared<SHCache>("nebulae200.shcache");

	// verification
	/*
	for( int l=0;l<=shcache.m_order;++l )
	{
		std::string filename = "sh_reconstruction_alt_$0.bgeo";
		filename = replace(filename, "$0", toString(l));

		Color3f* coeffs = shcache.getCoefficients(0);
		rasterizeSphericalFunctionSphere(filename, [&](double theta, double phi)->Color3f
		{
			return moexp::Y_real_sum<Color3f>(l, coeffs, theta, phi);
		}, 8.0);
	}
	*/


	// pn analysis and debugging ---
	/*
	{
		Integrator::Ptr integrator = integrators::dummy();
		scene.integrator = integrator.get();

		PNCache cache;
		cache.m_override = true;
		cache.m_doZeroScattering = false;
		int numMoments = 4;
		int numSamples = 1000000;
		int res = 1;
		cache.generate( outpath + "/analysis", &scene, numMoments, numSamples, res );


		Wedge wedge;
		std::vector<int> moment_values = {0, 1, 2, 3, 4};
		wedge.addParm("moment", moment_values);
		wedge.build();

		std::cout << "4pi=" << 4.0*M_PI << std::endl;
		std::cout << "4pi/3=" << 4.0*M_PI/3.0 << std::endl;

		std::vector<Wedge::Iteration> iterations = wedge.iterations();
		for( auto it : iterations )
		{
			int moment_index = it.getInt("moment");
			std::cout << "moment=" << moment_index << std::endl;

			tensor::Tensor<double> moment_tensor = cache.getMoment( 0, 0, 0, moment_index );

			for(auto component = moment_tensor.begin(), end=moment_tensor.end(); component!=end;++component)
			{
				std::cout << "\t" << component.index_str() << "=" << component.value() << std::endl;

//				houio::Geometry::Ptr geo = houio::Geometry::createSphere(50, 50, 1.0);
//				houio::Attribute::Ptr pattr = geo->getAttr("P");
//				int numPoints = pattr->numElements();
//				for(int i=0;i<numPoints;++i)
//				{
//					houio::math::V3f d = pattr->get<houio::math::V3f>(i).normalized();
//					d *= std::abs(component.weight(V3d(d.x, d.y, d.z)));
//					pattr->set<houio::math::V3f>(i, d);
//				}
//				houio::HouGeoIO::xport(it.expand_value( outpath+"/analysis_$0_" + component.index_str() + ".bgeo"), geo);
			}
		}
	}
	*/


	Eigen::Affine3d cameraToWorld;
	cameraToWorld = Eigen::Translation3d(V3d(0.0, 0.0, -2.5));

	//Camera::Ptr camera = std::make_shared<OrthographicCamera>(res.x(), res.y(), 2.0, 2.0, 1e-4f, 5000.0);
	Camera::Ptr camera = std::make_shared<PerspectiveCamera>(res.x(), res.y() );
	//camera->setCameraToWorld( Eigen::Affine3d(lookAt<double>(targetWS-viewDir*b, targetWS)) );
	camera->setCameraToWorld(Transformd(cameraToWorld));

	scene.camera = camera.get();


	// RENDERING -----------------------------------
	for( int i=0;i<=30;++i )
	{
		image_color.fill(Color3f(0.0f, 0.0f, 0.0f));
		image_transmittance.fill(Color3f(0.0f, 0.0f, 0.0f));
		shcache->m_evaluationOrder = i;

		//Integrator::Ptr integrator = integrators::raymarcher(0.005);
		//Integrator::Ptr integrator = integrators::volume_path_tracer(100, false);
		Integrator::Ptr integrator = integrators::volume_path_tracer_cached(shcache);
		//Integrator::Ptr integrator = integrators::cache_raymarcher(0.05, fluencecache_gpu);
		//Bitmap::Ptr pixel_estimates = readImage(outpath + "/nebulae_pixel_estimate.exr");
		//Integrator::Ptr integrator = integrators::adrrs_volume_path_tracer(pixel_estimates, fluence_field);
		scene.integrator = integrator.get();

		int numSamples = 10;

		// prepare render info ---
		RenderTaskInfo::g.scene = &scene;
		RenderTaskInfo::g.image = &image_color;
		RenderTaskInfo::g.image_transmittance = &image_transmittance;
		RenderTaskInfo::g.crop_window = Box2i( V2i(0,0), res );
		//RenderTaskInfo::g.debug_pixel = P2i(318, 209);
		//RenderTaskInfo::g.debug_pixel = P2i(256, 256);
		//RenderTaskInfo::g.debug_pixel = P2i(340, 340);
		//RenderTaskInfo::g.debug_pixel = P2i(195, 512-25);


		// execute multithreaded render ---
		Terminator terminator(numSamples);
		std::cout << "rendering image..."<< std::endl;std::flush(std::cout);
		runGenericTasks<RenderTaskInfo>( render_volume, terminator, ThreadPool::getNumSystemCores() );


		// save results ---
		flip(image_color);
		flip(image_transmittance);
		image_color.saveEXR( outpath + "/" + scene.id + "_" + integrator->getId() + "_color.exr" );
		std::string transmittance_exr = outpath + "/" + scene.id + "_" + integrator->getId() + "_transmittance.exr";
		image_transmittance.saveEXR(transmittance_exr);
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


