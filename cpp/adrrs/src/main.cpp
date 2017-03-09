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



	// generate sh cache -----------------------
	{
		Integrator::Ptr integrator = integrators::volume_path_tracer(10);
		scene.integrator = integrator.get();
		std::string filename = outpath + "/nebulae.shcache";


		RNGd rng;


		SHCache cache;
		cache.m_override = true;
		int order = 30;
		int numSamples = 50000;
		int res = 1;
		//cache.generate( filename, &scene, order, numSamples, res );
		cache.generate_images( filename, &scene, order, numSamples, res );
		/*
		P3d p = cache.m_localToWorld*P3d(0.5, 0.5, 0.5);

		for( int i=0;i<order;++i )
		{
			cache.m_evaluationOrder = i;

			moexp::SphericalFunction<Color3f> fun = [&](double theta, double phi)->Color3f
			{
				V3d d = sphericalDirection(theta, phi);
				Color3f c = cache.eval(p, d);
				return c;
			};

			std::string filename("shcache_reconstruction_$0.exr");
			filename = replace(filename, "$0", toString(i));
			rasterizeSphericalFunctionMap( filename, fun );
		}
		*/

		/*
		std::string groundtruth_map_filename = "groundtruth_spherical_function.exr";
		//std::string groundtruth_map_filename = "envmap.exr";

		double exposure = 8.0;

		std::vector<EnvMap> filtered_maps;
		auto filter_widths = linearSamples(0.0f, 10.0f, 5 );
		//std::vector<float> filter_widths = {0.0};
		for( int i=0;i<filter_widths.size();++i )
		{
			float filter_width = filter_widths[i];
			filtered_maps.push_back( EnvMap(groundtruth_map_filename) );
			EnvMap& map = filtered_maps[i];
			if(filter_width>0.0)
				map.bitmap() = filter(map.bitmap(), gaussFilter(filter_width));

			std::string filename = "groundtruth_$0.exr";
			filename = replace(filename, "$0", toString(i));
			std::string filename_geo = "groundtruth_$0.bgeo";
			filename_geo = replace(filename_geo, "$0", toString(i));

			map.bitmap().saveEXR(filename);
			map.saveGeo(filename_geo, exposure);
		}



		Wedge wedge;
		//std::vector<int> order_list = {30};
		std::vector<int> order_list = linearSamples(0, 30, 30);
		wedge.addParm( "order", order_list );
		wedge.addParm("filterwidth", filter_widths);
		for( auto it = wedge.begin(), end = wedge.end();it!=end;++it )
		{
			int order = it.getInt("order");
			int filter_width_index = it.getIndex("filterwidth");


			EnvMap& map = filtered_maps[filter_width_index];

			moexp::SphericalFunction<Color3f> fun = [&](double theta, double phi)->Color3f
			{
				return map.eval(theta, phi);
			};
			std::unique_ptr<std::vector<Color3f>> coeffs = moexp::project_Y_real( order, fun, numSamples );
			moexp::SphericalFunction<Color3f> fun_reconstruction = [&](double theta, double phi)->Color3f
			{
				return moexp::Y_real_sum(order,*coeffs, theta, phi);
			};

			rasterizeSphericalFunctionMap( it.expand_index( "sh_reconstruction_$0_$1.exr"), fun_reconstruction );
			rasterizeSphericalFunctionSphere( it.expand_index( "sh_reconstruction_$0_$1.bgeo"), fun_reconstruction, exposure );
		}
		*/
		return 0;


		/*
		rasterizeSphericalFunctionSphere("groundtruth.bgeo", [&](double theta, double phi)->Color3f
		{
			return map.eval(theta, phi)*std::pow(2.0, exposure);
		});
		*/

		//int order = 15;
		//for( int l=0;l<=order;++l )
		{
			/*
			std::unique_ptr<std::vector<Color3f>> coeffs = moexp::project_Y_real<Color3f>( order, fun, 100 );

			moexp::SphericalFunction<Color3f> fun_reconstruction = [&](double theta, double phi)->Color3f
			{
				return moexp::Y_real_sum<Color3f>( l, coeffs->data(), theta, phi );
			};
			*/



//			std::string filename("shcache_reconstruction_$0.exr");
//			filename = replace(filename, "$0", toString(l));
//			//rasterizeSphericalFunctionMap( filename, fun_reconstruction );

//			///*
//			EnvMap map(filename);
//			std::string filename_geo("shcache_reconstruction_$0.bgeo");
//			filename_geo = replace(filename_geo, "$0", toString(l));

//			rasterizeSphericalFunctionSphere(filename_geo, [&](double theta, double phi)->Color3f
//			{
//				return map.eval(theta, phi)*std::pow(2.0, exposure);
//			});
			//*/
		}




		/*
		moexp::SphericalFunction<Color3f> fun_groundtruth = [&](double theta, double phi)->Color3f
		{
			Color3f result(0.0);
			for( int i=0;i<numSamples;++i )
			{
				V3d d = sphericalDirection(theta, phi);
				RadianceQuery rq;
				rq.ray = Ray3d(p, d);
				Color3f sample = scene.integrator->Li(&scene, rq, rng);
				result+=(sample-result)/double(i+1);
			}
			return result;
		};
		rasterizeSphericalFunctionMap( "groundtruth.exr", fun_groundtruth );
		*/
	}

	return 0;
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


	/*
	// RENDERING -----------------------------------
	{
		//Integrator::Ptr integrator = integrators::raymarcher(0.005);
		PNCache::Ptr cache = std::make_shared<PNCache>(outpath + "/nebulae.moments");
		cache->eval(P3d(0.0f, 0.0, 0.0), normalize(V3d(1.0, 1.0, 1.0)), true);

		//Bitmap::Ptr pixel_estimates = readImage(outpath + "/nebulae_pixel_estimate.exr");
		Integrator::Ptr integrator = integrators::cache_raymarcher(0.005, cache);
		//Integrator::Ptr integrator = integrators::adrrs_volume_path_tracer(pixel_estimates, fluence_field);
		scene.integrator = integrator.get();

		int numSamples = 1;

		// prepare render info ---
		RenderTaskInfo::g.scene = &scene;
		RenderTaskInfo::g.image = &image_color;
		RenderTaskInfo::g.image_transmittance = &image_transmittance;
		RenderTaskInfo::g.crop_window = Box2i( V2i(0,0), res );
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
		image_color.saveEXR( outpath + "/" + scene.id + "_color.exr" );
		std::string transmittance_exr = outpath + "/test_transmittance.exr";
		image_transmittance.saveEXR(transmittance_exr);
	}
	*/




	return 0;
}


