#include <iostream>


#include <scene.h>
#include <integrator.h>
#include <util/threadpool.h>
#include <util/field.h>
#include <util/sh.h>

#include <houio/Geometry.h>

#include <pncache.h>





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




struct EnvMap
{
	typedef std::shared_ptr<EnvMap> Ptr;

	EnvMap( const std::string& filename )
	{
		m_bitmap = Bitmap(filename);
		m_transform = Transformd();
	}

	// evaluate environment map
	Color3f eval( double theta, double phi )const
	{
		V3d d = sphericalDirection<double>(theta, phi);
		P2d uv = directionToUV(d);
		return m_bitmap.eval(uv);
	}


	P2d directionToUV( const V3d& d )const
	{
		// using formulas given in http://gl.ict.usc.edu/Data/HighResProbes/
		// with the difference that u=[0,1] (instead of [0,2]) and we negate z
		P2d uv( (1+std::atan2(d.x(), d.z())/M_PI)/2,
				 safe_acos(d.y())/M_PI );
		return uv;
	}
	V3d uvToDirection( const P2d& uv )const
	{
		// using formulas given in http://gl.ict.usc.edu/Data/HighResProbes/
		// with the difference that u=[0,1] (instead of [0,2]) and we negate z
		// azimuthal angle
		double theta = M_PI*(uv.x()*2.0-1.0);
		// elevation angle
		double phi = M_PI*uv.y();
		return V3d( std::sin(phi)*std::sin(theta), std::cos(phi), std::sin(phi)*cos(theta) );
	}
	P2d uvToXY(const P2d& uv)const
	{
		P2d xy(
			(uv.x()*(m_bitmap.cols()-1)),
			(uv.y()*(m_bitmap.rows()-1))
			);
		return xy;
	}

	P2d xyToUV(const P2d& xy)const
	{
		return P2d(
			(xy.x())/double(m_bitmap.cols()-1),
			(xy.y())/double(m_bitmap.rows()-1)
			);
	}

	V3d xyToDirection( const P2d& xy )const
	{
		return uvToDirection( xyToUV(xy) );
	}
	P2d directionToXY( const V3d& d )const
	{
		return uvToXY(directionToUV(d));
	}

private:
	Transformd m_transform;
	Bitmap m_bitmap;
};


void writeSphericalFunction(const std::string& filename, sh::SphericalFunction<Color3f> func )
{
	houio::Geometry::Ptr geo = houio::Geometry::createSphere(120, 120, 1.0);
	houio::Attribute::Ptr pAttr = geo->getAttr("P");
	houio::Attribute::Ptr cdAttr = houio::Attribute::createV3f(pAttr->numElements());
	for( int i=0;i<pAttr->numElements();++i )
	{
		houio::math::V3f p = pAttr->get<houio::math::V3f>(i);
		P2d theta_phi = sphericalCoordinates<double>(V3d(p.x, p.y, p.z));
		double theta = theta_phi.x();
		double phi = theta_phi.y();
		Color3f col = func(theta, phi);
		cdAttr->set<houio::math::V3f>( i, houio::math::V3f(col.r(), col.g(), col.b()) );
	}
	geo->setAttr("Cd", cdAttr);
	houio::HouGeoIO::xport( filename, geo);
}


double test( int n, int l )
{
	int n2 = 2*n;
	int l2 = 2*l;
	return tensor::ipow(-1, n)*double(sh::factorial(l)*sh::doubleFactorial(l2-n2-1))/double(sh::factorial(l-n2)*sh::doubleFactorial(l2-1)*sh::doubleFactorial(n2));
}

int main()
{

	std::cout << "kuckuck\n";

	int l = 3;
	int n_end = std::ceil( double(l)/2.0 );

	for( int n=0;n<n_end;++n )
	{
		std::cout << test(n, l)<<std::endl;
	}



	/*
	EnvMap map("envmap.exr");

	sh::SphericalFunction<Color3f> fun = [&](double theta, double phi) -> Color3f
	{
		return map.eval(theta, phi);
	};

	writeSphericalFunction( "groundtruth.bgeo", fun );

	// project
	int order = 30;
	int numSamples = 50000;
	std::unique_ptr<std::vector<Color3f>> sh_coeffs;
	sh_coeffs = sh::project<Color3f>(order, fun, numSamples);

	// reconstruction
	for( int l=0;l<=order;++l )
	{
		sh::SphericalFunction<Color3f> reconstruction = [&](double theta, double phi) -> Color3f
		{
			return sh::evalSum(l, *sh_coeffs.get(), theta, phi);
		};

		std::string filename("testfit_reconstruction_$0.bgeo");
		filename = replace(filename, "$0", toString(l));

		writeSphericalFunction( filename, reconstruction );
	}
	*/


	/*
	int L = 6;
	for( int l=0;l<L;++l )
	{
		for( int m=-l;m<=l;++m )
		{

			// visualization for houdini
			{
				houio::Geometry::Ptr geo = houio::Geometry::createSphere(120, 120, 1.0);
				houio::Attribute::Ptr pAttr = geo->getAttr("P");
				houio::Attribute::Ptr cdAttr = houio::Attribute::createV3f(pAttr->numElements());
				for( int i=0;i<pAttr->numElements();++i )
				{
					houio::math::V3f p = pAttr->get<houio::math::V3f>(i);
					P2d theta_phi = sphericalCoordinates<double>(V3d(p.x, p.y, p.z));
					double theta = theta_phi.x();
					double phi = theta_phi.y();
					//Color3f col = map.eval(theta_phi.x(), theta_phi.y());

					float sh = EvalSHSlow(l, m, theta, phi);

					//Color3f col(sh);
					//cdAttr->set<houio::math::V3f>( i, houio::math::V3f(col.r(), col.g(), col.b()) );

					pAttr->set<houio::math::V3f>( i, p*sh );
				}
				//geo->setAttr("Cd", cdAttr);
				std::string filename("test_$0_$1.bgeo");
				filename = replace(filename, "$0", toString(l));
				filename = replace(filename, "$1", toString(m));
				houio::HouGeoIO::xport( filename, geo);
			}
		}
	}
	*/




	return 0;


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



	/*
	// generate pn cache -----------------------
	{
		Integrator::Ptr integrator = integrators::volume_path_tracer();
		scene.integrator = integrator.get();
		std::string filename = outpath + "/nebulae.moments";

		PNCache cache;
		int numMoments = 4;
		int numSamples = 100;
		int res = 128;
		cache.generate( outpath + "/nebulae.moments", &scene, numMoments, numSamples, res );
		return;
	}
	*/


	Eigen::Affine3d cameraToWorld;
	cameraToWorld = Eigen::Translation3d(V3d(0.0, 0.0, -2.5));

	//Camera::Ptr camera = std::make_shared<OrthographicCamera>(res.x(), res.y(), 2.0, 2.0, 1e-4f, 5000.0);
	Camera::Ptr camera = std::make_shared<PerspectiveCamera>(res.x(), res.y() );
	//camera->setCameraToWorld( Eigen::Affine3d(lookAt<double>(targetWS-viewDir*b, targetWS)) );
	camera->setCameraToWorld(Transformd(cameraToWorld));

	scene.camera = camera.get();

	///*
	// RENDERING -----------------------------------
	{
		//Integrator::Ptr integrator = integrators::raymarcher(0.005);
		PNCache::Ptr cache = std::make_shared<PNCache>(outpath + "/nebulae.moments");
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
	//*/




	return 0;
}

