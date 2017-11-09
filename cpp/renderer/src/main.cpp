


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
#include <integrators/simplept.h>
#include <fields/VoxelGridField.h>

#include <util/threadpool.h>

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
};

struct GlobalComputeFluenceInfo
{
	GlobalComputeFluenceInfo()
	{
	}

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
		V3d direction = sampleSphere(ti.rng);
		double direction_pdf = sampleSpherePDF();

		Ray3d rayWS( pWS, direction );

		// do raycast ---
		try
		{
			RadianceQuery rq;
			rq.ray = rayWS;
			rq.pixel = V2i(-1, -1);
			rq.debug = debug;
			if(gi->scene->volume->getBound().contains(pWS))
				rq.volume = gi->scene->volume;

			// lightsample is used to sample direct light
			LightSample ls;
			ls.refP = pWS;

			f = gi->scene->sample_attenuated_directlight( ls, ti.rng, debug )+ // direct light
				gi->scene->integrator->Li(gi->scene, rq, ti.rng)/direction_pdf; // indirect light
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



Eigen::VectorXd compute_fluence( Volume::Ptr volume, Light::Ptr light, Integrator::Ptr integrator, const Eigen::MatrixXd& points, int seed )
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

	MonteCarloTaskInfo::g_seed = seed;

	int numSamples = 10000;

	Terminator terminator(numSamples);
	std::cout << "computing fluence..."<< std::endl;std::flush(std::cout);
	runGenericTasks<MonteCarloTaskInfo, GlobalComputeFluenceInfo>( compute_fluence_thread, &gi, terminator, ThreadPool::getNumSystemCores() );

	Eigen::VectorXd result( fluence.size() );
	for( int i=0;i<numPoints;++i )
		result.coeffRef(i) = luminance(fluence[i]);
	return result;
}

/*
Image::Ptr render_radiance_distribution_at_point(Volume::Ptr volume, Light::Ptr light, Integrator::Ptr integrator, const P3d& pWS, int seed)
{
	std::cout << "render_radiance_distribution_at_point:\n";
	std::cout << volume->toString();
	std::cout << light->toString();
	std::cout << integrator->toString();

	V2i resolution(512, 256);
	Image::Ptr result = std::make_shared<Image>(resolution);

	Scene scene;
	scene.camera = 0;
	scene.integrator = integrator.get();
	scene.light = light.get();
	scene.volume = volume.get();

	GlobalRenderInfo gi;
	gi.scene = &scene;
	gi.image = result.get();
	gi.crop_window = Box2i( V2i(0,0), camera->getResolution() );
	//gi.debug_pixel = V2i(256, 256);

	int numSamples = 1;

	Terminator terminator(numSamples);
	std::cout << "rendering image..."<< std::endl;std::flush(std::cout);
	runGenericTasks<MonteCarloTaskInfo, GlobalRenderInfo>( render_thread, &gi, terminator, ThreadPool::getNumSystemCores() );
	return result;
}
*/




void render_thread( MonteCarloTaskInfo& ti, const GlobalRenderInfo* gi )
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
			//ti.g.scene->camera->sampleRay( P2d(x+ti.rng.next1D(), y+ti.rng.next1D()), rayWS );

			// do raycast ---
			try
			{
				RadianceQuery rq;
				rq.ray = rayWS;
				rq.pixel = V2i(x, y);
				rq.debug = debug;
				f = gi->scene->integrator->Li(gi->scene, rq, ti.rng);
				//T = rq.transmittance;
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
				c = V3d(1.0, 0.0, 0.0);
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
	//gi.debug_pixel = V2i(1, 0);

	//int numSamples = 1;

	Terminator terminator(numSamples);
	//Terminator terminator(90.0);
	std::cout << "rendering image..."<< std::endl;std::flush(std::cout);
	runGenericTasks<MonteCarloTaskInfo, GlobalRenderInfo>( render_thread, &gi, terminator, ThreadPool::getNumSystemCores() );


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


Light::Ptr create_point_light( const Eigen::Vector3d& position, const Eigen::Vector3d& power )
{
	return std::make_shared<PointLight>(position, power);
}
/*
Light::Ptr create_directional_light()
{
}
*/

Camera::Ptr create_perspective_camera( int xres, int yres, double hfov )
{
	return std::make_shared<PerspectiveCamera>(xres, yres, hfov);
}


Camera::Ptr create_sphere_camera( const Eigen::Vector3d& pWS, int xres, int yres )
{
	return std::make_shared<SphereCamera>( pWS, xres, yres );
}


Integrator::Ptr create_simplept_integrator()
{
	int maxDepth = std::numeric_limits<int>::max();
	bool doSingleScattering = true;
	return std::make_shared<SimplePT>(maxDepth, doSingleScattering);
}

Volume::Ptr create_volume()
{
	return std::make_shared<Volume>();
}

Field3d::Ptr create_constant_field3d( const Eigen::Vector3d& value )
{
	return std::make_shared<ConstantField3d>(value);
}


PNSolution::Ptr load_pnsolution( const std::string& filename )
{
	return std::make_shared<PNSolution>(filename);
}


Image::Ptr load_image( const std::string& filename )
{
	return std::make_shared<Image>(filename);
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
		PNSolution::SHVector::staticInit();

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
		return PNSolution::SHVector::eval( theta, phi, shv.m_coeffs.data(), m_order );
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




PYBIND11_MODULE(renderer, m)
{

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

	// Volume ============================================================
	py::class_<Volume, Volume::Ptr> class_volume(m, "Volume");
	class_volume
		.def( "setBound",
		[](Volume &m, const Eigen::Vector3d& min, const Eigen::Vector3d& max )
		{
			m.setBound(Box3d(min, max));
		})
		.def( "setExtinctionAlbedo", &Volume::setExtinctionAlbedo )
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

	// VoxelGrid ============================================================
	py::class_<VoxelGridFieldcd, VoxelGridFieldcd::Ptr> class_VoxelGridFieldcd(m, "VoxelGridFieldcd", class_fieldcd);
	class_VoxelGridFieldcd
	.def("__init__",
	[](VoxelGridFieldcd &m, py::array b)
	{
		typedef Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic> Strides;

		// location within voxel (center)
		Eigen::Matrix<double, 3, 1> offset(0.5, 0.5, 0.5);

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
		new (&m) VoxelGridFieldcd( data, resolution, offset );
	})
	//.def("test", &VoxelGridField::test)
	//.def("getSlice", &VoxelGridField::getSlice)
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
		const Eigen::VectorXcd& data
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
		return m.sample(pWS, pdf, sample);
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
	// lights
	m.def( "create_point_light", &create_point_light );
	//m.def( "create_directional_light", &create_directional_light );
	// cameras
	m.def( "create_perspective_camera", &create_perspective_camera );
	m.def( "create_sphere_camera", &create_sphere_camera );
	// integrators
	m.def( "create_simplept_integrator", &create_simplept_integrator );
	// volumes
	m.def( "create_volume", &create_volume );
	// fields
	m.def( "create_constant_field3d", &create_constant_field3d );
	// pnsolution
	m.def( "load_pnsolution", &load_pnsolution );
	m.def( "compute_fluence_pnsolution", &compute_fluence_pnsolution );
	// image
	m.def( "load_image", &load_image );
	m.def( "create_image_sampler", &create_image_sampler );
	m.def( "blur_image", &blur_image )
	;



}
