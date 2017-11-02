


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



Eigen::VectorXd compute_fluence( Volume::Ptr volume, Light::Ptr light, Camera::Ptr camera, Integrator::Ptr integrator, const Eigen::MatrixXd& points )
{
	std::cout << "computing fluence:\n";
	std::cout << volume->toString();
	std::cout << light->toString();
	std::cout << camera->toString();
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
	scene.camera = camera.get();
	scene.integrator = integrator.get();
	scene.light = light.get();
	scene.volume = volume.get();

	GlobalComputeFluenceInfo gi;
	gi.scene = &scene;
	gi.fluence = &fluence;
	gi.sensor_points = &sensor_points;

	int numSamples = 10000;

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


Image::Ptr render( Volume::Ptr volume, Light::Ptr light, Camera::Ptr camera, Integrator::Ptr integrator )
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
	//gi.debug_pixel = V2i(256, 256);

	int numSamples = 1;

	Terminator terminator(numSamples);
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


void test( Camera::Ptr cam )
{
	P2d pRS(256.5, 256.5);
	Ray3d ray;
	cam->sampleRay(pRS, ray);

	std::cout << "ray=" << ray.toString() << std::endl;
	std::cout << "ray(3.0)=" << ray(3.0).toString() << std::endl;
}


PNSolution::Ptr load_pnsolution( const std::string& filename )
{
	return std::make_shared<PNSolution>(filename);
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
	;

	// Image ============================================================
	py::class_<Image, Image::Ptr> class_image(m, "Image");
	class_image
		.def( "save", &Image::save )
	;

	// Field3d ============================================================
	py::class_<Field3d, Field3d::Ptr> class_field3d(m, "Field3d");
	class_field3d
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
	.def("evalCoefficient",
	 [](PNSolution &m,
		const Eigen::Matrix<double, 3, 1>& pWS,
		int coeff_index)
	 {
		return m.evalCoefficient(pWS, coeff_index);
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

	// main render function
	m.def( "render", &render );
	m.def( "compute_fluence", &compute_fluence );
	m.def( "test", &test );
	// lights
	m.def( "create_point_light", &create_point_light );
	//m.def( "create_directional_light", &create_directional_light );
	// cameras
	m.def( "create_perspective_camera", &create_perspective_camera );
	// integrators
	m.def( "create_simplept_integrator", &create_simplept_integrator );
	// volumes
	m.def( "create_volume", &create_volume );
	// fields
	m.def( "create_constant_field3d", &create_constant_field3d );
	// pnsolution
	m.def( "load_pnsolution", &load_pnsolution );
	m.def( "compute_fluence_pnsolution", &compute_fluence_pnsolution );


}
