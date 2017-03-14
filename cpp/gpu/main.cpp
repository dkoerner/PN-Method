#include <iostream>
#include <string>


#include <util/string.h>
#include <util/envmap.h>
#include <util/timer.h>
#include <util/field.h>
#include <util/moexp.h>

#include <volume.h>
#include <shcache.h>

#include <cuda.h>
#include <builtin_types.h>
//#include <drvapi_error_string.h>

#include <cuda_runtime.h>
#include <cuda/CudaData.h>
#include <cuda/CudaEnvMap.h>
#include <cuda/CudaPixelBuffer.h>
#include <cuda/CudaVolume.h>
#include <cuda/CudaLight.h>


// Forward declare the function in the .cu file
void vectorAddition(const float* a, const float* b, float* c, int n);
void pathtrace_envmap( cumath::V3f pWS, CudaPtr<CudaEnvMap>* envmap, CudaPtr<CudaVolume>* volume, CudaPtr<CudaLight>* light, int numSamples );
void compute_sh_coefficients( CudaPtr<CudaPixelBuffer>* coords_ptr,
							  CudaPtr<CudaPixelBuffer>* f_ptr,
							  cumath::V3f* coeffs_host,
							  int order,
							  int numCoeffs);


void getCoordFromIndex( const V3i& resolution, int index, int& x, int& y, int& z )
{
	div_t divresult;

	divresult = div( index, resolution.x()*resolution.y() );

	z = divresult.quot;
	divresult = div( divresult.rem, resolution.x() );
	y = divresult.quot;
	x = divresult.rem;
}


void pathtrace_fluencefield( CudaPtr<CudaVoxelGrid<float>>* fluencegrid, CudaPtr<CudaVolume>* volume, CudaPtr<CudaLight>* light, int numSamples );
void generate_fluencefield( const std::string& filename, Transformd volume_localToWorld, CudaPtr<CudaVolume> cuvolume_ptr, CudaPtr<CudaLight> culight_ptr)
{
	int numSamples = 1000;
	V3i resolution(64, 64, 64);

	CudaVoxelGrid<float> cufluencegrid;
	CudaPtr<CudaVoxelGrid<float>> cufluencegrid_ptr(&cufluencegrid);
	cufluencegrid.initialize( cumath::V3i(resolution.x(), resolution.y(), resolution.z()));


	pathtrace_fluencefield(&cufluencegrid_ptr, &cuvolume_ptr, &culight_ptr, numSamples);

	std::vector<float> data( resolution.x()*resolution.y()*resolution.z() );
	cufluencegrid_ptr.host()->download( (unsigned char*)data.data() );
	field::VoxelGridField<float>::Ptr fluencegrid = std::make_shared<field::VoxelGridField<float>>(resolution, data.data());

	Fieldf::Ptr fluencefield = field::xform<float>(fluencegrid, volume_localToWorld);
	field::write( filename, fluencefield, resolution, volume_localToWorld );
}


int main()
{
	int deviceCount = 0;
	int cudaDevice = 0;
	char cudaDeviceName [100];

	cuInit(0);
	cuDeviceGetCount(&deviceCount);
	cuDeviceGet(&cudaDevice, 0);
	cuDeviceGetName(cudaDeviceName, 100, cudaDevice);
	std::cout << "Number of devices: " << deviceCount << std::endl;
	std::cout << "Device name:" << cudaDeviceName << std::endl;




	// setup volume ---------------

	//
	CudaVoxelGrid<float> cuvoxelgrid;
	CudaPtr<CudaVoxelGrid<float>> cuvoxelgrid_ptr(&cuvoxelgrid);

	//Volume::Ptr volume;
	Box3d volume_boundWS;
	float volume_sigma_t_max;
	Transformd volume_localToWorld;
	{
		std::string basePath = "c:/projects/visus/data";

		// volume ----
		double stepSize;
		field::VoxelGridField<float, float>::Ptr density = field::bgeo<float>(basePath + "/datasets/nebulae200.bgeo", &volume_localToWorld, &stepSize);

		volume_sigma_t_max = field_maximum<float>(density->grid);

		// setup localspace bounding box
		Box3d bboxLS;
		bboxLS.reset();
		bboxLS.min = P3d(0.0f,0.0f,0.0f);
		bboxLS.max = P3d(1.0f,1.0f,1.0f);

		// compute worldspace boundingbox
		volume_boundWS.reset();
		for( int i=0;i<8;++i )
			volume_boundWS.expandBy( volume_localToWorld*bboxLS.getCorner(i) );

		// initialize cudavoxelgrid
		V3i resolution = density->grid.getResolution();

		cuvoxelgrid.initialize( cumath::V3i(resolution.x(), resolution.y(), resolution.z()),
								density->grid.getRawPointer());
	}


	// ...
	CudaVolume cuvolume;
	cuvolume.m_boundWS = cumath::Box3f( cumath::V3f(volume_boundWS.min.x(), volume_boundWS.min.y(), volume_boundWS.min.z()),
										cumath::V3f(volume_boundWS.max.x(), volume_boundWS.max.y(), volume_boundWS.max.z()) );

	cuvolume.m_sigma_t = cuvoxelgrid_ptr.device();
	cuvolume.m_sigma_t_max = cumath::V3f(volume_sigma_t_max);
	M44d l2w = volume_localToWorld.getMatrix();
	M44d l2w_inv = volume_localToWorld.getInverseMatrix();
	for( int j=0;j<4;++j )
		for( int i=0;i<4;++i )
		{
			// NB: we transpose the matrix here because cumath uses row vectors
			cuvolume.m_localToWorld.matrix.m[j][i] = float(l2w.coeff(i,j));
			cuvolume.m_localToWorld.inverse.m[j][i] = float(l2w_inv.coeff(i,j));
		}

	CudaPtr<CudaVolume> cuvolume_ptr(&cuvolume);


	// setup light ----------------
	CudaLight culight;
	CudaPtr<CudaLight> culight_ptr(&culight);
	culight.m_direction = cumath::V3f(0.0, -1.0, 0.0);
	culight.m_radiance = cumath::V3f(1.0, 1.0, 1.0);
	culight.m_distance = volume_boundWS.getExtents().norm();





	generate_fluencefield("nebulae_fluence.bgeo", volume_localToWorld, cuvolume_ptr, culight_ptr);
	//generate_shcache(volume_localToWorld, cuvolume_ptr, culight_ptr);







	return 0;
}





void generate_shcache( Transformd volume_localToWorld, CudaPtr<CudaVolume> cuvolume_ptr, CudaPtr<CudaLight> culight_ptr )
{
	V2i envmap_resolution(128, 64);
	V3i sampling_resolution(64,64,64);
	//V3i sampling_resolution(1);
	int order = 30;
	int numVoxels = sampling_resolution.x()*sampling_resolution.y()*sampling_resolution.z();


	// setup envmap which is being sampled -----------
	CudaEnvMap cumap;
	CudaPtr<CudaEnvMap> cumap_ptr(&cumap);
	cumap.initialize(envmap_resolution.x(), envmap_resolution.y());


	// Buffer for sh coefficients
	CudaPixelBuffer cushcoeffs;
	CudaPtr<CudaPixelBuffer> cushcoeffs_ptr(&cushcoeffs);


	// environment map of the current voxel...
	EnvMap map(envmap_resolution);

	// shcache which holds the final sh coefficients for all voxels...
	SHCache shcache(order, sampling_resolution, volume_localToWorld);


	// go for each voxel...
	Timer timer;
	timer.start();
	for( int voxel=0;voxel<numVoxels;++voxel )
	{
		// clear envmap on device to have a clean plate for the next voxel
		cumap_ptr.host()->clear();

		// find voxel center in world space
		int i, j, k;
		getCoordFromIndex(sampling_resolution, voxel, i, j, k);

		P3d pLS((i+0.5)/sampling_resolution.x(),
				(j+0.5)/sampling_resolution.y(),
				(k+0.5)/sampling_resolution.z());
		P3d pWS = volume_localToWorld*pLS;
		cumath::V3f p(pWS.x(), pWS.y(), pWS.z());

		// kick of rendering
		if(0)
		{
			int numSamples = 80;
			pathtrace_envmap( p, &cumap_ptr, &cuvolume_ptr, &culight_ptr, numSamples );
		}

		// download result and save to disk
		if(0)
		{
			cumap.download((unsigned char*)map.bitmap().data());
			//std::string filename = "nebulae_envmaps/envmap_$0.exr";
			std::string filename = "test_$0.exr";
			filename = replace(filename, "$0", toString(voxel));
			map.bitmap().saveEXR(filename);
		}

		// load precomputed envmap per voxel
		if(0)
		{
			//std::string filename = "test_$0.exr";
			std::string filename = "nebulae_envmaps/envmap_$0.exr";
			filename = replace(filename, "$0", toString(voxel));
			// load rendered images
			map.bitmap() = Bitmap(filename);
		}

		// apply filter
		if(0)
		{
			double filter_width = 3.0;
			map.filter(gaussFilter(filter_width));
		}

		// compute sh coefficients
		if(0)
		{
			int numCoeffs = moexp::numSHCoefficients(order);


			moexp::SphericalFunction<Color3f> fun = [&]( double theta, double phi ) -> Color3f
			{
				return map.eval(theta, phi);
			};

			//std::unique_ptr<std::vector<Color3f>> coeffs_gt = moexp::project_Y_real( order, fun, 150000 );


			double min_theta = 0;
			double max_theta = M_PI;
			double min_phi = 0;
			double max_phi = 2.0*M_PI;

			int resolution_phi=map.bitmap().cols(); // resolution along phi angle
			int resolution_theta=map.bitmap().rows(); // resolution along theta angle
			double dtheta = (max_theta-min_theta)/double(resolution_theta);
			double dphi = (max_phi-min_phi)/double(resolution_phi);

			//std::vector<Color3f> coeffs(numCoeffs, Color3f(0.0f, 0.0f, 0.0f));
			//Color3f* coeffs_ptr = coeffs.data();
			Color3f* coeffs_ptr = shcache.getCoefficients(voxel);

			/*
			// cpu version for computing moments ---
			double pixel_area = dtheta*dphi;

			for( int t=0;t<resolution_theta;++t )
				for( int p=0;p<resolution_phi;++p )
				{
					double phi = dphi*p;
					double theta = dtheta*t;

					Color3f f = map.eval(theta, phi);

					for( int l=0;l<=order;++l )
						for( int m=-l;m<=l;++m )
						{
							double sh = moexp::Y_real(l, m, theta, phi);
							Color3f sample = f*sh;
							coeffs[moexp::shIndex(l, m)]+=sample*pixel_area*std::sin(theta);
						}
				}
			*/

			// gpu version for computing sh coefficients ---
			EnvMap input_coords(V2i(resolution_theta, resolution_phi));
			EnvMap input_f(V2i(resolution_theta, resolution_phi));
			for( int t=0;t<resolution_theta;++t )
				for( int p=0;p<resolution_phi;++p )
				{
					double phi = dphi*p;
					double theta = dtheta*t;

					Color3f f = map.eval(theta, phi);
					input_coords.bitmap().coeffRef( p, t ) = Color3f(theta, phi, 0.0f);
					input_f.bitmap().coeffRef( p, t ) = f;
				}


			CudaPixelBuffer gpu_input_coords;
			gpu_input_coords.initialize(resolution_theta, resolution_phi, (unsigned char*)input_coords.bitmap().data());
			CudaPixelBuffer gpu_input_f;
			gpu_input_f.initialize(resolution_theta, resolution_phi, (unsigned char*)input_f.bitmap().data());

			CudaPtr<CudaPixelBuffer> gpu_input_coords_ptr(&gpu_input_coords);
			CudaPtr<CudaPixelBuffer> gpu_input_f_ptr(&gpu_input_f);
			compute_sh_coefficients( &gpu_input_coords_ptr, &gpu_input_f_ptr, (cumath::V3f*)coeffs_ptr, order, numCoeffs );

			/*
			for( int l=0;l<=order;++l )
			{
				std::string filename = "sh_reconstruction_alt_$0.bgeo";
				filename = replace(filename, "$0", toString(l));

				rasterizeSphericalFunctionSphere(filename, [&](double theta, double phi)->Color3f
				{
					//return moexp::Y_real_sum<Color3f>(order, coeffs.data(), theta, phi);
					return moexp::Y_real_sum<Color3f>(l, coeffs.data(), theta, phi);
				}, 8.0);
			}
			*/

		}
	}
	timer.stop();

	std::cout << "path tracing took " << timer.elapsedSeconds()/double(numVoxels) << "s/voxel in average" << std::endl;

	//shcache.save("nebulae200.shcache");
}
