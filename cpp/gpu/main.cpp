#include <iostream>
#include <string>


#include <util/string.h>
#include <util/envmap.h>
#include <util/timer.h>
#include <util/field.h>
#include <util/moexp.h>

#include <volume.h>

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


void getCoordFromIndex( const V3i& resolution, int index, int& x, int& y, int& z )
{
	div_t divresult;

	divresult = div( index, resolution.x()*resolution.y() );

	z = divresult.quot;
	divresult = div( divresult.rem, resolution.x() );
	y = divresult.quot;
	x = divresult.rem;
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



	//
	V2i envmap_resolution(128, 64);
	V3i sampling_resolution(64,64,64);
	int numVoxels = sampling_resolution.x()*sampling_resolution.y()*sampling_resolution.z();



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

	//P3d test = volume_localToWorld*P3d(0.5, 0.5, 0.5);
	//cumath::V3f test2 = cumath::V3f(0.5, 0.5, 0.5)*cuvolume.m_localToWorld.matrix;
	//return 0;
	CudaPtr<CudaVolume> cuvolume_ptr(&cuvolume);


	// setup envmap which is being sampled -----------
	CudaEnvMap cumap;
	CudaPtr<CudaEnvMap> cumap_ptr(&cumap);
	cumap.initialize(envmap_resolution.x(), envmap_resolution.y());


	// setup light ----------------
	CudaLight culight;
	CudaPtr<CudaLight> culight_ptr(&culight);
	culight.m_direction = cumath::V3f(0.0, -1.0, 0.0);
	culight.m_radiance = cumath::V3f(1.0, 1.0, 1.0);
	culight.m_distance = volume_boundWS.getExtents().norm();


	// Buffer for sh coefficients
	CudaPixelBuffer cushcoeffs;
	CudaPtr<CudaPixelBuffer> cushcoeffs_ptr(&cushcoeffs);


	// kick of rendering ------------
	EnvMap map(envmap_resolution);

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
		{
			int numSamples = 80;
			pathtrace_envmap( p, &cumap_ptr, &cuvolume_ptr, &culight_ptr, numSamples );
		}

		// download result and save to disk
		{
			cumap.download((unsigned char*)map.bitmap().data());
			std::string filename = "nebulae_envmaps/envmap_$0.exr";
			filename = replace(filename, "$0", toString(voxel));
			map.bitmap().saveEXR(filename);
		}

		// apply filter
		if(0)
		{
			double filter_width = 3.0;
			map.bitmap() = filter(map.bitmap(), gaussFilter(filter_width));
			std::string filename("test_$0_filtered.exr");
			filename = replace(filename, "$0", toString(voxel));
			map.bitmap().saveEXR(filename);

			// upload filtered image to device
			//cumap.upload( (unsigned char*)map.bitmap().data() );
		}

		// compute sh coefficients
		if(0)
		{
			int order = 30;
			int numCoeffs = moexp::numSHCoefficients(order);
			std::cout << numCoeffs << std::endl;

			//double min_theta = std::numeric_limits<double>::max();
			//double max_theta = std::numeric_limits<double>::min();

			moexp::SphericalFunction<Color3f> fun = [&]( double theta, double phi ) -> Color3f
			{
				return map.eval(theta, phi);
				//return 1.0;
			};

			//std::unique_ptr<std::vector<Color3f>> coeffs_gt = moexp::project_Y_real( order, fun, 150000 );


			{
				V2i res(1024, 1024);
				EnvMap map2(res);

				for( int i=0;i<res.y();++i )
					for( int j=0;j<res.x();++j )
					{
						double u = double(j)/double(res.x());
						double v = double(i)/double(res.y());
						map2.bitmap().coeffRef( i, j ) = Color3f(u, 0.0f, 0.0f);
					}
				map2.bitmap().saveEXR("uv.exr");
				rasterizeSphericalFunctionSphere("uv.bgeo", [&](double theta, double phi)->Color3f
				{
					return map2.eval(theta, phi);
				});
			}


			/*
			std::vector<Color3f> coeffs(numCoeffs, Color3f(0.0f, 0.0f, 0.0f));

			double min_theta = 0;
			double max_theta = M_PI;
			double min_phi = 0;
			double max_phi = 2.0*M_PI;

			V2i res(1024, 1024);
			EnvMap map2(res);
			rasterizeSphericalFunctionMap2(map2, fun);
			map2.bitmap().saveEXR("test2.exr");
			double dtheta = (max_theta-min_theta)/double(res.y());
			double dphi = (max_phi-min_phi)/double(res.x());
			double pixel_area = dtheta*dphi;

			for( int i=0;i<res.y();++i )
				for( int j=0;j<res.x();++j )
				{
					//P2d theta_phi = sphericalCoordinates( map2.xyToDirection(P2d(j+0.5, i+0.5)));
					//double theta = theta_phi.x();
					//double phi = theta_phi.y();
					double phi = dphi*j;
					double theta = dtheta*i;



					Color3f f = map2.bitmap().coeffRef(i, j);

					for( int l=0;l<=order;++l )
						for( int m=-l;m<=l;++m )
						{
							double sh = moexp::Y_real(l, m, theta, phi);
							Color3f sample = f*sh;
							coeffs[moexp::shIndex(l, m)]+=sample*pixel_area*std::sin(theta);
						}
				}
			*/



			for( int l=0;l<=order;++l )
			{
				for( int m=-l;m<=l;++m )
				{
					//std::cout << "l=" << l  << " " << (*coeffs_gt)[moexp::shIndex(l,m)].r() << " " << coeffs[moexp::shIndex(l,m)].r() << std::endl;
					//std::cout << "l=" << l  << " " << coeffs[moexp::shIndex(l,m)].r() << std::endl;
				}
			}

			//std::cout << "check=" << sum << " " << 4.0*M_PI << std::endl;


			//rasterizeSphericalFunctionMap( "test.exr", fun, res );
			//rasterizeSphericalFunctionSphere("groundtruth.bgeo", fun, 8.0);
			/*
			rasterizeSphericalFunctionSphere("sh_reconstruction_gt.bgeo", [&](double theta, double phi)->Color3f
			{
				return moexp::Y_real_sum<Color3f>(order, coeffs_gt->data(), theta, phi);
			}, 8.0);
			*/
			/*
			rasterizeSphericalFunctionSphere("sh_reconstruction_alt.bgeo", [&](double theta, double phi)->Color3f
			{
				return moexp::Y_real_sum<Color3f>(order, coeffs.data(), theta, phi);
			}, 8.0);
			*/

			rasterizeSphericalFunctionSphere("uv.bgeo", [&](double theta, double phi)->Color3f
			{
				return Color3f(phi/(2.0*M_PI), 0.0f, 0.0f);
			});

			//std::cout << integral << " " << 4.0*M_PI << std::endl;


			/*
			// check
			{
				V2i res(1024, 1024);
				sh::DefaultImage env_sh(res.x(), res.y());
				for( int i=0;i<res.y();++i )
					for( int j=0;j<res.x();++j )
					{
						Eigen::Array3f p;
						p.coeffRef(0) = 1.0;
						p.coeffRef(1) = 1.0;
						p.coeffRef(2) = 1.0;
						env_sh.SetPixel( j, i, p );
					}
				std::unique_ptr<std::vector<Eigen::Array3f>> coeff_sh = sh::ProjectEnvironment( order, env_sh );
				for( int l=0;l<=order;++l )
					for( int m=-l;m<=l;++m )
						std::cout << "l=" <<" " << (*coeff_sh)[sh::GetIndex(l,m)].coeff(0) << " " << std::endl;
			}
			*/

			/*
			//int numSamples = 100;
			//compute_sh_coefficients( &cumap_ptr, &cushcoeffs_ptr, numSamples );
			for( int l=0;l<=order;++l )
				for( int m=-l;m<=l;++m )
				{
					//int l = 1;
					//int m = -1;
					double phi = 0.4*M_PI/2.0;
					double theta = 0.6*M_PI;
					double sh_gt = sh::EvalSHSlow(l, m, phi, theta);
					double sh = moexp::Y_real(l, m, theta, phi);
					std::cout << "l=" << l << " m=" << m << " " << sh_gt << " " << sh << std::endl;
				}
			*/
		}
	}
	timer.stop();


	std::cout << "path tracing took " << timer.elapsedSeconds()/double(numVoxels) << "s/voxel in average" << std::endl;

	/*
	rasterizeSphericalFunctionSphere( "test.bgeo", [&](double theta, double phi) -> Color3f
	{
		return map.eval(theta, phi);
	});
	*/


	return 0;
}

