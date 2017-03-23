#include <cuda.h>
#include <builtin_types.h>

#include <stdio.h>

#include <cuda/CudaData.h>
#include <cuda/CudaEnvMap.h>
#include <cuda/CudaPixelBuffer.h>
#include <cuda/CudaVolume.h>
#include <cuda/CudaLight.h>



// returns sigma_t at sampled position (is invalid when we exceeded maxt)
DEVICE float delta_tracking_2d( const CudaVolume2D* volume, const cumath::Ray2f& ray, double maxt, int component, cumath::RNG& rng )
{
	float sigma_t_max = volume->getMaxExtinction()[component];

	float t = 0.0;
	while(true)
	{
		float step = -log( 1.0-rng.randomFloat() )/sigma_t_max;
		t += step;

		if(t>= maxt)
			break;

		cumath::V3f sigma_t = volume->evalExtinction(ray.getPosition(t), ray.d);

		// russian roulette
		if(rng.randomFloat()<sigma_t[component]/sigma_t_max)
			break;
	}

	return t;
}

struct TraceInfo2D
{
	CudaVolume2D* volume;
	bool debug;

	DEVICE cumath::V3f integrator_sample_transmittance( const cumath::Ray2f& ray, float maxt, cumath::RNG& rng )
	{
		float distance = delta_tracking_2d( volume, ray, maxt, 0, rng );

		if( distance < maxt )
			// since we return 0, we dont need to produce the pdf in the denominator
			return cumath::V3f(0.0, 0.0, 0.0);

		//NB: transmittance sampling pdf cancels out with the transmittance term

		return cumath::V3f(1.0, 1.0, 1.0);
	}
};


extern "C"
__global__ void pathtrace_fluencefield_2d_sample(   int sample,
													CudaPixelBuffer* fluencegrid,
													CudaVolume2D* volume)
{
	int resx = fluencegrid->getWidth();
	int resy = fluencegrid->getHeight();

	int x = blockIdx.x * blockDim.x + threadIdx.x;
	int y = blockIdx.y * blockDim.y + threadIdx.y;

	if((x<resx)&&(y<resy))
	{
		bool debug = false;

		//if((x == 256)&&(y==256))
		//	debug = true;

		cumath::RNG rng( (sample+1)*cumath::randhash(x*y) );


		// work out voxelcenter in worldspace
		cumath::V2f pVS( x+0.5, y+0.5 );
		cumath::V2f pLS = cumath::V2f( pVS.x/resx,
									   pVS.y/resy);
		cumath::V2f pWS = pLS*volume->getLocalToWorld();


		// compute direct light contribution ---
		// == scene_sample_attenuated_directlight in 2d
		// for now assume unit power point light source at 0,0,0
		cumath::V2f lightPos(0.0f, 0.0f);
		float lightDistance;
		cumath::V2f lightDir = (lightPos-pWS).normalized(lightDistance);
		cumath::Ray2f lightRay(pWS, lightDir, 0.0f, lightDistance);
		cumath::V3f light_over_pdf = 1.0/(2.0*MATH_PI*lightDistance);

		// sample transmittance between pWS and light position
		TraceInfo2D ti;
		ti.volume = volume;
		ti.debug = debug;
		cumath::V3f transmittance_over_pdf = ti.integrator_sample_transmittance( lightRay, lightDistance, rng );

		// final sample
		cumath::V3f attenuated_light_over_pdf = light_over_pdf*transmittance_over_pdf;


		cumath::V3f& c = fluencegrid->lvalue(x, y);
		c += (attenuated_light_over_pdf-c)/float(sample+1);
	}
}




void pathtrace_fluencefield_2d( CudaPtr<CudaPixelBuffer>* fluencegrid,
								CudaPtr<CudaVolume2D>* volume,
								//CudaPtr<CudaLight>* light,
								int numSamples )
{
	printf("pathtrace_fluencefield_2d\n");
	int resx = fluencegrid->host()->getWidth();
	int resy = fluencegrid->host()->getHeight();

	dim3 threadsPerBlock = dim3( 32, 32 );
	dim3 numBlocks = dim3(  int(ceil(double(resx)/double(threadsPerBlock.x))),
							int(ceil(double(resy)/double(threadsPerBlock.y))));


	for( int i=0;i<numSamples;++i )
	{
		std::cout << "pathtrace_fluencefield sample=" << i << std::endl;std::flush(std::cout);
		pathtrace_fluencefield_2d_sample<<<numBlocks, threadsPerBlock>>>(i, fluencegrid->device(), volume->device());

		cudaError_t err = cudaGetLastError();
		if (err != cudaSuccess)
			printf("Error: %s\n", cudaGetErrorString(err));

		CudaSafeCall( cudaDeviceSynchronize() );
	}
	cudaError_t err = cudaGetLastError();
	if (err != cudaSuccess)
		printf("Error: %s\n", cudaGetErrorString(err));

	CudaSafeCall( cudaDeviceSynchronize() );
}


