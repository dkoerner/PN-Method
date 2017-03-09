#include <cuda.h>
#include <builtin_types.h>

#include <stdio.h>

#include <cuda/CudaData.h>
#include <cuda/CudaEnvMap.h>
#include <cuda/CudaPixelBuffer.h>
#include <cuda/CudaVolume.h>
#include <cuda/CudaLight.h>

// has to come last
#include <cuda/pathtracing.cu.h>


extern "C"
__global__ void vectorAdditionCUDA(const float* a, const float* b, float* c, int n)
{
	int ii = blockDim.x * blockIdx.x + threadIdx.x;
	printf("test\n");
	if (ii < n)
		c[ii] = a[ii] + b[ii];
}

void vectorAddition(const float* a, const float* b, float* c, int n)
{
	float *a_cuda, *b_cuda, *c_cuda;
	unsigned int nBytes = sizeof(float) * n;
	int threadsPerBlock = 256;
	int blocksPerGrid   = (n + threadsPerBlock - 1) / threadsPerBlock;

	// allocate and copy memory into the device
	cudaMalloc((void **)& a_cuda, nBytes);
	cudaMalloc((void **)& b_cuda, nBytes);
	cudaMalloc((void **)& c_cuda, nBytes);
	cudaMemcpy(a_cuda, a, nBytes, cudaMemcpyHostToDevice);
	cudaMemcpy(b_cuda, b, nBytes, cudaMemcpyHostToDevice);

	vectorAdditionCUDA<<<blocksPerGrid, threadsPerBlock>>>(a_cuda, b_cuda, c_cuda, n);

	cudaError_t err = cudaGetLastError();
	if (err != cudaSuccess)
		printf("Error: %s\n", cudaGetErrorString(err));

	// load the answer back into the host
	cudaMemcpy(c, c_cuda, nBytes, cudaMemcpyDeviceToHost);

	cudaFree(a_cuda);
	cudaFree(b_cuda);
	cudaFree(c_cuda);
}



extern "C"
__global__ void pathtrace_envmap_sample( int sample,
										 cumath::V3f pWS,
										 CudaEnvMap* envmap,
										 CudaVolume* volume,
										 CudaLight* light /*int seed, CudaField<float> *grid_density, float density_max, float kcs, float albedo, cumath::V3f lightDir, float lightIntensity, CudaCamera<float> *cam, CudaPixelBuffer *pixels, int maxScatteringOrder*/ )
{
	int x = blockIdx.x * blockDim.x + threadIdx.x;
	int y = blockIdx.y * blockDim.y + threadIdx.y;
	int width = envmap->getWidth();
	int height = envmap->getHeight();

	if( (x<width)&&(y<height) )
	{
		bool debug = false;
		int index = y*width + x;

		if((x==128)&&(y==64))
		{
			//debug = true;
			//printf("pixels_device: %0xd     pixels_device: %0xd\n", pixels->m_data, pixels2);
		}

		cumath::RNG rng( (sample+1)*cumath::randhash(x*y) );
		cumath::Ray3<float> ray( pWS, envmap->xyToDirection(cumath::V2f(x, y)), 0.0f, 1000.0f );
		//cumath::V3f color = trace( ray );

		cumath::V3f f_over_pdf(0.0f);
		//double mint, maxt;
		//if( scene->volume->intersectBound(rq.ray, mint, maxt, rq.debug) )
		{
			// start tracing
			TraceInfo ti;
			ti.debug = debug;
			ti.volume = volume;
			ti.light = light;
			ti.depth = 0;
			ti.maxDepth = 50;
			ti.current_vertex = Vertex();
			ti.current_vertex.setPosition(pWS, cumath::V3f(0.0, 0.0, 0.0), cumath::V3f(0.0, 0.0, 0.0));
			ti.current_direction = ray.d;
			if(debug)
				printf("ti.current_direction=%f %f %f\n", ti.current_direction.x, ti.current_direction.y, ti.current_direction.z );

			//ti.scene = scene;
			//ti.debug = rq.debug;
			f_over_pdf = trace( ti, rng );
		}//else
		{
			// no intersection with the medium boundary
		}

		/*

		float L = pt_ms<float>( rng, grid_density, density_max, kcs, albedo, lightDir, lightIntensity, ray, maxScatteringOrder );
		*/

		// for testing: set a gradient
		//sample = cumath::V3f( float(x)/float(width), float(y)/float(height), 0.0f );
		//cumath::V3f rgb( rng.randomFloat() );
		//cumath::V3f sample( ray.d.x, ray.d.y, ray.d.z*0.5 );

		//cumath::V3f rgb( float(x)/float(width), float(y)/float(height), 0.0f );
		//cumath::V3f rgb( float(x)/float(width), 0.0f, 0.0f );
		cumath::V3f& c = envmap->lvalue(index);
		c += (f_over_pdf-c)/float(sample+1);
	}
}



void pathtrace_envmap( cumath::V3f pWS, CudaPtr<CudaEnvMap>* envmap, CudaPtr<CudaVolume>* volume, CudaPtr<CudaLight>* light, int numSamples )
{
	int width = envmap->host()->getWidth();
	int height= envmap->host()->getHeight();

	printf("pathtrace_envmap\n");
	dim3 threadsPerBlock = dim3( 16, 16 );
	dim3 numBlocks = dim3( int(ceil(double(width)/double(threadsPerBlock.x))),
							int(ceil(double(height)/double(threadsPerBlock.y))));

	for( int i=0;i<numSamples;++i )
		pathtrace_envmap_sample<<<numBlocks, threadsPerBlock>>>(i, pWS, envmap->device(), volume->device(), light->device());

	cudaError_t err = cudaGetLastError();
	if (err != cudaSuccess)
		printf("Error: %s\n", cudaGetErrorString(err));

	CudaSafeCall( cudaDeviceSynchronize() );
}











/*
void compute_sh_coefficients( CudaPtr<CudaEnvMap>* envmap, CudaPtr<CudaPixelBuffer>* coeffs, int numSamples )
{
	int width = envmap->host()->getWidth();
	int height= envmap->host()->getHeight();

	printf("compute_sh_coefficients\n");
	dim3 threadsPerBlock = dim3( 16, 16 );
	dim3 numBlocks = dim3( int(ceil(double(width)/double(threadsPerBlock.x))),
							int(ceil(double(height)/double(threadsPerBlock.y))));

	for( int i=0;i<numSamples;++i )
		pathtrace_envmap_sample<<<numBlocks, threadsPerBlock>>>(i, pWS, envmap->device(), volume->device(), light->device());

	cudaError_t err = cudaGetLastError();
	if (err != cudaSuccess)
		printf("Error: %s\n", cudaGetErrorString(err));

	CudaSafeCall( cudaDeviceSynchronize() );
}
*/
