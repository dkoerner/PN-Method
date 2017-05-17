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





DEVICE double factorial( int x )
{
	double s = 1.0;
	for (int n = 2; n <= x; ++n)
	{
		s *= n;
	}
	return s;
}

DEVICE double doubleFactorial(int x)
{
	const double dbl_factorial_cache[16] = {1, 1, 2, 3, 8, 15, 48, 105,
												  384, 945, 3840, 10395, 46080,
												  135135, 645120, 2027025};

	if (x < 16)
	{
		return dbl_factorial_cache[x];
	}else
	{
		double s = 1.0;
		double n = x;
		while (n > 1.0)
		{
			s *= n;
			n -= 2.0;
		}
		return s;
	}
}

DEVICE double P_real( int l, int m, double x )
{
	// Compute Pmm(x) = (-1)^m(2m - 1)!!(1 - x^2)^(m/2), where !! is the double
	// factorial.
	double pmm = 1.0;
	// P00 is defined as 1.0, do don't evaluate Pmm unless we know m > 0
	if (m > 0)
	{
		double sign = (m % 2 == 0 ? 1 : -1);
		pmm = sign * doubleFactorial(2 * m - 1) * pow(1 - x * x, m / 2.0);
	}

	if (l == m) {
	  // Pml is the same as Pmm so there's no lifting to higher bands needed
	  return pmm;
	}

	// Compute Pmm+1(x) = x(2m + 1)Pmm(x)
	double pmm1 = x * (2 * m + 1) * pmm;
	if (l == m + 1) {
	  // Pml is the same as Pmm+1 so we are done as well
	  return pmm1;
	}

	// Use the last two computed bands to lift up to the next band until l is
	// reached, using the recurrence relationship:
	// Pml(x) = (x(2l - 1)Pml-1 - (l + m - 1)Pml-2) / (l - m)
	for (int n = m + 2; n <= l; n++)
	{
	  double pmn = (x * (2 * n - 1) * pmm1 - (n + m - 1) * pmm) / (n - m);
	  pmm = pmm1;
	  pmm1 = pmn;
	}
	// Pmm1 at the end of the above loop is equal to Pml
	return pmm1;
}

DEVICE double compute_Ylm( int l, int m, double theta, double phi )
{
	double kml = sqrt((2.0 * l + 1) * factorial(l - abs(m)) /
					  (4.0 * MATH_PI * factorial(l + abs(m))));
	if (m > 0)
	{
		return sqrt(2.0) * kml * cos(m * phi) * P_real(l, m, cos(theta));
	}else
	if (m < 0)
	{
		return sqrt(2.0) * kml * sin(-m * phi) * P_real(l, -m, cos(theta));
	}else
	{
		return kml * P_real(l, 0, cos(theta));
	}
}


extern "C"
__global__ void compute_sh_coefficients_kernel_evalYlm( int l, int m,
														CudaPixelBuffer* f,
														CudaPixelBuffer* coords,
														cumath::V3f* coeffmap)
{
	int width = f->getWidth();
	int height= f->getHeight();
	int numPixels = width*height;
	float dtheta = MATH_PI/float(height);
	float dphi = 2.0f*MATH_PI/float(width);
	float pixel_area = dtheta*dphi;

	int index = blockDim.x * blockIdx.x + threadIdx.x;

	if( index < numPixels )
	{
		const cumath::V3f& coord = coords->lvalue(index);
		float theta = coord.x;
		float phi = coord.y;
		cumath::V3f value = f->lvalue(index);
		double Ylm = compute_Ylm(l, m, theta, phi);
		coeffmap[index] = cumath::V3f(float(Ylm*value.x),
									  float(Ylm*value.y),
									  float(Ylm*value.z))*pixel_area*sin(theta);
	}
}

extern "C"
__global__ void compute_sh_coefficients_kernel_reduce( int width, int height,
													   cumath::V3f** coeffmap_list,
													   cumath::V3f* coeff_list,
													   int numCoeffs)
{
	int coeff_index = blockDim.x * blockIdx.x + threadIdx.x;
	if( coeff_index < numCoeffs )
	{
		cumath::V3f sum(0.0f, 0.0f, 0.0f);
		int pix_index = 0;
		for( int i=0;i<height;++i )
			for( int j=0;j<width;++j, ++pix_index)
			{
				sum = sum+coeffmap_list[coeff_index][pix_index];
			}
		coeff_list[coeff_index] = sum;
	}
}

// coords_ptr -> contains for each pixel the theta, phi coordinate
// f_ptr -> contains for each pixel the function f, evaluated at theta, phi
// coeffs_host -> host list of final coefficients
// order -> sh order
// numCoeffs -> number of coefficients
void compute_sh_coefficients( CudaPtr<CudaPixelBuffer>* coords_ptr,
							  CudaPtr<CudaPixelBuffer>* f_ptr,
							  cumath::V3f* coeffs_host,
							  int order,
							  int numCoeffs)
{
	printf("compute_sh_coefficients\n");
	int width = f_ptr->host()->getWidth();
	int height= f_ptr->host()->getHeight();
	int numPixels = width*height;

	cumath::V3f* coeffs_device=0;
	cudaMalloc( &coeffs_device, numCoeffs*sizeof(cumath::V3f) );

	cumath::V3f** coeffmap_list_host = (cumath::V3f**)malloc( numCoeffs*sizeof(cumath::V3f*) );
	cumath::V3f** coeffmap_list_device=0;
	cudaMalloc( &coeffmap_list_device, numCoeffs*sizeof(cumath::V3f*) );

	for( int i=0;i<numCoeffs;++i )
		CudaSafeCall( cudaMalloc( &coeffmap_list_host[i], sizeof(cumath::V3f)*numPixels ) );
	CudaSafeCall( cudaMemcpy( coeffmap_list_device, coeffmap_list_host, numCoeffs*sizeof(cumath::V3f*), cudaMemcpyHostToDevice ));


	for( int l=0;l<=order;++l )
		for( int m=-l;m<=l;++m )
		{
			int shIndex =  l * (l + 1) + m;
			int threadsPerBlock = 256;
			int blocksPerGrid  = (numPixels + threadsPerBlock - 1) / threadsPerBlock;
			compute_sh_coefficients_kernel_evalYlm<<<blocksPerGrid, threadsPerBlock>>>( l, m,
																						f_ptr->device(),
																						coords_ptr->device(),
																						coeffmap_list_host[shIndex] );
		}

	// now reduce each coefficient map into a single coefficient (accumulate the sums)
	{
		int threadsPerBlock = 256;
		int blocksPerGrid  = (numCoeffs + threadsPerBlock - 1) / threadsPerBlock;
		compute_sh_coefficients_kernel_reduce<<<blocksPerGrid, threadsPerBlock>>>( f_ptr->host()->getWidth(), f_ptr->host()->getHeight(),
																					coeffmap_list_device,
																				   coeffs_device,
																				   numCoeffs);
	}

	// TODO: copy coeffs from device to host
	CudaSafeCall( cudaMemcpy( coeffs_host, coeffs_device, numCoeffs*sizeof(cumath::V3f), cudaMemcpyDeviceToHost));

	cudaError_t err = cudaGetLastError();
	if (err != cudaSuccess)
		printf("Error: %s\n", cudaGetErrorString(err));

	CudaSafeCall( cudaDeviceSynchronize() );

	for( int i=0;i<numCoeffs;++i )
		CudaSafeCall( cudaFree(coeffmap_list_host[i]) );
	CudaSafeCall( cudaFree(coeffmap_list_device) );

	CudaSafeCall( cudaFree(coeffs_device) );

	free(coeffmap_list_host);
}


// ==================================================================


extern "C"
__global__ void pathtrace_fluencefield_sample(  int sample,
												CudaVoxelGrid<float>* fluencegrid,
												CudaVolume* volume,
												CudaLight* light )
{
	cumath::V3f resolution = fluencegrid->m_resolution;
	int x = blockIdx.x * blockDim.x + threadIdx.x;
	int y = blockIdx.y * blockDim.y + threadIdx.y;
	int z = blockIdx.z * blockDim.z + threadIdx.z;

	if( (x<resolution.x)&&(y<resolution.y)&&(z<resolution.z) )
	{
		bool debug = false;
		//int index = z*resolution.x*resolution.y + y*resolution.x + x;

		cumath::RNG rng( (sample+1)*cumath::randhash(x*y*z) );


		// work out voxelcenter in worldspace
		cumath::V3f pVS( x+0.5, y+0.5, z+0.5 );
		cumath::V3f pLS = fluencegrid->voxelToLocal(pVS);
		cumath::V3f pWS = pLS*volume->getLocalToWorld();


		// pick random direction
		float d_pdf;
		cumath::V3f d = cumath::sampleSphere(rng, &d_pdf);

		// assemble ray
		cumath::Ray3<float> ray( pWS, d, 0.0f, 1000.0f );

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
			ti.maxDepth = 100;
			ti.current_vertex = Vertex();
			ti.current_vertex.setPosition(pWS, cumath::V3f(0.0, 0.0, 0.0), cumath::V3f(0.0, 0.0, 0.0));
			ti.current_direction = ray.d;
			if(debug)
				printf("ti.current_direction=%f %f %f\n", ti.current_direction.x, ti.current_direction.y, ti.current_direction.z );

			//ti.scene = scene;
			//ti.debug = rq.debug;
			f_over_pdf = trace( ti, rng )/d_pdf;
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
		//cumath::V3f& c = envmap->lvalue(index);
		//c += (f_over_pdf-c)/float(sample+1);
		float& c = fluencegrid->lvalue(x, y, z);
		c += (f_over_pdf.x-c)/float(sample+1);
	}
}




void pathtrace_fluencefield( CudaPtr<CudaVoxelGrid<float>>* fluencegrid,
							 CudaPtr<CudaVolume>* volume,
							 CudaPtr<CudaLight>* light,
							 int numSamples )
{
	cumath::V3i resolution = fluencegrid->host()->m_resolution;

	printf("pathtrace_fluencefield\n");
	dim3 threadsPerBlock = dim3( 8, 8, 8 );
	dim3 numBlocks = dim3(  int(ceil(double(resolution.x)/double(threadsPerBlock.x))),
							int(ceil(double(resolution.y)/double(threadsPerBlock.y))),
							int(ceil(double(resolution.z)/double(threadsPerBlock.z))));


	for( int i=0;i<numSamples;++i )
	{
		std::cout << "pathtrace_fluencefield sample=" << i << std::endl;std::flush(std::cout);
		pathtrace_fluencefield_sample<<<numBlocks, threadsPerBlock>>>(i, fluencegrid->device(), volume->device(), light->device());

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
// ===================================================================

DEVICE cumath::V3f complex_mul( const cumath::V3f& a, const cumath::V3f& b )
{
	return cumath::V3f( a.x*b.x - a.y*b.y, a.x*b.y + a.y*b.x, 0.0f );
}


extern "C"
__global__ void compute_sh_coefficients_phase_kernel( CudaPhaseFunctionSGGX3D* phase, // phase function to be projected
													  CudaVoxelGrid<cumath::V3f>* Ylm, // precomputed Ylm
													  CudaVoxelGrid<cumath::V3f>* result // projected sh coefficients for each wi-pixel
													  )
{
	int phi_index = blockIdx.x * blockDim.x + threadIdx.x;
	int theta_index = blockIdx.y * blockDim.y + threadIdx.y;
	int sh_index = blockIdx.z * blockDim.z + threadIdx.z;

	int resolution_phi = Ylm->m_resolution.x;
	int resolution_theta = Ylm->m_resolution.y;
	int numCoeffs = Ylm->m_resolution.z;

	float dtheta = MATH_PIf/float(resolution_theta);
	float dphi = 2.0f*MATH_PIf/float(resolution_phi);
	if( (phi_index<resolution_phi)&&(theta_index<resolution_theta)&&(sh_index<numCoeffs) )
	{
		// compute angle theta_i, phi_i associated with the current pixel position from discretized wi
		float theta_i = dtheta*theta_index;
		float phi_i = dphi*phi_index;
		cumath::V3d wi = cumath::sphericalDirection(theta_i, phi_i);

		// iterate all pixels wich discretize wo
		cumath::V3f sum(0.0f, 0.0f, 0.0f);
		for( int i=0;i<resolution_phi;++i )
			for( int j=0;j<resolution_theta;++j )
			{
				// compute theta_o, phi_o for current pixel
				float theta_o = dtheta*j;
				float phi_o = dphi*i;
				cumath::V3d wo = cumath::sphericalDirection(theta_o, phi_o);

				// evaluate phase function at wi, wo
				cumath::V3f p(phase->eval(wi, wo), 0.0f, 0.0f);
				//printf("p(%f %f %f,  %f %f %f)=%f %f\n", wi.x, wi.y, wi.z, wo.x, wo.y, wo.z, p.x, p.y );

				// multiply with precomputed Ylm at that position (and jacobian and delta area)
				// and add to the current sum
				sum += complex_mul( p, Ylm->sample(i, j, sh_index) )*sin(theta_o)*dtheta*dphi;
			}

		// store result in buffer
		result->lvalue(phi_index, theta_index, sh_index) = sum;
	}
}

DEVICE void getCoordFromIndex( int index, int resX, int &x, int &y )
{
	int numer = index;
	int denom = resX;
	int quot = numer / denom;
	int rem = numer % denom;


	if(numer >= 0 && rem < 0)
	{
		++quot;
		rem -= denom;
	}

	y = quot;
	x = rem;
}

extern "C"
__global__ void compute_sh_coefficients_phase_kernel_final_project(CudaVoxelGrid<cumath::V3f>* Ylm, // precomputed Ylm
																   CudaVoxelGrid<cumath::V3f>* shCoeffs, // computed sh coefficients for each wi-pixel
																   CudaPixelBuffer* result )
{
	int index = blockDim.x * blockIdx.x + threadIdx.x;
	int numCoeffs = Ylm->m_resolution.z;
	if( index < numCoeffs*numCoeffs )
	{
		// work out 2d i,j coordinate within matrix from flattened 1d index
		int i, j;
		getCoordFromIndex( index, numCoeffs, i, j );

		// we have discretized the sh coefficients of the phase function on a theta-phi grid,
		// where theta and phi were the wi vectors and wo was used for the projection
		//
		// here we now compute the sh projection of the discretized spherical function established
		// by the i-th coefficient (which varies over theta,phi) onto the sh basis function represented by j
		cumath::V3f sum(0.0f, 0.0f, 0.0f);
		int resolution_phi = Ylm->m_resolution.x;
		int resolution_theta = Ylm->m_resolution.y;
		float dtheta = MATH_PIf/float(resolution_theta);
		float dphi = 2.0f*MATH_PIf/float(resolution_phi);
		for( int phi_index=0;phi_index<resolution_phi;++phi_index )
			for( int theta_index=0;theta_index<resolution_theta;++theta_index )
			{
				// compute theta_i, phi_i for current pixel
				float theta_i = dtheta*theta_index;
				//float phi_i = dphi*phi_index;

				// retrieve i-th coefficient (which varies across pixels)
				cumath::V3f coeff = shCoeffs->sample( phi_index, theta_index, i );

				//cumath::V3f coeff(0.079577, 0.0f, 0.0f);
				//printf("coeff=%f %f %f\n", coeff.x, coeff.y, coeff.z);

				// multiply with precomputed Y_j at that position (and jacobian and delta area)
				// and add to the current sum
				sum += complex_mul(coeff, Ylm->sample(phi_index, theta_index, j) )*sin(theta_i)*dtheta*dphi;
			}

		result->lvalue(i, j) = sum;
		//result->lvalue(i, j) = cumath::V3f(1.0f, 2.0f, 3.0f);
	}
}

void compute_sh_coefficients_phase( CudaPtr<CudaPhaseFunctionSGGX3D>* cuPhase, // the phase function to be projected
									CudaPtr<CudaVoxelGrid<cumath::V3f>>* cuYlm, // precomputed sh basis coefficients
									CudaPtr<CudaPixelBuffer>* cuP // the resulting coefficient matrix
									)
{
	cumath::V3i resolution = cuYlm->host()->m_resolution;
	int numCoeffs = resolution.z;

	// create temporary datastructure which will hold the sh coefficients of the phase function for all wi
	// the final result will project each coefficient into sh again (giving rise to a matrix of coefficients)
	CudaVoxelGrid<cumath::V3f> temp;
	temp.initialize( cuYlm->host()->m_resolution );
	CudaPtr<CudaVoxelGrid<cumath::V3f>> temp_ptr(&temp);


	// for each wi (given by a pixel coordinate) compute the sh coefficients for the phase function evaluated
	// at wi and the angle over which the sh projection is carried out
	// this gives a set of sh coefficients for each theta,phi-pixel
	{
		dim3 threadsPerBlock = dim3( 8, 8, 8 );
		dim3 numBlocks = dim3(  int(ceil(double(resolution.x)/double(threadsPerBlock.x))),
								int(ceil(double(resolution.y)/double(threadsPerBlock.y))),
								int(ceil(double(resolution.z)/double(threadsPerBlock.z))));


		std::cout << "computing sh coefficients" << std::endl;std::flush(std::cout);
		compute_sh_coefficients_phase_kernel<<<numBlocks, threadsPerBlock>>>(cuPhase->device(),
																			 cuYlm->device(),
																			 temp_ptr.device());

		cudaError_t err = cudaGetLastError();
		if (err != cudaSuccess)
			printf("Error: %s\n", cudaGetErrorString(err));

		CudaSafeCall( cudaDeviceSynchronize() );
	}

	// The n-th coefficient now forms a spherical function over theta-phi. We project this function into
	// a set of spherical harmonics coefficients. This gives rise to a nxn coefficient matrix where n
	// is the number of SH coefficients for given order
	// we launch n x n threads where each thread computes the sh coefficient i of the expanded coefficient
	// function j
	{
		int threadsPerBlock = 256;
		// NB:launching numCoeffs x numCoeffs threads...a thread or each matrix component
		int blocksPerGrid  = (numCoeffs*numCoeffs + threadsPerBlock - 1) / threadsPerBlock;
		std::cout << "computing final coefficient matrix" << std::endl;std::flush(std::cout);
		compute_sh_coefficients_phase_kernel_final_project<<<blocksPerGrid, threadsPerBlock>>>(cuYlm->device(),
																							   temp_ptr.device(),
																							   cuP->device());

		cudaError_t err = cudaGetLastError();
		if (err != cudaSuccess)
			printf("Error: %s\n", cudaGetErrorString(err));

		CudaSafeCall( cudaDeviceSynchronize() );
	}

}













