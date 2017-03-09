#pragma once
#include <cuda.h>
#include "cumath/cumath.h"


struct CudaEnvMap
{
	HOST CudaEnvMap() : m_width(0), m_height(0), m_data(0)
	{
	}
	HOST virtual ~CudaEnvMap()
	{
		if(m_data)
		{
			CudaSafeCall( cudaDeviceSynchronize() );
			CudaSafeCall( cudaFree( m_data ) );
		}
	}

	// initializes fields and uploads host data if needed
	HOST void initialize( int width, int height, unsigned char *host_ptr = 0 )
	{
		m_width = width;
		m_height = height;

		CudaSafeCall( cudaDeviceSynchronize() );
		CudaSafeCall( cudaFree( m_data ) );
		int size =  m_width*m_height*sizeof(cumath::V3f);

		CudaSafeCall( cudaMalloc( &m_data, size ) );

		if(host_ptr)
			CudaSafeCall( cudaMemcpy( m_data, host_ptr, size, cudaMemcpyHostToDevice ));
		else
			CudaSafeCall( cudaMemset( m_data, 0, size ));
		CudaSafeCall( cudaDeviceSynchronize() );
	}

	HOST void download( unsigned char *host_ptr )
	{
		int size =  m_width*m_height*sizeof(cumath::V3f);
		CudaSafeCall( cudaMemcpy( host_ptr, m_data, size, cudaMemcpyDeviceToHost ));
	}

	HOST void upload( unsigned char *host_ptr )
	{
		int size =  m_width*m_height*sizeof(cumath::V3f);
		CudaSafeCall( cudaMemcpy( m_data, host_ptr, size, cudaMemcpyHostToDevice ));
	}

	DEVICE cumath::V3f& lvalue( int index )
	{
		return m_data[index];
	}



	DEVICE cumath::V2f directionToUV( const cumath::V3f& d )const
	{
		// using formulas given in http://gl.ict.usc.edu/Data/HighResProbes/
		// with the difference that u=[0,1] (instead of [0,2]) and we negate z
		cumath::V2f uv((1.0+atan2(d.x, d.z)/CUDART_PI_F/2.0f),
					   cumath::safe_acos(d.y)/CUDART_PI_F);
		return uv;
	}
	DEVICE cumath::V3f uvToDirection( const cumath::V2f& uv )const
	{
		// using formulas given in http://gl.ict.usc.edu/Data/HighResProbes/
		// with the difference that u=[0,1] (instead of [0,2]) and we negate z
		// azimuthal angle
		double theta = CUDART_PI_F*(uv.x*2.0-1.0);
		// elevation angle
		double phi = CUDART_PI_F*uv.y;
		return cumath::V3f( sin(phi)*sin(theta), cos(phi), sin(phi)*cos(theta) );
	}

	DEVICE cumath::V2f uvToXY(const cumath::V2f& uv)const
	{
		cumath::V2f xy(
			(uv.x*(m_width-1)),
			(uv.y*(m_height-1))
			);
		return xy;
	}

	DEVICE cumath::V2f xyToUV(const cumath::V2f& xy)const
	{
		return cumath::V2f(
			(xy.x)/double(m_width-1),
			(xy.y)/double(m_height-1)
			);
	}

	DEVICE cumath::V3f xyToDirection( const cumath::V2f& xy )const
	{
		return uvToDirection( xyToUV(xy) );
	}

	DEVICE cumath::V2f directionToXY( const cumath::V3f& d )const
	{
		return uvToXY(directionToUV(d));
	}

	void clear()
	{
		CudaSafeCall( cudaDeviceSynchronize() );
		CudaSafeCall( cudaMemset( m_data, 0, m_width*m_height*sizeof(cumath::V3f) ));
	}

	HOST_DEVICE int getWidth()const
	{
		return m_width;
	}

	HOST_DEVICE int getHeight()const
	{
		return m_height;
	}

	int                         m_width;
	int                         m_height;
	cumath::V3f*                m_data; // device ptr
	int                         m_pixelSize;
};

