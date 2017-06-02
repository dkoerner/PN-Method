#pragma once


#include <memory>
#include <math/vector.h>
#include <math/transform.h>

#include <Eigen/Geometry>

struct Camera
{

	Camera( int xres, int yres, double near, double far ) :
		m_resolution(xres, yres),
		m_invResolution(V2d(1.0/xres, 1.0/yres)),
		m_nearClip(near),
		m_farClip(far)
	{
		// world space transform of camera
		setCameraToWorld(Transformd());
	}

	typedef std::shared_ptr<Camera> Ptr;
	void setCameraToWorld( const Transformd& cameraToWorld )
	{
		m_cameraToWorld = cameraToWorld;
	}


	int getResolutionX()const
	{
		return m_resolution.x();
	}

	int getResolutionY()const
	{
		return m_resolution.y();
	}

	virtual void sampleRay( const Point2d& rasterP, Ray3d& ray, bool debug = false )const =0;


	V2i         m_resolution;
	V2d         m_invResolution;
	Transformd  m_cameraToWorld;  // camera space to world space
	Transformd  m_rasterToCamera; // normalized raster space to camera
	double      m_nearClip;
	double      m_farClip;
};

struct OrthographicCamera : public Camera
{

	OrthographicCamera( int xres=1, int yres=1, double width=2.0, double height=2.0, double near=0.0, double far=10.0):
		Camera(xres, yres, near, far)
	{
		float aspect = float(xres)/float(yres);

		// build transform from camera to raster space
		Eigen::Affine3d cameraToRaster = Eigen::Affine3d::Identity();

		// move negative near distance (bring points on the near plane into the xy plane)
		cameraToRaster = Eigen::Translation3d(V3d(0.0, 0.0, -near)) * cameraToRaster;
		// scale near-far range into 0-1 range
		cameraToRaster = Eigen::DiagonalMatrix<double, 3>(V3d(1.0, 1.0, far-near)) * cameraToRaster;
		// translate to bring xy from -1,1 range into 0,2 range
		cameraToRaster = Eigen::Translation3d(V3d(width*0.5, height*0.5*(1.0f/aspect), 0.0)) * cameraToRaster;
		// scale to bring xy from 0,2 range into 0,1 range (normalized raster space)
		cameraToRaster = Eigen::DiagonalMatrix<double, 3>(V3d(1.0/width, (1.0/height)*aspect, 1.0)) * cameraToRaster;
		// scale to bring xy samples from 01 range to xres/yres resolution range (raster space)
		cameraToRaster = Eigen::DiagonalMatrix<double, 3>(V3d(m_resolution.x(), m_resolution.y(), 1.0)) * cameraToRaster;


		// transforms from normalized raster space into camera space
		m_rasterToCamera = Transformd(cameraToRaster).inverse();
	}



	virtual void sampleRay( const Point2d& rasterP, Ray3d& ray, bool debug = false )const override
	{
		// transform from raster space to camera space
		P3d cameraNearPlaneP = m_rasterToCamera*P3d(rasterP.x(), rasterP.y(), 0.0);
		// transform camera space position to world space
		ray = Ray3d( m_cameraToWorld*cameraNearPlaneP,
					 m_cameraToWorld.getMatrix().topLeftCorner<3,3>()*V3d(0.0, 0.0, 1.0));

		if(debug)
		{
			std::cout << "OrthographicCamera::sampleRay: rasterP=" << rasterP.toString() << std::endl;
			std::cout << "OrthographicCamera::sampleRay: cameraNearPlaneP=" << cameraNearPlaneP.toString() << std::endl;
			std::cout << "OrthographicCamera::sampleRay: ray=" << ray.toString() << std::endl;
			std::cout << "OrthographicCamera::sampleRay: m_cameraToWorld=\n" << m_cameraToWorld.getMatrix() << std::endl;
		}

	}
};


struct PerspectiveCamera : public Camera
{
	PerspectiveCamera( int xres=1, int yres=1, double hfov=30.0, double near=1e-4f, double far=1000.0):
		Camera(xres, yres, near, far),
		m_hfov(hfov)
	{
		double aspect = double(xres)/double(yres);
		m_apertureRadius = 0.0;

		// Project vectors in camera space onto a plane at z=1:
		//
		//  xProj = cot * x / z
		//  yProj = cot * y / z
		//  zProj = (far * (z - near)) / (z * (far-near))
		//  The cotangent factor ensures that the field of view is
		//  mapped to the interval [-1, 1].
		//
		double recip = 1.0f / (m_farClip - m_nearClip),
			  cot = 1.0f / std::tan(degToRad(m_hfov / 2.0f));
		Eigen::Matrix4d perspective;
		perspective <<
			cot, 0,   0,   0,
			0, -cot,   0,   0, // the minus sign is because image plane is flipped vertically
			0,   0,   m_farClip * recip, -m_nearClip * m_farClip * recip,
			0,   0,   1,   0;
		m_rasterToCamera = Transformd( Eigen::DiagonalMatrix<double, 3>(Vector3d(-0.5, -0.5 * aspect, 1.0))*
							Eigen::Translation<double, 3>(-1.0, -1.0/aspect, 0.0)*perspective ).inverse();
	}



	virtual void sampleRay( const Point2d& rasterP, Ray3d& ray, bool debug = false )const override
	{
		P2d tmp(0.0, 0.0);
		double m_focusDistance = m_farClip;


		/* Compute the corresponding position on the
		   near plane (in local camera space) */
		P3d nearP = m_rasterToCamera * P3d(
			rasterP.x() * m_invResolution.x(),
			rasterP.y() * m_invResolution.y(), 0.0f);

		P3d apertureP(tmp.x(), tmp.y(), 0.0f);

		/* Sampled position on the focal plane */
		P3d focusP = nearP * (m_focusDistance / nearP.z());

		/* Aperture position */
		/* Turn these into a normalized ray direction, and
		   adjust the ray interval accordingly */
		V3d d = (focusP - apertureP).normalized();
		float invZ = 1.0f / d.z();

		ray.o = m_cameraToWorld * apertureP;
		ray.d = m_cameraToWorld * d;
		ray.mint = m_nearClip * invZ;
		ray.maxt = m_farClip * invZ;
		ray.update();
	}


	double m_hfov;
	double m_apertureRadius;

};



Camera::Ptr createView(const V3d& viewDir, const Box3d& bound, V2i res, double hfov );
Camera::Ptr createView( const V3d& viewDir, const Box3d& bound, V2i res );
Camera::Ptr createView( const V3d& viewDir, const Box3d& bound, V2i res, const Box2d& region );
