#pragma once


#include <memory>
#include <math/vector.h>
#include <math/transform.h>

#include <Eigen/Geometry>

struct Camera
{
	typedef std::shared_ptr<Camera> Ptr;

	Camera( int xres, int yres, double near, double far ) :
		m_resolution(xres, yres),
		m_invResolution(V2d(1.0/xres, 1.0/yres)),
		m_nearClip(near),
		m_farClip(far)
	{
		// world space transform of camera
		setCameraToWorld(Transformd());
	}

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

	V2i getResolution()const
	{
		return m_resolution;
	}

	virtual std::string toString()const=0;

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

	virtual std::string toString()const override
	{
		std::ostringstream ss;
		ss << "OrthographicCamera " << std::endl;
		return ss.str();
	}
};


struct PerspectiveCamera : public Camera
{
	typedef std::shared_ptr<PerspectiveCamera> Ptr;

	PerspectiveCamera( int xres=1, int yres=1, double hfov=30.0, double near=1e-4f, double far=1000.0):
		Camera(xres, yres, near, far),
		m_hfov(hfov),
		m_z_direction(-1.0)
	{
		double aspect = double(xres)/double(yres);
		m_apertureRadius = 0.0;

		// this matrix transforms a homogeneous point p=[x, y, z, 1] from camera space
		// into clip space, such that points in range [nearClip, farClip] are mapped to [-1, 1]
		// (after normalization of homogeneous coordinates)
		// in eigen this is done by: P3d p_clip = (m_clipToCamera * p_camera.colwise().homogeneous()).colwise().hnormalized();
		// NB: that in camera space near->far is in -z direction while in clipspace near->far is in +z direction
		// I found this video helpfull: https://www.youtube.com/watch?v=dul0mui292Q ()
		double cot =  1.0f / std::tan(degToRad(m_hfov / 2.0f));
		Eigen::Matrix4d cameraToClip;
		cameraToClip <<
			cot, 0,   0,   0,
			0, cot,   0,   0,
			0,   0,   m_z_direction*-(m_nearClip+m_farClip)/(m_nearClip-m_farClip), m_z_direction*-2.0*m_nearClip*m_farClip/(m_nearClip-m_farClip),
			0,   0,   m_z_direction*1,   0;


		Transformd cameraToRaster = Transformd(
					Eigen::DiagonalMatrix<double, 3>(Vector3d(xres*0.5, yres*0.5 * aspect, 1.0))* // 3. scale clipspace, such that x and y coordinates range [0, 1] and further [0, xres/yres]
					Eigen::Translation<double, 3>(1.0, 1.0/aspect, 0.0)*                          // 2. move clipspace (where x and y coordinates go from [-1, 1]), such that x and y coordinates range [0, 2]
					cameraToClip                                                                // 1. project from camera space into clip-space (see above)
					);

		m_rasterToCamera = cameraToRaster.inverse();
	}



	virtual void sampleRay( const Point2d& rasterP, Ray3d& ray, bool debug = false )const override
	{
		P2d tmp(0.0, 0.0);
		double m_focusDistance = m_farClip;


		/* Compute the corresponding position on the
		   near plane (in local camera space) */
		P3d nearP = m_rasterToCamera * P3d( rasterP.x(), rasterP.y(), 0.0f);


		P3d apertureP(tmp.x(), tmp.y(), 0.0f);

		/* Sampled position on the focal plane */
		P3d focusP = nearP * (m_z_direction*m_focusDistance / nearP.z());

		/* Aperture position */
		/* Turn these into a normalized ray direction, and
		   adjust the ray interval accordingly */
		V3d d = (focusP - apertureP).normalized();
		float invZ = 1.0f / d.z();

		ray.o = m_cameraToWorld * apertureP;
		ray.d = m_cameraToWorld * d;
		ray.mint = m_z_direction*m_nearClip * invZ;
		ray.maxt = m_z_direction*m_farClip * invZ;
		ray.update();
	}

	virtual std::string toString()const override
	{
		std::ostringstream ss;
		ss << "PerspectiveCamera res=" << m_resolution.toString() << " hfov=" <<m_hfov << std::endl;
		return ss.str();
	}

	double m_hfov;
	double m_apertureRadius;
	double m_z_direction; // this is -1 if camera points into negative z, or 1 if camera points into positive z direction

};


/*
Camera::Ptr createView(const V3d& viewDir, const Box3d& bound, V2i res, double hfov );
Camera::Ptr createView( const V3d& viewDir, const Box3d& bound, V2i res );
Camera::Ptr createView( const V3d& viewDir, const Box3d& bound, V2i res, const Box2d& region );
*/
