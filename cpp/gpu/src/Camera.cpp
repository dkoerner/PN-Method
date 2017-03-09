#include <camera.h>

#include <math/bbox.h>









Camera::Ptr createView( const V3d& viewDir, const Box3d& bound, V2i res, double hfov )
{
	P3d targetWS = bound.getCenter();
	double aspect = double(res.x())/double(res.y());
	double fov[2] = {hfov*aspect, hfov};

	// computed 2d bounding box of volume in image plane
	Transformd cameraToWorld = Eigen::Affine3d(lookAt<double>(targetWS-viewDir, targetWS));

	BoundingBox3d boundCS;
	boundCS.reset();
	for( int i=0;i<8;++i )
	{
		// project point into plane given by direction normal vector
		P3d pWS = bound.getCorner(i);
		P3d pCS = cameraToWorld.inverse()*pWS;

		// we get rid of the z dimension
		boundCS.expandBy( P3d(pCS.x(), pCS.y(), 0.0) );
	}

	// dimension of maximum extend
	int axis = boundCS.getMajorAxis();

	// now we know the image plane extend of our object
	// find the distance along z at which the whole object will be contained
	// in our frustum
	double a = boundCS.getExtents()[axis]*0.5*1.5;
	double theta = degToRad(fov[axis]*0.5);
	double b = a/std::tan(theta);

	// TODO: compute near and far values

	// setup camera
	Camera::Ptr camera = std::make_shared<PerspectiveCamera>(res.x(), res.y(), hfov, 1e-4f, 5000.0);
	camera->setCameraToWorld( Eigen::Affine3d(lookAt<double>(targetWS-viewDir*b, targetWS)) );
	return camera;
}


Camera::Ptr createView( const V3d& viewDir, const Box3d& bound, V2i res )
{
	P3d targetWS = bound.getCenter();
	double aspect = double(res.x())/double(res.y());
	double width[2] = {aspect, 1.0};

	// computed 2d bounding box of volume in image plane
	Transformd cameraToWorld = Eigen::Affine3d(lookAt<double>(targetWS-viewDir, targetWS));

	BoundingBox3d boundCS;
	boundCS.reset();
	for( int i=0;i<8;++i )
	{
		// project point into plane given by direction normal vector
		P3d pWS = bound.getCorner(i);
		P3d pCS = cameraToWorld.inverse()*pWS;

		boundCS.expandBy( P3d(pCS.x(), pCS.y(), pCS.z()) );
	}
	//std::cout << "bound=" << bound.toString() << std::endl;
	//std::cout << "boundCS=" << boundCS.toString() << std::endl;

	// dimension of maximum extend
	double b = boundCS.getExtents().z();
	//std::cout << "b=" << b << std::endl;

	// TODO: compute near and far values

	// setup camera
	//Camera::Ptr camera = std::make_shared<OrthographicCamera>(res.x(), res.y(), 1.5, 1.5, 1e-4f, 5000.0);
	Camera::Ptr camera = std::make_shared<OrthographicCamera>(res.x(), res.y(), boundCS.getExtents().x()*1.3, boundCS.getExtents().y()*1.3, 1e-4f, 5000.0);
	camera->setCameraToWorld( Eigen::Affine3d(lookAt<double>(targetWS-viewDir*b, targetWS)) );
	return camera;
}


Camera::Ptr createView( const V3d& viewDir, const Box3d& bound, V2i res, const Box2d& region )
{
	P3d targetWS = bound.getCenter();
	double aspect = double(res.x())/double(res.y());
	double width[2] = {aspect, 1.0};

	// computed 2d bounding box of volume in image plane
	Transformd cameraToWorld = Eigen::Affine3d(lookAt<double>(targetWS-viewDir, targetWS));

	BoundingBox3d boundCS;
	boundCS.reset();
	for( int i=0;i<8;++i )
	{
		// project point into plane given by direction normal vector
		P3d pWS = bound.getCorner(i);
		P3d pCS = cameraToWorld.inverse()*pWS;

		boundCS.expandBy( P3d(pCS.x(), pCS.y(), pCS.z()) );
	}
	//std::cout << "bound=" << bound.toString() << std::endl;
	//std::cout << "boundCS=" << boundCS.toString() << std::endl;

	// dimension of maximum extend
	double b = boundCS.getExtents().z();
	//std::cout << "b=" << b << std::endl;

	// TODO: compute near and far values

	// setup camera
	//Camera::Ptr camera = std::make_shared<OrthographicCamera>(res.x(), res.y(), 1.5, 1.5, 1e-4f, 5000.0);
	Camera::Ptr camera = std::make_shared<OrthographicCamera>(res.x(), res.y(), boundCS.getExtents().x()*1.3, boundCS.getExtents().y()*1.3, 1e-4f, 5000.0);
	camera->setCameraToWorld( Eigen::Affine3d(lookAt<double>(targetWS-viewDir*b, targetWS)) );
	return camera;
}
