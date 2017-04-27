#include <util/envmap.h>
#include <houio/Geometry.h>
#include <houio/HouGeoIO.h>




void EnvMap::saveGeo( const std::string& filename, double exposure )
{
	double scale = std::pow(2.0,exposure );
	houio::Geometry::Ptr geo = houio::Geometry::createSphere(120, 120, 1.0);
	houio::Attribute::Ptr pAttr = geo->getAttr("P");
	houio::Attribute::Ptr cdAttr = houio::Attribute::createV3f(pAttr->numElements());
	for( int i=0;i<pAttr->numElements();++i )
	{
		houio::math::V3f p = pAttr->get<houio::math::V3f>(i);
		P2d theta_phi = sphericalCoordinates<double>(V3d(p.x, p.y, p.z));
		double theta = theta_phi.x();
		double phi = theta_phi.y();
		Color3f col = this->eval(theta, phi)*scale;
		cdAttr->set<houio::math::V3f>( i, houio::math::V3f(col.r(), col.g(), col.b()) );
	}
	geo->setAttr("Cd", cdAttr);
	houio::HouGeoIO::xport( filename, geo);
}


void rasterizeSphericalFunctionSphere(const std::string& filename, std::function<Color3f (double, double)> func, double exposure, int res )
{
	double scale = std::pow(2.0, exposure);
	houio::Geometry::Ptr geo = houio::Geometry::createSphere(res, res, 1.0);
	houio::Attribute::Ptr pAttr = geo->getAttr("P");
	houio::Attribute::Ptr cdAttr = houio::Attribute::createV3f(pAttr->numElements());
	int numPoints = pAttr->numElements();
	for( int i=0;i<numPoints;++i )
	{
		std::cout << i << "/" <<numPoints << std::endl;
		houio::math::V3f& p = pAttr->get<houio::math::V3f>(i);
		P2d theta_phi = sphericalCoordinates<double>(V3d(p.x, p.y, p.z));
		double theta = theta_phi.x();
		double phi = theta_phi.y();
		Color3f col = func(theta, phi)*scale;

		/*
		p.x = col.r();
		p.y = col.g();
		p.z = col.b();
		*/

		//pAttr->set<houio::math::V3f>( i, p*float(col.r()) );
		cdAttr->set<houio::math::V3f>( i, houio::math::V3f(col.r(), col.g(), col.b()) );
	}
	geo->setAttr("Cd", cdAttr);
	houio::HouGeoIO::xport( filename, geo);
}



void rasterizeSphericalFunctionSphere(const std::string& filename, std::function<Color3f (double, double)> func, Transformd transform, double exposure, int res = 120)
{
	double scale = std::pow(2.0, exposure);
	houio::Geometry::Ptr geo = houio::Geometry::createSphere(res, res, 1.0);
	houio::Attribute::Ptr pAttr = geo->getAttr("P");
	houio::Attribute::Ptr cdAttr = houio::Attribute::createV3f(pAttr->numElements());
	int numPoints = pAttr->numElements();
	for( int i=0;i<numPoints;++i )
	{
		std::cout << i << "/" <<numPoints << std::endl;
		houio::math::V3f& p = pAttr->get<houio::math::V3f>(i);
		V3d d = V3d(p.x, p.y, p.z);
		P2d theta_phi = sphericalCoordinates<double>(d);
		double theta = theta_phi.x();
		double phi = theta_phi.y();
		Color3f col = func(theta, phi)*scale;

		P3d transformed = transform*d;
		pAttr->set<houio::math::V3f>( i, houio::math::V3f(transformed.x(), transformed.y(), transformed.z()) );
		cdAttr->set<houio::math::V3f>( i, houio::math::V3f(col.r(), col.g(), col.b()) );
	}
	geo->setAttr("Cd", cdAttr);
	houio::HouGeoIO::xport( filename, geo);
}

void displaceSphere(const std::string& filename, std::function<double (double, double)> func, double exposure)
{
	double scale = std::pow(2.0, exposure);
	houio::Geometry::Ptr geo = houio::Geometry::createSphere(120, 120, 1.0);
	houio::Attribute::Ptr pAttr = geo->getAttr("P");
	houio::Attribute::Ptr cdAttr = houio::Attribute::createV3f(pAttr->numElements());
	int numPoints = pAttr->numElements();
	for( int i=0;i<numPoints;++i )
	{
		std::cout << i << "/" <<numPoints << std::endl;
		houio::math::V3f& p = pAttr->get<houio::math::V3f>(i);
		P2d theta_phi = sphericalCoordinates<double>(V3d(p.x, p.y, p.z));
		double theta = theta_phi.x();
		double phi = theta_phi.y();
		double value = func(theta, phi)*scale;

		pAttr->set<houio::math::V3f>( i, p*float(value) );
		//cdAttr->set<houio::math::V3f>( i, houio::math::V3f(col.r(), col.g(), col.b()) );
	}
	//geo->setAttr("Cd", cdAttr);
	houio::HouGeoIO::xport( filename, geo);
}

