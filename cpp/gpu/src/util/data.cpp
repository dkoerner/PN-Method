#include <util/data.h>
#include <houio/HouGeoIO.h>
#include <houio/Geometry.h>


void initializeData( const Slice& slice, int resX, int resY, Data& data )
{
	data.clear();
	data.m_type = Data::Type::EImage;
	data.m_xres = resX;
	data.m_yres = resY;

//	for( int j=resY/2;j<resY;++j )
//		for( int i=resX/2;i<resX;++i )
//		{
//			int x = i;
//			int y = j;
//			V3d o, d;
//			camera.sampleRay(x+0.5, y+0.5, o, d);
//			data.addSensor(o);
//		}


	double dx = slice.width/resX;
	double dy = slice.height/resY;
	for( int j=0;j<resY;++j )
		for( int i=0;i<resX;++i )
		{
			V3d o = slice.o + slice.right*(i*dx + dx*0.5) + slice.up*(j*dy + dy*0.5);
			data.addSensor(o);
		}

}


void writeImage( const std::vector<double>& samples, int resX, int resY, const std::string& filename )
{
	std::shared_ptr<Bitmap> bitmap = toImage(samples, resX, resY);
	bitmap->saveEXR(filename);
	//bitmap->saveTXT(filename+".txt");
}

std::shared_ptr<Bitmap> toImage( const std::vector<double> &samples, int resX, int resY)
{
	std::shared_ptr<Bitmap> bitmap = std::make_shared<Bitmap>(Vector2i(resX, resY));

	int index =0;
	for(int y=0;y<resY;++y)
		for(int x=0;x<resX;++x, ++index)
		{
			double val = samples[index];
			bitmap->coeffRef(y, x) = Color3f(val, val, val);
		}

	return bitmap;
}




/*






void natlogSamples( double minRadius, double maxRadius, int num, std::vector<double>& samples )
{
	linearSamples( minRadius, maxRadius, num, samples );
	for( auto& s:samples )
		s = std::exp(s);
}

void saveRadialData( const Data& data, const std::string& filename, const P3d& refP  )
{
	std::vector<double> radii;
	for(auto sensor : data.sensors)
		radii.push_back((refP-sensor).norm());
	writeSamples( filename, radii, data.estimates );
}















void initializeData(Camera::Ptr camera, Data &data)
{
	int resX = camera->getResolutionX();
	int resY = camera->getResolutionY();
	data.clear();

// // quarter image
//	for( int j=resY/2;j<resY;++j )
//		for( int i=resX/2;i<resX;++i )
//		{
//			int x = i;
//			int y = j;
//			Ray3d ray;
//			camera->sampleRay(P2d(x+0.5, y+0.5), ray);
//			data.addRay(ray);
//		}


	for( int j=0;j<resY;++j )
		for( int i=0;i<resX;++i )
		{
			int x = i;
			int y = j;
			Ray3d ray;
			camera->sampleRay(P2d(x+0.5, y+0.5), ray);
			data.addRay(ray);
		}
}

*/


// type dependent save function
void Data::save(const std::string &filename) const
{
	switch(m_type)
	{
		default:
		case Type::EUnknown:
		{
			std::cout << "Data::save: unknown type of sensor data\n";
			if( estimates.size() == 1 )
				std::cout << "Data::save: single estimate is: " << estimates[0] << std::endl;
		}break;
		case Type::ERadial:
		{
			saveRadialData( sensors, estimates, filename, m_o );
		}break;
		case Type::EImage:
		{
			writeImage( estimates, m_xres, m_yres, filename + ".exr" );
		}break;
	};
}











void xport_observations(const std::string &filename, std::vector<std::pair<P3d, double>>& obervations)
{
	houio::Geometry::Ptr geo = houio::Geometry::createPointGeometry();

	houio::Attribute::Ptr pAttr = geo->getAttr("P");
	houio::Attribute::Ptr valueAttr = houio::Attribute::createFloat();
	houio::Attribute::Ptr cAttr = houio::Attribute::createV3f(); // color (for houdini)
	//houio::Attribute::Ptr albedoAttr = houio::Attribute::createV3f();
	//houio::Attribute::Ptr extinctionAttr = houio::Attribute::createV3f();
	//houio::Attribute::Ptr extinctionAttr = houio::Attribute::createV3f();


	int i = 0;
	for( auto& o:obervations )
	{
		//std::cout << "observation " << i++ << " " << o.first.x() << " " << o.first.y() << " " << o.first.z() << " " << o.second << std::endl;
		pAttr->appendElement( float(o.first.x()), float(o.first.y()), float(o.first.z()) );
		valueAttr->appendElement( float(o.second) );
		cAttr->appendElement( float(o.second), float(o.second), float(o.second) );
	}
	geo->setAttr("value", valueAttr);
	geo->setAttr("Cd", cAttr);

	std::cout << "xport_observations: writing " << filename << std::endl;
	houio::HouGeoIO::xport(filename, geo);
}

void xport_line(const std::string &filename, P3d& a, P3d& b)
{
	houio::Geometry::Ptr geo = houio::Geometry::createLineGeometry();

	houio::Attribute::Ptr pAttr = geo->getAttr("P");
	pAttr->appendElement( float(a.x()), float(a.y()), float(a.z()) );
	pAttr->appendElement( float(b.x()), float(b.y()), float(b.z()) );

	geo->addLine(0, 1);

	std::cout << "xport_line: writing " << filename << std::endl;
	houio::HouGeoIO::xport(filename, geo);
}
