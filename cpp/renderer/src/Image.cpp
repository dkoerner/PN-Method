#include <Image.h>



#include <OpenEXR/ImfInputFile.h>
#include <OpenEXR/ImfOutputFile.h>
#include <OpenEXR/ImfChannelList.h>
#include <OpenEXR/ImfStringAttribute.h>
#include <OpenEXR/ImfVersion.h>
#include <OpenEXR/ImfHeader.h>
#include <OpenEXR/ImfIO.h>









Image::Image(const Vector2i &size)
	: m_data(size.y(), size.x())
{
	m_data.fill( Eigen::Vector3d(0.0, 0.0, 0.0) );
}

Image::Image(const std::string& filename)
{
	//
	// OpenEXR reading code has been commented to get rid of openexr dependency in this
	// stripped down version of nori
	//
	Imf::InputFile file(filename.c_str());
	const Imf::Header &header = file.header();
	const Imf::ChannelList &channels = header.channels();

	Imath::Box2i dw = file.header().dataWindow();
	V2i res(dw.max.x - dw.min.x + 1, dw.max.y - dw.min.y + 1);

	//m_data.resize(dw.max.y - dw.min.y + 1, dw.max.x - dw.min.x + 1);
	//std::cout << "Image::Image: Reading a " << m_data.cols() << "x" << m_data.rows() << " OpenEXR file from \"" << filename << "\"" << endl;

	const char *ch_r = NULL, *ch_g = NULL, *ch_b = NULL;
	for (Imf::ChannelList::ConstIterator it = channels.begin(); it != channels.end(); ++it) {
		std::string name = toLower(it.name());

		if (it.channel().xSampling != 1 || it.channel().ySampling != 1) {
			// Sub-sampled layers are not supported
			continue;
		}

		if (!ch_r && (name == "r" || name == "red" ||
				endsWith(name, ".r") || endsWith(name, ".red"))) {
			ch_r = it.name();
		} else if (!ch_g && (name == "g" || name == "green" ||
				endsWith(name, ".g") || endsWith(name, ".green"))) {
			ch_g = it.name();
		} else if (!ch_b && (name == "b" || name == "blue" ||
				endsWith(name, ".b") || endsWith(name, ".blue"))) {
			ch_b = it.name();
		}
	}

//	Eigen::Array<Eigen::Vector3f, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> img_data( res[1], res[0] );
//	if (!ch_r || !ch_g || !ch_b)
//		throw NoriException("This is not a standard RGB OpenEXR file!");

	size_t compStride = sizeof(float),
		   pixelStride = 3 * compStride,
		   rowStride = pixelStride * res[0];

	// image data container
	int size = res[0]*res[1]*3;
	std::vector<float> exr_data_float( size );


	char *ptr = reinterpret_cast<char *>(exr_data_float.data());

	Imf::FrameBuffer frameBuffer;
	frameBuffer.insert(ch_r, Imf::Slice(Imf::FLOAT, ptr, pixelStride, rowStride)); ptr += compStride;
	frameBuffer.insert(ch_g, Imf::Slice(Imf::FLOAT, ptr, pixelStride, rowStride)); ptr += compStride;
	frameBuffer.insert(ch_b, Imf::Slice(Imf::FLOAT, ptr, pixelStride, rowStride));
	file.setFrameBuffer(frameBuffer);
	file.readPixels(dw.min.y, dw.max.y);

	// convert to double
	Base m(res[1], res[0]);
	double* ptr_dest = (double*)m.data();
	for( int i=0;i<size;++i )
	{
		ptr_dest[i] = double(exr_data_float[i]);
	}

	//
	//m_data = m.colwise().reverse();
	m_data = m;
}

Image::Image(const Image& image)
	: m_data(image.getResolution()[1], image.getResolution()[0])
{
	//memcpy( m_data.data(), image.data(), sizeof(V3d)*m_data.rows()*m_data.cols() );
	memcpy( m_data.data(), image.data(), sizeof(V3d)*m_data.rows()*m_data.cols() );
}


V3d Image::eval( const P2d& pRaster )const
{
	///*
	P2d xy = pRaster-P2d(0.5, 0.5);

	// lower left corner
	P2i ij( int(floor(xy.x())),
			int(floor(xy.y())));

	// fractional part
	double tx = xy.x()-ij.x();
	double ty = xy.y()-ij.y();

	// upper right corner
	V2i ij2 = ij+V2i(1, 1);

	// clamp indexing coordinates
	ij[0] = std::max(0, std::min(ij.x(), int(m_data.cols()-1)));
	ij2[0] = std::max(0, std::min(ij2.x(), int(m_data.cols()-1)));
	ij[1] = std::max(0, std::min(ij.y(), int(m_data.rows()-1)));
	ij2[1] = std::max(0, std::min(ij2.y(), int(m_data.rows()-1)));

	V3d result(0.0, 0.0, 0.0);
	result += pixel(ij.x(), ij.y())*(1.0-tx)*(1.0-ty);
	result += pixel(ij2.x(), ij.y())*tx*(1.0-ty);
	result += pixel(ij.x(), ij2.y())*(1.0-tx)*ty;
	result += pixel(ij2.x(), ij2.y())*tx*ty;
	return result;
}



const Eigen::Vector3d* Image::data()const
{
	return m_data.data();
}

Eigen::Vector3d* Image::data()
{
	return m_data.data();
}

Image::Base& Image::getArray()
{
	return m_data;
}

Eigen::Vector3d& Image::pixel( int x, int y )
{
	return m_data.coeffRef(y, x);
}

const Eigen::Vector3d& Image::pixel( int x, int y )const
{
	return m_data.coeffRef(y, x);
}

void Image::save(const std::string& filename)
{
	std::cout << "saving image to " << filename << " resolution=" << getResolution().toString() << std::endl;

	V2i res = getResolution();
	Imf::Header header(res[0], res[1]);
	header.insert("comments", Imf::StringAttribute("generated by renderer"));

	Imf::ChannelList &channels = header.channels();
	channels.insert("R", Imf::Channel(Imf::FLOAT));
	channels.insert("G", Imf::Channel(Imf::FLOAT));
	channels.insert("B", Imf::Channel(Imf::FLOAT));

	Imf::FrameBuffer frameBuffer;
	size_t compStride = sizeof(float),
		   pixelStride = 3 * compStride,
		   rowStride = pixelStride * res[0];

	// exr coordinate system has its origin at the top-left corner of the image,
	// while we have the origin at the bottom left of the image
	// therefore we need to flip the image data vertically
	//Base img_flipped_vertically = m_data.colwise().reverse();
	Base img_flipped_vertically = m_data;

	// we convert to float
	int size = res[0]*res[1]*3;
	std::vector<float> converted( size );
	double* org = img_flipped_vertically.data()->data();
	for( int i=0;i<size;++i )
		converted[i] = float(org[i]);


	char *ptr = reinterpret_cast<char *>(converted.data());

	frameBuffer.insert("R", Imf::Slice(Imf::FLOAT, ptr, pixelStride, rowStride)); ptr += compStride;
	frameBuffer.insert("G", Imf::Slice(Imf::FLOAT, ptr, pixelStride, rowStride)); ptr += compStride;
	frameBuffer.insert("B", Imf::Slice(Imf::FLOAT, ptr, pixelStride, rowStride));

	Imf::OutputFile file(filename.c_str(), header);
	file.setFrameBuffer(frameBuffer);
	file.writePixels(res[1]);
}

V2i Image::getResolution()const
{
	return V2i(m_data.cols(), m_data.rows());
}

P2d Image::uvToRaster(const P2d& uv)const
{
	return P2d( uv[0]*(m_data.cols()-1),
				uv[1]*(m_data.rows()-1));
}

P2d Image::rasterToUV(const P2d& pRaster)const
{
	return P2d( pRaster[0]/(m_data.cols()-1),
				pRaster[1]/(m_data.rows()-1) );
}

std::string Image::toString()const
{
	std::ostringstream ss;
	ss << "Image resolution=" << getResolution().toString() << std::endl;
	return ss.str();
}

double Image::luminance( const V3d& c )
{
	return c[0] * 0.212671 + c[1] * 0.715160 + c[2] * 0.072169;
}

Image::Ptr blur_image( Image::Ptr image, double stddev )
{
	V2i res = image->getResolution();
	Image::Ptr image_blurred = std::make_shared<Image>( res );

	double variance = stddev*stddev;
	int filter_size = int(6*stddev - 1);
	auto weight = [&](int i, int j, int ki, int kj)
	{
		double x = i+0.5;
		double y = j+0.5;
		double kx = ki+0.5;
		double ky = kj+0.5;
		double dx = kx-x;
		double dy = ky-y;

		// we dont have a normalization constant since we normalize weights
		// over the pixelneighbourhood anyways...
		return std::exp( -(dx*dx + dy*dy)/(2*variance) );
	};

	int kshift = filter_size/2; // used to center the kernel at current pixel
	int resx = res[0];
	int resy = res[1];
	for( int j=0;j<resy;++j )
		for( int i=0;i<resx;++i )
		{
			Eigen::Vector3d c(0.0, 0.0, 0.0);
			// iterate over kernel
			int j_kstart = std::max( j-kshift, 0 );
			int j_kend = std::min( j+kshift, int(res[1])-1);
			int i_kstart = std::max( i-kshift, 0 );
			int i_kend = std::min(i+kshift, int(res[0])-1);
			double sum = 0.0;
			for( int kj=j_kstart;kj<=j_kend;++kj )
				for( int ki=i_kstart;ki<=i_kend;++ki )
				{
					double w = weight( i, j, ki, kj );
					c += image->pixel(ki, kj)*w;
					sum += w;
				}

			// set center pixel
			image_blurred->pixel(i, j) = c/sum;
		}

	return image_blurred;
}


// ImageSampler -------------------------------

P2d ImageSampler::sample(double &pdf, RNGd &rng)
{
	double xPdf;
	double yPdf;

	//sample xy pixel position based on the cdf of the envmaps radiance
	double r1 = rng.next1D();
	int x = m_colDpdf.sample(r1,xPdf);
	double r2 = rng.next1D();
	int y = m_rowDpdf[x].sample(r2,yPdf);

	pdf = xPdf*yPdf;

	return P2d(x, y) + P2d(rng.next1D(), rng.next1D());
}

P2d ImageSampler::sample(double &pdf, P2d sample)
{
	double xPdf;
	double yPdf;

	//sample xy pixel position based on the cdf of the envmaps radiance
	double r1 = sample[0];
	int x = m_colDpdf.sample(r1,xPdf);
	double r2 = sample[1];
	int y = m_rowDpdf[x].sample(r2,yPdf);

	pdf = xPdf*yPdf;

	return P2d(x, y);
}

P2d ImageSampler::rasterToUV(const P2d &xy) const
{
	return P2d( (xy.x())/double(m_resolution[0]-1),
				(xy.y())/double(m_resolution[1]-1));
}

void ImageSampler::buildPDF( Image::Ptr image )
{
	m_resolution = image->getResolution();
	int xres = m_resolution[0];
	int yres = m_resolution[1];

	//std::cout << "xres=" << xres << std::endl;
	//d::cout << "yres=" << yres << std::endl;

	m_rowDpdf.resize(xres);

	for( int i=0;i<xres;++i )
	{
		float sy = 0.0;
		for( int j=0;j<yres;++j )
		{
			P2d uv = rasterToUV(P2d(i+0.5, j+0.5));
			double theta = M_PI*uv.y();
			double c = Image::luminance(image->pixel(i,j))*std::sin(theta);
			sy += c;
			m_rowDpdf[i].append(c);
		}

		m_rowDpdf[i].normalize();
		m_colDpdf.append(sy);
	}
	m_colDpdf.normalize();
}
