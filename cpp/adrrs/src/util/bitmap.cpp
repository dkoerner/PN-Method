#include <util/bitmap.h>
#include <fstream>

#include <OpenEXR/ImfInputFile.h>
#include <OpenEXR/ImfOutputFile.h>
#include <OpenEXR/ImfChannelList.h>
#include <OpenEXR/ImfStringAttribute.h>
#include <OpenEXR/ImfVersion.h>
#include <OpenEXR/ImfHeader.h>
#include <OpenEXR/ImfIO.h>

#include <util/string.h>

Bitmap::Bitmap(const std::string &filename)
{
//	if (!QFile(filename).exists())
//		throw NoriException(QString("EXR file \"%1\" does not exist!").arg(filename));

	//
	// OpenEXR reading code has been commented to get rid of openexr dependency in this
	// stripped down version of nori
	//
	Imf::InputFile file(filename.c_str());
	const Imf::Header &header = file.header();
	const Imf::ChannelList &channels = header.channels();

	Imath::Box2i dw = file.header().dataWindow();
	resize(dw.max.y - dw.min.y + 1, dw.max.x - dw.min.x + 1);

	cout << "Reading a " << cols() << "x" << rows() << " OpenEXR file from \"" << filename << "\"" << endl;
	
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

//	if (!ch_r || !ch_g || !ch_b)
//		throw NoriException("This is not a standard RGB OpenEXR file!");

	size_t compStride = sizeof(float),
	       pixelStride = 3 * compStride,
	       rowStride = pixelStride * cols();

	char *ptr = reinterpret_cast<char *>(data());

	Imf::FrameBuffer frameBuffer;
	frameBuffer.insert(ch_r, Imf::Slice(Imf::FLOAT, ptr, pixelStride, rowStride)); ptr += compStride;
	frameBuffer.insert(ch_g, Imf::Slice(Imf::FLOAT, ptr, pixelStride, rowStride)); ptr += compStride;
	frameBuffer.insert(ch_b, Imf::Slice(Imf::FLOAT, ptr, pixelStride, rowStride)); 
	file.setFrameBuffer(frameBuffer);
	file.readPixels(dw.min.y, dw.max.y);

}



Color3f Bitmap::eval( const P2d& uv )const
{
	P2d xy(
		uv.x()*(cols()-1)-0.5,
		uv.y()*(rows()-1)-0.5
		);

	// lower left corner
	P2i ij( int(floor(xy.x())),
			int(floor(xy.y())));

	// fractional part
	P2d st( xy.x()-ij.x(),
			xy.y()-ij.y());

	// upper right corner
	V2i ij2 = ij+V2i(1, 1);

	// clamp indexing coordinates
	ij[0] = std::max(0, std::min(ij.x(), int(cols()-1)));
	ij2[0] = std::max(0, std::min(ij2.x(), int(cols()-1)));
	ij[1] = std::max(0, std::min(ij.y(), int(rows()-1)));
	ij2[1] = std::max(0, std::min(ij2.y(), int(rows()-1)));


	Color3f c = lerp(
				lerp(coeff(ij.y(), ij.x()), coeff(ij.y(), ij2.x()), st.x()),
				lerp(coeff(ij2.y(), ij.x()), coeff(ij2.y(), ij2.x()), st.x()),
				st.y());
	return c;

	// lower left
	//return coeff(ij.y(), ij.x());
}


/*
Color3f Bitmap::eval( const P2d& uv )const
{
	P2d xy(
		uv.x()*cols()-0.5,
		uv.y()*rows()-0.5
		);

	// lower left corner
	P2i ij( int(floor(xy.x())),
			int(floor(xy.y())));

	// fractional part
	P2d st( xy.x()-floor(xy.x()),
			xy.y()-floor(xy.y()));

	// upper right corner
	V2i ij2 = ij+V2i(1, 1);

	// clamp indexing coordinates
	ij[0] = std::max(0, std::min(ij.x(), int(cols()-1)));
	ij2[0] = std::max(0, std::min(ij2.x(), int(cols()-1)));
	ij[1] = std::max(0, std::min(ij.y(), int(rows()-1)));
	ij2[1] = std::max(0, std::min(ij2.y(), int(rows()-1)));

	Color3f c = lerp(
				lerp(coeff(ij.y(), ij.x()), coeff(ij.y(), ij2.x()), st.x()),
				lerp(coeff(ij2.y(), ij.x()), coeff(ij2.y(), ij2.x()), st.x()),
				st.y());
	return c;
}
*/
void Bitmap::saveEXR(const std::string &filename)
{
	std::cout << "Writing a " << cols() << "x" << rows() << " OpenEXR file to \"" << filename << "\"" << std::endl;

	Imf::Header header(cols(), rows());
	header.insert("comments", Imf::StringAttribute("generated by svpt"));

	Imf::ChannelList &channels = header.channels();
	channels.insert("R", Imf::Channel(Imf::FLOAT));
	channels.insert("G", Imf::Channel(Imf::FLOAT));
	channels.insert("B", Imf::Channel(Imf::FLOAT));

	Imf::FrameBuffer frameBuffer;
	size_t compStride = sizeof(float),
	       pixelStride = 3 * compStride,
	       rowStride = pixelStride * cols();

	char *ptr = reinterpret_cast<char *>(data());

	frameBuffer.insert("R", Imf::Slice(Imf::FLOAT, ptr, pixelStride, rowStride)); ptr += compStride;
	frameBuffer.insert("G", Imf::Slice(Imf::FLOAT, ptr, pixelStride, rowStride)); ptr += compStride;
	frameBuffer.insert("B", Imf::Slice(Imf::FLOAT, ptr, pixelStride, rowStride)); 

	Imf::OutputFile file(filename.c_str(), header);
	file.setFrameBuffer(frameBuffer);
	file.writePixels(rows());
}

void Bitmap::makeAbs()
{
	int numPixels = cols()*rows();
	for( int i=0;i<numPixels;++i )
	{
		Color3f c = coeffRef(i);
		c[0] = std::abs(c[0]);
		c[1] = std::abs(c[1]);
		c[2] = std::abs(c[2]);
		coeffRef(i) = c;
	}
}

void Bitmap::saveTXT( const std::string& filename )const
{
	std::cout << "writing pixelbuffer to: " << filename <<std::endl;
	std::ofstream out_values( filename.c_str(), std::ios_base::out | std::ios_base::trunc );

	int resX = cols();
	int resY = rows();

	out_values << resX << " " << resY << std::endl;
	for(int y=0;y<resY;++y)
		for(int x=0;x<resX;++x)
			out_values << coeffRef(y, x).x() << std::endl;
}


// computes mean square error
double Bitmap::computeMSE( const Bitmap& groundtruth )const
{
	int xres = cols();
	int yres = rows();

	if( (xres != groundtruth.cols())||
		(yres != groundtruth.rows()))
		throw std::runtime_error("Bitmap::computeMSE image dimensions need to match");

	const Color3f *aptr = data();
	const Color3f *bptr = groundtruth.data();

	double mean = 0.0;
	int count =0;
	for( int j=0;j<yres;++j )
		for( int i=0;i<xres;++i )
		{
			int index = j*xres +i;
			Color3f aa = aptr[index];
			Color3f bb = bptr[index];
			double error = aa[0]-bb[0];
			error*=error; // squared error

			mean += (error-mean)/(count+1);
			++count;
		}
	return mean;
}


void Bitmap::resizeInterpolated( int xres, int yres )
{
	// since we resize in place, we need to temporary hold the new pixels
	std::vector<Color3f> pixels(xres*yres);
	int c = 0;
	for( int j=0;j<yres;++j )
		for( int i=0;i<xres;++i,++c )
		{
			P2d uv( double(i)/double(xres-1), double(j)/double(yres-1) );
			pixels[c] = eval(uv);
		}

	// now resize array
	resize(yres, xres);

	// copy data
	memcpy( data(), &pixels[0], xres*yres*sizeof(Color3f) );
}

void Bitmap::resizeInterpolated( double scale )
{
	int newxres = int( scale*cols() );
	int newyres = int( scale*rows() );
	resizeInterpolated( newxres, newyres );
}



void compositeImages( std::vector<Bitmap::Ptr>& images, const std::string& filename, int maxCols, int spacing )
{
	if(images.empty())
		return;
	// now composite all images together
	{
		int numImages = int(images.size());
		int cols = numImages<maxCols ? numImages : maxCols;
		div_t divResult = div(numImages, cols);
		int rows = std::max(1, divResult.rem != 0? divResult.quot+1 : divResult.quot);
		int resx = images[0]->cols();
		int resy = images[0]->rows();

		Bitmap all( Vector2i(resx*cols+(cols-1)*spacing, resy*rows+(rows-1)*spacing) );
		all.fill(-1.0);
		int index = 0;

		for( int row = 0;row<rows;++row )
			for( int col = 0;col<cols;++col, ++index )
			{
				if(index < numImages)
				{
					Bitmap::Ptr img = images[index];
					all.block(row*resy+row*spacing,col*resx+col*spacing, resy, resx) = img->block(0,0,resy,resx);
				}else
					all.block(row*resy+row*spacing,col*resx+col*spacing, resy, resx).fill(0);
			}

		//rasterizeDot( all, V2d(0.05, 0.5), 0.04, Color3f(0.34182677047, 0.194753921622, 0.107302837237 ) );

		all.saveEXR(filename);
	}
}

Bitmap::Ptr compositeImages2( std::vector<Bitmap::Ptr>& images, int maxCols, int spacing )
{
	if(images.empty())
		return Bitmap::Ptr();
	// now composite all images together
	{
		int numImages = int(images.size());
		int cols = numImages<maxCols ? numImages : maxCols;
		div_t divResult = div(numImages, cols);
		int rows = std::max(1, divResult.rem != 0? divResult.quot+1 : divResult.quot);
		int resx = images[0]->cols();
		int resy = images[0]->rows();

		Bitmap::Ptr all = std::make_shared<Bitmap>( Vector2i(resx*cols+(cols-1)*spacing, resy*rows+(rows-1)*spacing) );
		all->fill(-1.0);
		int index = 0;

		for( int row = 0;row<rows;++row )
			for( int col = 0;col<cols;++col, ++index )
			{
				if(index < numImages)
				{
					Bitmap::Ptr img = images[index];
					all->block(row*resy+row*spacing,col*resx+col*spacing, resy, resx) = img->block(0,0,resy,resx);
				}else
					all->block(row*resy+row*spacing,col*resx+col*spacing, resy, resx).fill(0);
			}

		return all;
	}
}



void rasterizeDot( Bitmap& map, const V2d& center, double radius, Color3f value )
{
	int xres = map.cols();
	int yres = map.rows();

	for( int j=0;j<yres;++j )
		for( int i=0;i<xres;++i )
		{
			V2d uv( (i+0.5)/double(xres-1), (j+0.5)/double(yres-1) );
			if( (uv-center).norm() < radius )
				map.coeffRef( j, i ) = value;
		}
}


void expose( Bitmap& map, float exposure )
{
	int xres = map.cols();
	int yres = map.rows();

	for( int j=0;j<yres;++j )
		for( int i=0;i<xres;++i )
		{
			map.coeffRef( j, i ) *= std::pow(2.0, exposure);
		}
}

void putBitmap( Bitmap& dest, Bitmap& src, const V2d& dest_uv_min, const V2d& dest_uv_max )
{
	// find start and end indices in dest image
	P2i dest_xy_min = dest.uvToRaster( dest_uv_min );
	P2i dest_xy_max = dest.uvToRaster( dest_uv_max );

	//std::cout << "min: " << dest_xy_min.x() << " " << dest_xy_min.y() << std::endl;
	//std::cout << "max: " << dest_xy_max.x() << " " << dest_xy_max.y() << std::endl;

	for( int i=dest_xy_min.x();i<dest_xy_max.x();++i )
	{
		for( int j=dest_xy_min.y();j<dest_xy_max.y();++j )
		{
			if( (i>=0)&&(i<dest.cols())&&
				(j>=0)&&(j<dest.rows()) )
			{
				// work out uv coordinates in src image
				V2d uv( double(i-dest_xy_min.x())/double(dest_xy_max.x()-dest_xy_min.x()),
						double(j-dest_xy_min.y())/double(dest_xy_max.y()-dest_xy_min.y()));
				Color3f color = src.eval(uv);
				dest.coeffRef( j, i ) = color;
				//dest.coeffRef( j, i ) = Color3f(1.0, 0.0, 0.0);
			}
		}
	}
}


Bitmap::Ptr makeLeftRightImage(Bitmap& left, Bitmap& right)
{
	Bitmap::Ptr bitmap_combined = std::make_shared<Bitmap>(V2i( left.cols(), left.rows() ));
	for( int j=0;j<left.rows();++j )
		for( int i=0;i<left.cols();++i )
		{
			if( i<left.cols()/2 )
				bitmap_combined->coeffRef(j, i) = left.coeff(j, i);
			else
				bitmap_combined->coeffRef(j, i) = right.coeff(j, i);
		}
	return bitmap_combined;
}

void flip(Bitmap &map)
{
	Eigen::Array<Color3f, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> A = map.colwise().reverse();
	for( int j=0;j<map.cols();++j )
		for( int i=0;i<map.rows();++i )
			map.coeffRef(i, j) = A.coeffRef(i, j);
}

Color3f variance(Bitmap &map)
{
	std::vector<double> samples_r;
	std::vector<double> samples_g;
	std::vector<double> samples_b;
	int numPixels = map.cols()*map.rows();
	for( int i=0;i<numPixels;++i )
	{
		Color3f& c = map.coeffRef(i);
		samples_r.push_back(c[0]);
		samples_g.push_back(c[1]);
		samples_b.push_back(c[2]);
	}

	return Color3f(variance(samples_r), variance(samples_g), variance(samples_b));
}

Color3f variance(Bitmap &map, Color3f mean)
{
	std::vector<double> samples_r;
	std::vector<double> samples_g;
	std::vector<double> samples_b;
	int numPixels = map.cols()*map.rows();
	for( int i=0;i<numPixels;++i )
	{
		Color3f& c = map.coeffRef(i);
		samples_r.push_back(c[0]);
		samples_g.push_back(c[1]);
		samples_b.push_back(c[2]);
	}

	return Color3f(variance(samples_r, double(mean[0])), variance(samples_g, double(mean[1])), variance(samples_b, double(mean[2])));
}


Color3f mean(Bitmap &map)
{
	std::vector<double> samples_r;
	std::vector<double> samples_g;
	std::vector<double> samples_b;
	int numPixels = map.cols()*map.rows();
	for( int i=0;i<numPixels;++i )
	{
		Color3f& c = map.coeffRef(i);
		samples_r.push_back(c[0]);
		samples_g.push_back(c[1]);
		samples_b.push_back(c[2]);
	}

	return Color3f(mean(samples_r), mean(samples_g), mean(samples_b));
}




Bitmap filter( Bitmap& bm, Filter filter )
{
	Bitmap tmp( V2i( bm.cols(), bm.rows()) );


	int kshift = filter.size/2; // used to center the kernel at current pixel
	int resx = bm.cols();
	int resy = bm.rows();
	for( int j=0;j<resy;++j )
		for( int i=0;i<resx;++i )
		{
			Color3f c(0.0, 0.0, 0.0);
			// iterate over kernel
			int j_kstart = std::max( j-kshift, 0 );
			int j_kend = std::min( j+kshift, int(bm.rows())-1);
			int i_kstart = std::max( i-kshift, 0 );
			int i_kend = std::min(i+kshift, int(bm.cols())-1);
			double sum = 0.0;
			for( int kj=j_kstart;kj<=j_kend;++kj )
				for( int ki=i_kstart;ki<=i_kend;++ki )
				{
					double w = filter.weight( bm, i, j, ki, kj );
					c += bm.coeffRef(kj, ki)*w;
					sum += w;
				}

			// set center pixel
			tmp.coeffRef(j,i) = c/sum;
		}


	return tmp;
}



Filter gaussFilter(double stddev)
{
	Filter f;
	double variance = stddev*stddev;
	auto weight = [&](Bitmap& input, int i, int j, int ki, int kj)
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
	auto weight2 = [&](int i, int j, int ki, int kj)
	{
		double x = i+0.5;
		double y = j+0.5;
		double kx = ki+0.5;
		double ky = kj+0.5;
		double dx = kx-x;
		double dy = ky-y;

		// we dont have a normalization constant since we normalize weights
		// over the pixelneighbourhood anyways...
		//return std::exp( -(dx*dx + dy*dy)/(2*variance) );
		return 1.0;
	};

	f.weight = weight;
	f.weight2 = weight2;

	// kernelsize for gaussian:
	f.size = int(6*stddev - 1);
	return f;
}


// alpha and beta are smoothing parameters per feature (alpha=distance, beta=color similarity)
Filter bilateralFilter(double alpha, double beta)
{
	Filter f;
	auto weight = [alpha, beta](Bitmap& input, int i, int j, int ki, int kj)
	{
		double x = i+0.5;
		double y = j+0.5;
		double kx = ki+0.5;
		double ky = kj+0.5;
		double dx = kx-x;
		double dy = ky-y;
		double la = input.coeffRef(j,i).getLuminance();
		double lb = input.coeffRef(kj,ki).getLuminance();
		double dc = std::abs( la - lb );

		// we dont have a normalization constant since we normalize weights
		// over the pixelneighbourhood anyways...
		return std::exp( -(dx*dx + dy*dy)/(2*alpha*alpha) -(dc*dc)/(2*beta*beta) );
	};
	f.weight = weight;

	// kernelsize for gaussian:
	f.size = int(6*std::max(alpha, beta) - 1);
	return f;
}


