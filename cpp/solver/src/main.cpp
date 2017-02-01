#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/eigen.h>
#include <iostream>
#include <math/common.h>
#include <math/vector.h>
#include <math/ray.h>
#include <util\dda2d.h>


namespace py = pybind11;


int index_2d(int i, int j, int res_x)
{
	return i*res_x + j;
}

int index_M_2d(int i, int j, int res_x)
{
	return i*res_x*4 + j*4;
}

void iterate_2d(py::array_t<double> phi_array, // phi grid
	py::array_t<double> D_array, // diffusion coefficient grid
	py::array_t<double> Q_array, // grid of emission values
	py::array_t<double> phi_boundary_array, // boundary values for phi
	py::array_t<double> D_boundary_array, // boundary values for D
	double h, // voxelsize (in x and y)
	bool debug = false
)
{
	py::buffer_info phi_array_info = phi_array.request();
	py::buffer_info D_array_info = D_array.request();
	py::buffer_info Q_array_info = Q_array.request();
	py::buffer_info phi_boundary_array_info = phi_boundary_array.request();
	py::buffer_info D_boundary_array_info = D_boundary_array.request();

	auto phi = static_cast<double *>(phi_array_info.ptr);
	auto D = static_cast<double *>(D_array_info.ptr);
	auto Q = static_cast<double *>(Q_array_info.ptr);
	auto phi_boundary = static_cast<double *>(phi_boundary_array_info.ptr);
	auto D_boundary = static_cast<double *>(D_boundary_array_info.ptr);

	int res_y = phi_array_info.shape[0];
	int res_x = phi_array_info.shape[1];

	// update boundary conditions ---
	// here we copy over just the ghost cell values from the boundary grids
	for (int i = 0; i < res_y; ++i)
	{
		int idx_0 = index_2d(i, 0, res_x);
		int idx_1 = index_2d(i, res_x - 1, res_x);
		phi[idx_0] = phi_boundary[idx_0];
		phi[idx_1] = phi_boundary[idx_1];
		D[idx_0] = D_boundary[idx_0];
		D[idx_1] = D_boundary[idx_1];
	}
	for (int i = 0; i < res_x; ++i)
	{
		int idx_0 = index_2d(0, i, res_x);
		int idx_1 = index_2d(res_y - 1, i, res_x);
		phi[idx_0] = phi_boundary[idx_0];
		phi[idx_1] = phi_boundary[idx_1];
		D[idx_0] = D_boundary[idx_0];
		D[idx_1] = D_boundary[idx_1];
	}

	// now iterate over all inner cells and update phi and D
	for (int i=1; i < res_y-1; ++i)
		for (int j=1; j < res_x-1; ++j)
		{
			int idx = index_2d(i, j, res_x);
			int idx_xp = index_2d(i, j+1, res_x);
			int idx_xm = index_2d(i, j-1, res_x);
			int idx_yp = index_2d(i-1, j, res_x);
			int idx_ym = index_2d(i+1, j, res_x);

			// compute phi gradient at cell center-- -
			V2d grad_phi;
			grad_phi[0] = (phi[idx_xp] - phi[idx_xm]) / (2.0*h);
			grad_phi[1] = (phi[idx_yp] - phi[idx_ym]) / (2.0*h);

			// Here is where we need to handle the case where D becomes zero. This would cause a division by zero during update of phi.
			// NB: the threshold value for phi has an effect on the accuracy of the solution. It has been hand-picked for the pointsource test.
			double dc = std::max(phi[idx], 0.0004) / std::max(grad_phi.norm(), 1.0e-4);

			// compute diffusion coefficients at cell faces-- -
			double D_xph = (dc + D[idx_xp])*0.5;
			double D_xmh = (dc + D[idx_xm])*0.5;
			double D_yph = (dc + D[idx_yp])*0.5;
			double D_ymh = (dc + D[idx_ym])*0.5;

			double numerator = 0.0;
			numerator += D_xph*phi[idx_xp];
			numerator += D_xmh*phi[idx_xm];
			numerator += D_yph*phi[idx_yp];
			numerator += D_ymh*phi[idx_ym];
			numerator /= h;
			numerator += h*Q[idx];

			double denominator = 0.0;
			denominator += (D_xph + D_xmh + D_yph + D_ymh) /h;

			phi[idx] = numerator / denominator;
			D[idx] = dc;

			if( i==76 && j==93 && debug == true )
			{
				std::cout << "debug!!!\n";
				std::cout << "phi=" << phi[idx];
				std::cout << "grad_phi=" << grad_phi.toString() << std::endl;
			}
		}
}



void iterate_2d_anisotropic(py::array_t<double> phi_array, // phi grid
	py::array_t<double> D_array, // diffusion coefficient grid
	py::array_t<double> Q_array, // grid of emission values
	py::array_t<double> M_x_array, // grid of 2x2 matrices which represent anisotropy on faces perpendicular to x direction
	py::array_t<double> M_y_array, // grid of 2x2 matrices which represent anisotropy on faces perpendicular to y direction
	py::array_t<double> phi_boundary_array, // boundary values for phi
	py::array_t<double> D_boundary_array, // boundary values for D
	double h // voxelsize (in x and y)
)
{
	py::buffer_info phi_array_info = phi_array.request();
	py::buffer_info D_array_info = D_array.request();
	py::buffer_info Q_array_info = Q_array.request();
	py::buffer_info M_x_array_info = M_x_array.request();
	py::buffer_info M_y_array_info = M_y_array.request();
	py::buffer_info phi_boundary_array_info = phi_boundary_array.request();
	py::buffer_info D_boundary_array_info = D_boundary_array.request();

	auto phi = static_cast<double *>(phi_array_info.ptr);
	auto D = static_cast<double *>(D_array_info.ptr);
	auto Q = static_cast<double *>(Q_array_info.ptr);
	auto M_x = static_cast<double *>(M_x_array_info.ptr);
	auto M_y = static_cast<double *>(M_y_array_info.ptr);
	auto phi_boundary = static_cast<double *>(phi_boundary_array_info.ptr);
	auto D_boundary = static_cast<double *>(D_boundary_array_info.ptr);

	int res_y = phi_array_info.shape[0];
	int res_x = phi_array_info.shape[1];


	// update boundary conditions ---
	// here we copy over just the ghost cell values from the boundary grids
	for (int i = 0; i < res_y; ++i)
	{
		int idx_0 = index_2d(i, 0, res_x);
		int idx_1 = index_2d(i, res_x - 1, res_x);
		phi[idx_0] = phi_boundary[idx_0];
		phi[idx_1] = phi_boundary[idx_1];
		D[idx_0] = D_boundary[idx_0];
		D[idx_1] = D_boundary[idx_1];
	}
	for (int i = 0; i < res_x; ++i)
	{
		int idx_0 = index_2d(0, i, res_x);
		int idx_1 = index_2d(res_y - 1, i, res_x);
		phi[idx_0] = phi_boundary[idx_0];
		phi[idx_1] = phi_boundary[idx_1];
		D[idx_0] = D_boundary[idx_0];
		D[idx_1] = D_boundary[idx_1];
	}

	// now iterate over all inner cells and update phi and D
	for (int i = 1; i < res_y - 1; ++i)
		for (int j = 1; j < res_x - 1; ++j)
		{
			int idx = index_2d(i, j, res_x);
			int idx_xp = index_2d(i, j + 1, res_x);
			int idx_xm = index_2d(i, j - 1, res_x);
			int idx_yp = index_2d(i - 1, j, res_x);
			int idx_ym = index_2d(i + 1, j, res_x);

			// we will need the diagonal voxels as well
			int idx_xpyp = index_2d(i-1, j+1, res_x);
			int idx_xpym = index_2d(i+1, j+1, res_x);
			int idx_xmyp = index_2d(i-1, j-1, res_x);
			int idx_xmym = index_2d(i+1, j-1, res_x);

			// compute phi gradient at cell center-- -
			V2d grad_phi;
			grad_phi[0] = (phi[idx_xp] - phi[idx_xm]) / (2.0*h);
			grad_phi[1] = (phi[idx_yp] - phi[idx_ym]) / (2.0*h);

			// retrieve M matrices at cell faces ---
			int idx_M_xph = index_M_2d(i, j+1, res_x+1);
			int idx_M_xmh = index_M_2d(i, j, res_x+1);
			int idx_M_yph = index_M_2d(i, j, res_x+1);
			int idx_M_ymh = index_M_2d(i+1, j, res_x+1);
			M22d M_xph;M_xph << M_x[idx_M_xph+0], M_x[idx_M_xph+1], M_x[idx_M_xph+2], M_x[idx_M_xph+3];
			M22d M_xmh;M_xmh << M_x[idx_M_xmh+0], M_x[idx_M_xmh+1], M_x[idx_M_xmh+2], M_x[idx_M_xmh+3];
			M22d M_yph;M_yph << M_y[idx_M_yph+0], M_y[idx_M_yph+1], M_y[idx_M_yph+2], M_y[idx_M_yph+3];
			M22d M_ymh;M_ymh << M_y[idx_M_ymh+0], M_y[idx_M_ymh+1], M_y[idx_M_ymh+2], M_y[idx_M_ymh+3];

			// compute anisotropic diffusion coefficient matrices at cell faces
			// Here is where we need to handle the case where D becomes zero. This would cause a division by zero during update of phi.
			// NB: the threshold value for phi has an effect on the accuracy of the solution. It has been hand-picked for the pointsource test.
			double dc = std::max(phi[idx], 0.0004) / std::max(grad_phi.norm(), 1.0e-4);

			// compute diffusion coefficient matrices at cell faces ---
			M22d D_xph = (dc + D[idx_xp])*0.5*M_xph;
			M22d D_xmh = (dc + D[idx_xm])*0.5*M_xmh;
			M22d D_yph = (dc + D[idx_yp])*0.5*M_yph;
			M22d D_ymh = (dc + D[idx_ym])*0.5*M_ymh;

			// update rule for phi ---
			double hinv = 1.0/h;
			double h4inv = 1.0/(4.0*h);

			double numerator = 0.0;
			numerator += h4inv*( D_xmh(0,1) + D_ymh(1,0) )*phi[idx_xmym];
			numerator += hinv*(D_xmh(0,0) - 0.25*D_yph(1,0) + 0.25*D_ymh(1,0))*phi[idx_xm];
			numerator += -h4inv*(D_xmh(0,1)+D_yph(1,0))*phi[idx_xmyp];
			numerator += hinv*(-0.25*D_xph(0,1)+0.25*D_xmh(0,1) + D_ymh(1, 1))*phi[idx_ym];
			numerator += hinv*(0.25*D_xph(0,1)-0.25*D_xmh(0,1)+D_yph(1,1))*phi[idx_yp];
			numerator += -h4inv*(D_xph(0,1) + D_ymh(1,0))*phi[idx_xpym];
			numerator += hinv*(D_xph(0,0)+0.25*D_yph(1,0)-0.25*D_ymh(1,0))*phi[idx_xp];
			numerator += h4inv*(D_xph(0,1) + D_yph(1,0))*phi[idx_xpyp];
			numerator += h*Q[idx];

			double denominator = 0.0;
			denominator += hinv*(D_xph(0,0)+D_xmh(0,0)+D_yph(1,1)+D_ymh(1,1));
						
	
			phi[idx] = numerator / denominator;
			D[idx] = dc;
		}
}

// p0 and p1 are given in voxelspace
void rasterize_line( py::array_t<double> pixel_array, Eigen::Vector2d& p0, Eigen::Vector2d& p1)
{
	py::buffer_info pixel_array_info = pixel_array.request();
	auto pixel = static_cast<double *>(pixel_array_info.ptr);
	V2i res( pixel_array_info.shape[0], pixel_array_info.shape[1] );

	V2i pix0;
	pix0.x() = std::min(res.x(), std::max(int(p0.x()), 0));
	pix0.y() = std::min(res.y(), std::max(int(p0.y()), 0));
	V2i pix1;
	pix1.x() = std::min(res.x(), std::max(int(p1.x()), 0));
	pix1.y() = std::min(res.y(), std::max(int(p1.y()), 0));

	DrawPixelFunction dpf = [&](short i, short j, short baseColor)
	{
		int pixel_index = index_2d( res.y()-j-1, i, res.x() );
		pixel[pixel_index]	= 1.0-baseColor/256.0;
	};
	short BaseColor = 0;
	short NumLevels = 256;
	short IntensityBits = 8;
	draw_wu_line(pix0.x(), pix0.y(), pix1.x(), pix1.y(), BaseColor, NumLevels, IntensityBits, dpf);
}



P2i toVoxelIndex( P2d pWS, const Box2d& bound, const V2i& res )
{
	// convert from world to localspace
	P2d pLS = (pWS-bound.min).cwiseQuotient(bound.getExtents());
	// convert from local to voxelspace
	P2d pVS( pLS.x()*res.x(), pLS.y()*res.y() );
	P2i index;
	index.x() = std::min(res.x()-1, std::max(int(pVS.x()), 0));
	index.y() = std::min(res.y()-1, std::max(int(pVS.y()), 0));
	return index;
}

P2d voxelToWorld( const P2d& pVS, const Box2d& bound, const V2i& res  )
{
	P2d pLS( pVS.x()/res.x(), pVS.y()/res.y() );
	return pLS.cwiseProduct(bound.getExtents()) + bound.min;
}

void trace_2d( py::array_t<double> pixel_array, Eigen::Vector2d& o, Eigen::Vector2d& d, Eigen::Vector2d& bound_min, Eigen::Vector2d& bound_max )
{
	py::buffer_info pixel_array_info = pixel_array.request();
	auto pixel = static_cast<double *>(pixel_array_info.ptr);
	V2i res( pixel_array_info.shape[0], pixel_array_info.shape[1] );

	// this is for line drawing ---
	//double energy = 0.0;
	DrawPixelFunction dpf = [&](short i, short j, short baseColor)
	{
		int pixel_index = index_2d( res.y()-j-1, i, res.x() );
		double nrg = 1.0-baseColor/256.0;
		pixel[pixel_index]	+= nrg;
		//energy += nrg;
	};


	// now do the tracing ---
	int depth = 0;

	Box2d bound( bound_min, bound_max );
	Ray2d ray(o, d);

	while(depth < 2)
	{
		double tnear, tfar;
		if( bound.rayIntersect(ray, tnear, tfar) )
		{
			if(tnear > 0)
				// we are outside the domain-invalid sample
				break;
			//std::cout << "tnear=" << tnear << std::endl;

			// now splat the line segment --
			P2i pix0 = toVoxelIndex(ray.o, bound, res);
			P2i pix1 = toVoxelIndex(ray(tfar), bound, res);
			draw_wu_line(pix0.x(), pix0.y(), pix1.x(), pix1.y(), 0, 256, 8, dpf);

			// reflect on the right border
			if(pix1.x() == res.x()-1)
			{
				ray = Ray2d(voxelToWorld( P2d(pix1.x()+0.5, pix1.y()+0.5), bound, res ), V2d(-ray.d.x(), ray.d.y()));
			}else
				// we left the domain
				break;
		}else
		{
			std::cout << "trace_2d: error: no intersection\n";
		}
		++depth;
	}
}


void toSRGB( py::array_t<double> pixel_array )
{
	py::buffer_info pixel_array_info = pixel_array.request();
	auto pixel = static_cast<double *>(pixel_array_info.ptr);
	V2i res( pixel_array_info.shape[0], pixel_array_info.shape[1] );

	int numPixels = res.x()*res.y();
	for( int i=0;i<numPixels;++i )
	{
		double value = pixel[i];
		double new_value;
		if (value <= 0.0031308f)
			new_value = 12.92f * value;
		else
			new_value = (1.0f + 0.055f)*std::pow(value, 1.0f/2.4f) -  0.055f;

		// clamp
		new_value = std::max(std::min(new_value, 1.0), 0.0);
		
		pixel[i] = new_value;
	}

}

PYBIND11_PLUGIN(solver)
{
	py::module m("solver", "some inner loops in c++");

	m.def("iterate_2d", &iterate_2d);
	m.def("iterate_2d_anisotropic", &iterate_2d_anisotropic);
	m.def("rasterize_line", &rasterize_line);
	m.def("trace_2d", &trace_2d);
	m.def("toSRGB", &toSRGB);

	return m.ptr();
}