#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/eigen.h>
#include <iostream>
#include <math/common.h>
#include <math/vector.h>
#include <math/ray.h>
#include <util\dda2d.h>


namespace py = pybind11;


int index_2d(int i, int j, int res_x, int stride = 1)
{
	return i*res_x*stride + j*stride;
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
		}
}



// this function is used to solve isotropic absorption equation for delta radiance distribution
void iterate_2d_isotropic_absorption(py::array_t<double> phi_array, // phi grid
	py::array_t<double> D_array, // diffusion coefficient grid
	py::array_t<double> sigma_t_array, // extinction coefficient grid
	py::array_t<double> Q_array, // grid of emission values
	py::array_t<double> phi_boundary_array, // boundary values for phi
	py::array_t<double> D_boundary_array, // boundary values for D
	double h, // voxelsize (in x and y)
	bool debug = false
)
{
	py::buffer_info phi_array_info = phi_array.request();
	py::buffer_info D_array_info = D_array.request();
	py::buffer_info sigma_t_array_info = sigma_t_array.request();
	py::buffer_info Q_array_info = Q_array.request();
	py::buffer_info phi_boundary_array_info = phi_boundary_array.request();
	py::buffer_info D_boundary_array_info = D_boundary_array.request();

	auto phi = static_cast<double *>(phi_array_info.ptr);
	auto D = static_cast<double *>(D_array_info.ptr);
	auto sigma_t = static_cast<double *>(sigma_t_array_info.ptr);
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

			double denominator = h*sigma_t[idx];
			denominator += (D_xph + D_xmh + D_yph + D_ymh) /h;

			phi[idx] = numerator / denominator;
			D[idx] = dc;
		}
}


void iterate_2d_diffusion(py::array_t<double> phi_array, // phi grid
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

			

			// compute diffusion coefficients at cell faces-- -
			double dc = D[idx];
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
			denominator += (D_xph + D_xmh + D_yph + D_ymh)/h;

			phi[idx] = numerator / denominator;
		}
}

void iterate_2d_p1(
	py::array_t<double> phi_array, // phi grid
	py::array_t<double> E_array, // E grid
	py::array_t<double> D_array, // diffusion coefficient grid
	py::array_t<double> Q_array, // grid of emission values
	py::array_t<double> phi_boundary_array, // boundary values for phi
	py::array_t<double> D_boundary_array, // boundary values for D
	double h, // voxelsize (in x and y)
	bool debug = false
)
{
	py::buffer_info phi_array_info = phi_array.request();
	py::buffer_info E_array_info = E_array.request();
	py::buffer_info D_array_info = D_array.request();
	py::buffer_info Q_array_info = Q_array.request();
	py::buffer_info phi_boundary_array_info = phi_boundary_array.request();
	py::buffer_info D_boundary_array_info = D_boundary_array.request();

	auto phi = static_cast<double *>(phi_array_info.ptr);
	auto E = static_cast<double *>(E_array_info.ptr);
	auto D = static_cast<double *>(D_array_info.ptr);
	auto Q = static_cast<double *>(Q_array_info.ptr);
	auto phi_boundary = static_cast<double *>(phi_boundary_array_info.ptr);
	auto D_boundary = static_cast<double *>(D_boundary_array_info.ptr);

	int res_y = phi_array_info.shape[0];
	int res_x = phi_array_info.shape[1];

	double hInv = 1.0/h;

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

			int idx_E = index_2d(i, j, res_x, 2);
			int idx_E_xp = index_2d(i, j+1, res_x, 2);
			int idx_E_xm = index_2d(i, j-1, res_x, 2);
			int idx_E_yp = index_2d(i-1, j, res_x, 2);
			int idx_E_ym = index_2d(i+1, j, res_x, 2);


			// first we update Ex
			double Ex = E[idx_E+0];
			double Ey = E[idx_E+1];

			//for( int k=0;k<10;++k )
			{
				Ey = h*Q[idx] + E[idx_E_xm+0] - E[idx_E+0] + E[idx_E_ym+1];
				Ex = h*Q[idx] + E[idx_E_xm+0] - E[idx_E+1] + E[idx_E_ym+1];
			}


			// then we update phi
			//double phi_ij = h*Ex + phi[idx_xm];

			// finally we upadte Ey
			//double Ey = hInv*( phi_ij - phi[idx_ym] );

			// write new values
			E[idx_E+0] = Ex;
			E[idx_E+1] = Ey;
			//phi[idx] = phi_ij;
		}
}

std::string iterate_2d_anisotropic_absorption(
	py::array_t<double> phi_array, // phi grid
	py::array_t<double> mu0_sigma_t_array, // zero moment of the extinction coefficient at cell centers
	py::array_t<double> Q_array, // grid of emission values
	py::array_t<double> M_x_array, // grid of 2x2 matrices which represent anisotropy on faces perpendicular to x direction
	py::array_t<double> M_y_array, // grid of 2x2 matrices which represent anisotropy on faces perpendicular to y direction
	py::array_t<double> phi_boundary_array, // boundary values for phi
	double h, // voxelsize (in x and y)
	bool debug_global = false
)
{
	std::ostringstream oss;
	py::buffer_info phi_array_info = phi_array.request();
	py::buffer_info mu0_sigma_t_array_info = mu0_sigma_t_array.request();
	py::buffer_info M_x_array_info = M_x_array.request();
	py::buffer_info M_y_array_info = M_y_array.request();
	py::buffer_info Q_array_info = Q_array.request();
	py::buffer_info phi_boundary_array_info = phi_boundary_array.request();

	auto phi = static_cast<double *>(phi_array_info.ptr);
	auto mu0_sigma_t = static_cast<double *>(mu0_sigma_t_array_info.ptr);
	auto M_x = static_cast<double *>(M_x_array_info.ptr);
	auto M_y = static_cast<double *>(M_y_array_info.ptr);
	auto Q = static_cast<double *>(Q_array_info.ptr);
	auto phi_boundary = static_cast<double *>(phi_boundary_array_info.ptr);

	int res_y = phi_array_info.shape[0];
	int res_x = phi_array_info.shape[1];

	double hinv = 1.0/h;
	double h4inv = 1.0/(4.0*h);

	// update boundary conditions ---
	// here we copy over just the ghost cell values from the boundary grids
	for (int i = 0; i < res_y; ++i)
	{
		int idx_0 = index_2d(i, 0, res_x);
		int idx_1 = index_2d(i, res_x - 1, res_x);
		phi[idx_0] = phi_boundary[idx_0];
		phi[idx_1] = phi_boundary[idx_1];
	}
	for (int i = 0; i < res_x; ++i)
	{
		int idx_0 = index_2d(0, i, res_x);
		int idx_1 = index_2d(res_y - 1, i, res_x);
		phi[idx_0] = phi_boundary[idx_0];
		phi[idx_1] = phi_boundary[idx_1];
	}

	// now iterate over all inner cells and update phi and D
	for (int i=1; i < res_y-1; ++i)
		for (int j=1; j < res_x-1; ++j)
		{
			bool debug = false;
			//if( (j>=80)&&(j<=81)&&(i==80) )
			if( (j>=0)&&(j<=81)&&(i==80) && debug_global )
				debug = true;

			if(debug)
			{
				oss << "-------- i=" << i << " j=" << j << std::endl;
			}

			int idx = index_2d(i, j, res_x);
			int idx_xp = index_2d(i, j+1, res_x);
			int idx_xpyp = index_2d(i-1, j+1, res_x);
			int idx_xm = index_2d(i, j-1, res_x);
			int idx_xmyp = index_2d(i-1, j-1, res_x);
			int idx_yp = index_2d(i-1, j, res_x);
			int idx_ym = index_2d(i+1, j, res_x);
			int idx_xpym = index_2d(i+1, j+1, res_x);
			int idx_xmym = index_2d(i+1, j-1, res_x);

			// compute phi at cell faces ---
			double phi2[4];
			phi2[0] = (phi[idx] + phi[idx_xm])/2.0;
			phi2[1] = (phi[idx_xp] + phi[idx])/2.0;
			phi2[2] = (phi[idx] + phi[idx_ym])/2.0;
			phi2[3] = (phi[idx_yp] + phi[idx])/2.0;

			// compute phi gradient at cell faces ---
			// index 0 == left face (i-1/2)
			// index 1 == right face (i+1/2)
			// index 2 == bottom face (j-1/2)
			// index 3 == top face (j+1/2)
			V2d grad_phi[4];
			grad_phi[0].x() = (phi[idx] - phi[idx_xm])/h;
			grad_phi[0].y() = (phi[idx_xmyp] + phi[idx_yp] - phi[idx_xmym] - phi[idx_ym])/(4.0*h);

			grad_phi[1].x() = (phi[idx_xp] - phi[idx])/h;
			grad_phi[1].y() = (phi[idx_yp] + phi[idx_xpyp] - phi[idx_ym] - phi[idx_xpym])/(4.0*h);

			grad_phi[2].x() = (phi[idx_xp] + phi[idx_xpym] - phi[idx_xm] - phi[idx_xmym])/(4.0*h);
			grad_phi[2].y() = (phi[idx] - phi[idx_ym])/h;

			grad_phi[3].x() = (phi[idx_xpyp] + phi[idx_xp] - phi[idx_xmyp] - phi[idx_xm])/(4.0*h);
			grad_phi[3].y() = (phi[idx_yp] - phi[idx])/h;

			// retrieve M matrices at cell faces ---
			int idx_M[4];
			idx_M[0] = index_M_2d(i, j, res_x+1);
			idx_M[1] = index_M_2d(i, j+1, res_x+1);
			idx_M[2] = index_M_2d(i+1, j, res_x+1);
			idx_M[3] = index_M_2d(i, j, res_x+1);

			M22d M[4];
			M[0] << M_x[idx_M[0]+0], M_x[idx_M[0]+1], M_x[idx_M[0]+2], M_x[idx_M[0]+3];
			M[1] << M_x[idx_M[1]+0], M_x[idx_M[1]+1], M_x[idx_M[1]+2], M_x[idx_M[1]+3];
			M[2] << M_y[idx_M[2]+0], M_y[idx_M[2]+1], M_y[idx_M[2]+2], M_y[idx_M[2]+3];
			M[3] << M_y[idx_M[3]+0], M_y[idx_M[3]+1], M_y[idx_M[3]+2], M_y[idx_M[3]+3];

			// compute D at cell faces ---
			M22d D[4];
			for(int k=0;k<4;++k)
			{
				double coeff = std::max(phi2[k], 0.0004) / std::max((M[k]*grad_phi[k]).norm(), 1.0e-4);
				//double coeff = std::max(phi2[k], 0.0004) / std::max((grad_phi[k]).norm(), 1.0e-4);
				D[k] = M[k]*coeff;
			}

			double numerator = 0.0;
			/*
			numerator += -D[1].coeffRef(0,0)*phi[idx_xp];
			numerator += -D[0].coeffRef(0,0)*phi[idx_xm];
			numerator += -D[3].coeffRef(1,1)*phi[idx_yp];
			numerator += -D[2].coeffRef(1,1)*phi[idx_ym];
			*/

			if(debug)
			{





			}

			numerator += -0.25*(D[0].coeffRef(0,1)+D[2].coeffRef(1,0))*phi[idx_xmym];
			if(debug)
			{
				oss << "\tnumerator=" << numerator << std::endl;
				oss << "\tphi[idx_xmym]" << phi[idx_xmym] << std::endl;
			}
			numerator += -(D[0].coeffRef(0,0)-0.25*D[3].coeffRef(1,0) + 0.25*D[2].coeffRef(1,0))*phi[idx_xm];
			if(debug)
			{
				oss << "\tnumerator=" << numerator << std::endl;
				oss << "\tphi[idx_xm]" << phi[idx_xm] << std::endl;
			}
			numerator += 0.25*(D[0].coeffRef(0,1)+D[3].coeffRef(1,0))*phi[idx_xmyp];
			if(debug)
			{
				oss << "\tnumerator=" << numerator << std::endl;
				oss << "\tphi[idx_xmyp]" << phi[idx_xmyp] << std::endl;
			}
			numerator += -(-0.25*D[1].coeffRef(0,1)+0.25*D[0].coeffRef(0,1)+D[2].coeffRef(1,1))*phi[idx_ym];
			if(debug)
			{
				oss << "\tnumerator=" << numerator << std::endl;
				oss << "\tphi[idx_ym]" << phi[idx_ym] << std::endl;
			}
			numerator += -(0.25*D[1].coeffRef(0,1)-0.25*D[0].coeffRef(0,1)+D[3].coeffRef(1,1))*phi[idx_yp];
			if(debug)
			{
				oss << "\tnumerator=" << numerator << std::endl;
				oss << "\tphi[idx_yp]" << phi[idx_yp] << std::endl;
			}
			numerator += 0.25*(D[1].coeffRef(0,1) + D[2].coeffRef(1,0))*phi[idx_xpym];
			if(debug)
			{
				oss << "\tnumerator=" << numerator << std::endl;
				oss << "\tphi[idx_xpym]" << phi[idx_xpym] << std::endl;
			}
			numerator += -(D[1].coeffRef(0,0) + 0.25*D[3].coeffRef(1,0) - 0.25*D[2].coeffRef(1,0))*phi[idx_xp];
			if(debug)
			{
				oss << "\tnumerator=" << numerator << std::endl;
				oss << "\tphi[idx_xp]" << phi[idx_xp] << std::endl;
			}
			numerator += -0.25*(D[1].coeffRef(0,1)+D[3].coeffRef(1,0))*phi[idx_xpyp];
			if(debug)
			{
				oss << "\tnumerator=" << numerator << std::endl;
				oss << "\tphi[idx_xpyp]" << phi[idx_xpyp] << std::endl;
			}

			numerator /= h;
			if(debug)
				oss << "\tnumerator=" << numerator << std::endl;

			numerator += -h*Q[idx];
			if(debug)
				oss << "\tnumerator=" << numerator << std::endl;

			double denominator = -h*mu0_sigma_t[idx];
			denominator += -(D[1].coeffRef(0,0) + D[0].coeffRef(0,0) + D[3].coeffRef(1,1) + D[2].coeffRef(1,1))/h;

			phi[idx] = numerator / denominator;
			if(debug)
			{
				oss << "\tphi=" << phi[idx] << std::endl;
				oss << "\tnumerator=" << numerator << std::endl;
				oss << "\tdenominator=" << denominator << std::endl;
				for(int k=0;k<4;++k)
				{
					double coeff = std::max(phi2[k], 0.0004) / std::max((M[k]*grad_phi[k]).norm(), 1.0e-4);
					oss << "\tD["<< k << "]=" << D[k].coeffRef(0,0) << " " << D[k].coeffRef(0,1) << " " << D[k].coeffRef(1,0) << " " << D[k].coeffRef(1,1) << std::endl;
				}
			}
		}


	return oss.str();

	/*
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
	*/
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
	m.def("iterate_2d_isotropic_absorption", &iterate_2d_isotropic_absorption);
	m.def("iterate_2d_anisotropic_absorption", &iterate_2d_anisotropic_absorption);
	m.def("iterate_2d_diffusion", &iterate_2d_diffusion);
	m.def("iterate_2d_p1", &iterate_2d_p1);
	m.def("toSRGB", &toSRGB);

	return m.ptr();
}