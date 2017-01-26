#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <iostream>
#include <math/common.h>
#include <math/vector.h>


namespace py = pybind11;


int index_2d(int i, int j, int res_x)
{
	return i*res_x + j;
}


void iterate_2d(py::array_t<double> phi_array, // phi grid
	py::array_t<double> D_array, // diffusion coefficient grid
	py::array_t<double> Q_array, // grid of emission values
	py::array_t<double> phi_boundary_array, // boundary values for phi
	py::array_t<double> D_boundary_array, // boundary values for D
	double h // voxelsize (in x and y)
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




PYBIND11_PLUGIN(solver)
{
	py::module m("solver", "some inner loops in c++");

	m.def("iterate_2d", &iterate_2d);

	return m.ptr();
}