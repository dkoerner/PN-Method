#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <complex>
#include <iostream>

#include <math/common.h>
#include <math/vector.h>
#include <math/ray.h>

#include<common/Domain.h>
#include<common/GridLocation.h>
#include<field/VoxelGrid.h>
#include<field/Constant.h>
#include<field/SHEXP.h>


namespace py = pybind11;





V3d to_V3d( py::array array )
{
	auto data = static_cast<double *>(array.request().ptr);
	return V3d(data[0], data[1], data[2]);
}

P3d to_P3d( py::array array )
{
	auto data = static_cast<double *>(array.request().ptr);
	return P3d(data[0], data[1], data[2]);
}

P2d to_P2d( py::array array )
{
	auto data = static_cast<double *>(array.request().ptr);
	return P2d(data[0], data[1]);
}



PYBIND11_MODULE(pnbuilder_cpp, m)
{
	// Domain ============================================================
	py::class_<Domain> class_domain(m, "Domain");
	class_domain
	.def("__init__",
	[](Domain &m, py::array size_array, py::array resolution_array, py::array offset_array)
	{
		py::buffer_info size_info = size_array.request();
		py::buffer_info resolution_info = resolution_array.request();
		py::buffer_info offset_info = offset_array.request();

		// Some sanity checks ...

		if (resolution_info.ndim != 1 || offset_info.ndim != 1|| size_info.ndim != 1)
			throw std::runtime_error("Incompatible dimension!");

		if (resolution_info.shape[0] != 2 || offset_info.shape[0] != 2 || size_info.shape[0] != 2)
			throw std::runtime_error("Incompatible vector size!");

		auto size_data = static_cast<double *>(size_info.ptr);
		auto resolution_data = static_cast<int *>(resolution_info.ptr);
		auto offset_data = static_cast<double *>(offset_info.ptr);

		new (&m) Domain( V2d(size_data[0], size_data[1]),
						 V2i(resolution_data[0], resolution_data[1]),
						 V2d(offset_data[0], offset_data[1]));
	})
	.def("resolution",
	[](Domain &m)
	{
		V2i resolution = m.resolution();
		py::array_t<int> a({ 2 });
		a.mutable_data()[0] = resolution[0];
		a.mutable_data()[1] = resolution[1];
		return a;

	})
	.def("numVoxels",
	[](Domain &m)
	{
		return m.numVoxels();
	});

	// GridLocation ============================================================
	py::class_<GridLocation> class_gridlocation(m, "GridLocation");
	class_gridlocation
	.def("__init__",
	[](GridLocation &m, Domain &domain, py::array voxel_array, py::array offset_array)
	{
		py::buffer_info voxel_info = voxel_array.request();
		py::buffer_info offset_info = offset_array.request();

		// Some sanity checks ...
		if (voxel_info.ndim != 1 || offset_info.ndim != 1)
			throw std::runtime_error("Incompatible dimension!");

		if (voxel_info.shape[0] != 2 || offset_info.shape[0] != 2 )
			throw std::runtime_error("Incompatible vector size!");

		auto voxel_data = static_cast<int *>(voxel_info.ptr);
		auto offset_data = static_cast<int *>(offset_info.ptr);

		new (&m) GridLocation(  domain,
								V2i(voxel_data[0], voxel_data[1]),
								V2i(offset_data[0], offset_data[1]));
	})
	.def("getVoxel",
	[](GridLocation &m)
	{
		V2i voxel = m.getVoxel();
		py::array_t<int> a({ 2 });
		a.mutable_data()[0] = voxel[0];
		a.mutable_data()[1] = voxel[1];
		return a;

	})
	.def("getOffset",
	[](GridLocation &m)
	{
		V2i offset = m.getOffset();
		py::array_t<int> a({ 2 });
		a.mutable_data()[0] = offset[0];
		a.mutable_data()[1] = offset[1];
		return a;
	})
	.def("getPWS", &GridLocation::getPWS )
	.def("getShiftedLocation",
	[](GridLocation &m, py::array offset_array)
	{
		auto offset_data = static_cast<int *>(offset_array.request().ptr);
		//std::cout << "offset_data: " << offset_data[0] << " " << offset_data[1]  << std::endl;
		return m.getShiftedLocation(V2i(offset_data[0], offset_data[1]));
	});


	// Field ============================================================
	py::class_<Field, Field::Ptr> class_field(m, "Field");
	class_field
	.def("__call__",[](Field &m, py::array pWS_array){return m.eval(to_P2d(pWS_array));})
	.def("dx", [](Field &m, py::array pWS_array){return m.dx(to_P2d(pWS_array));})
	.def("dxdx", [](Field &m, py::array pWS_array){return m.dxdx(to_P2d(pWS_array));})
	.def("dxdy", [](Field &m, py::array pWS_array){return m.dxdy(to_P2d(pWS_array));})
	.def("dy", [](Field &m, py::array pWS_array){return m.dy(to_P2d(pWS_array));})
	.def("dydy", [](Field &m, py::array pWS_array){return m.dydy(to_P2d(pWS_array));})
	.def("dydx", [](Field &m, py::array pWS_array){return m.dydx(to_P2d(pWS_array));})
	.def("dz", [](Field &m, py::array pWS_array){return m.dz(to_P2d(pWS_array));})
	;

	// VoxelGrid ============================================================
	py::class_<VoxelGrid, VoxelGrid::Ptr> class_VoxelGrid(m, "VoxelGrid", class_field);
	class_VoxelGrid
	.def("__init__",
	[](VoxelGrid &m, py::array b, Domain& domain, py::array offset_array)
	{
		typedef Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic> Strides;

		// Request a buffer descriptor from Python
		py::buffer_info info = b.request();
		py::buffer_info offset_info = offset_array.request();

		// Some sanity checks ...
		//std::cout << "format=" << info.format << std::endl;
		if (info.format != py::format_descriptor<std::complex<double>>::format())
			throw std::runtime_error("Incompatible format: expected a complex array!");

		if (info.ndim != 2)
			throw std::runtime_error("Incompatible buffer dimension!");

		int res_x = int(info.shape[0]);
		int res_y = int(info.shape[1]);

		//auto strides = Strides(
		//	info.strides[rowMajor ? 0 : 1] / (py::ssize_t)sizeof(Scalar),
		//	info.strides[rowMajor ? 1 : 0] / (py::ssize_t)sizeof(Scalar));

		//auto map = Eigen::Map<Matrix, 0, Strides>( static_cast<Scalar *>(info.ptr), info.shape[0], info.shape[1], strides);

		auto data = static_cast<std::complex<double> *>(info.ptr);

		auto offset_data = static_cast<int*>(offset_info.ptr);

		new (&m) VoxelGrid( data, domain, V2i(offset_data[0], offset_data[1]) );
	})
	;

	// Constant ============================================================
	py::class_<Constant, Constant::Ptr> class_constant(m, "Constant", class_field);
	class_constant
	.def("__init__",
	[](Constant &m, std::complex<double> value)
	{
		new (&m) Constant(value);
	});

	// RadianceField ============================================================
	py::class_<RadianceField, RadianceField::Ptr> class_radiancefield(m, "RadianceField");
	class_radiancefield
	.def("__call__",[](RadianceField &m, py::array pWS_array, py::array omega_array)
	{
		return m.eval(to_P2d(pWS_array), to_V3d(omega_array));
	})
	.def("dx", [](RadianceField &m, py::array pWS_array, py::array omega_array)
	{
		return m.dx(to_P2d(pWS_array), to_V3d(omega_array));
	})
	.def("dxdx", [](RadianceField &m, py::array pWS_array, py::array omega_array)
	{
		return m.dxdx(to_P2d(pWS_array), to_V3d(omega_array));
	})
	.def("dxdy", [](RadianceField &m, py::array pWS_array, py::array omega_array)
	{
		return m.dxdy(to_P2d(pWS_array), to_V3d(omega_array));
	})
	.def("dy", [](RadianceField &m, py::array pWS_array, py::array omega_array)
	{
		return m.dy(to_P2d(pWS_array), to_V3d(omega_array));
	})
	.def("dydy", [](RadianceField &m, py::array pWS_array, py::array omega_array)
	{
		return m.dydy(to_P2d(pWS_array), to_V3d(omega_array));
	})
	.def("dydx", [](RadianceField &m, py::array pWS_array, py::array omega_array)
	{
		return m.dydx(to_P2d(pWS_array), to_V3d(omega_array));
	})
	.def("dz", [](RadianceField &m, py::array pWS_array, py::array omega_array)
	{
		return m.dz(to_P2d(pWS_array), to_V3d(omega_array));
	})
	;

	// SHEXP ============================================================
	py::class_<SHEXP, SHEXP::Ptr> class_shexp(m, "SHEXP", class_radiancefield);
	class_shexp
	.def("__init__",
	[](SHEXP &m, int order)
	{
		new (&m) SHEXP(order);
	})
	.def("set_coeff", &SHEXP::set_coeff)
	;
}
