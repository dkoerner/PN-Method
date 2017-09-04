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
#include<field/VoxelGrid.h>
#include<field/Constant.h>
#include<field/SHEXP.h>
#include <PNSystem.h>

namespace py = pybind11;

















V2d to_V2d( py::array array )
{
	auto data = static_cast<double *>(array.request().ptr);
	return V2d(data[0], data[1]);
}


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
	// PNSystem ==============================
	py::class_<PNSystem> class_pnsystem(m, "PNSystem");
	class_pnsystem
	.def("__init__",
	 [](PNSystem &m, Domain &domain, int order)
	 {
		new (&m) PNSystem(domain, order);
	 })
	.def("getGlobalIndex", &PNSystem::getGlobalIndex )
	.def("getVoxelAndCoefficient",
	[](PNSystem &m, int global_index)
	{
		V2i voxel;
		int coeff;
		m.getVoxelAndCoefficient(global_index, voxel, coeff);
		return py::make_tuple(voxel, coeff);
	})
	.def("getNumCoefficients", &PNSystem::getNumCoefficients )
	.def("getNumVoxels", &PNSystem::getNumVoxels )
	.def("build", &PNSystem::build )
	.def("solve", &PNSystem::solve )
	.def("setField", &PNSystem::setField )
	.def("get_A", &PNSystem::get_A )
	.def("get_b", &PNSystem::get_b )
	.def("get_A_real", &PNSystem::get_A_real )
	.def("get_b_real", &PNSystem::get_b_real )
	.def("get_S", &PNSystem::get_S )
	.def("get_S_inv", &PNSystem::get_S_inv )
	.def("get_M", &PNSystem::get_M )
	.def("setDebugVoxel", &PNSystem::setDebugVoxel )
	;

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
	/*
	.def("resolution",
	[](Domain &m)
	{
		V2i resolution = m.resolution();
		py::array_t<int> a({ 2 });
		a.mutable_data()[0] = resolution[0];
		a.mutable_data()[1] = resolution[1];
		return a;
	})
	*/
	.def("resolution", &Domain::resolution )
	.def("setResolution", &Domain::setResolution )
	/*
	.def("voxelSize",
	[](Domain &m)
	{
		V2d voxelSize = m.voxelSize();
		py::array_t<double> a({ 2 });
		a.mutable_data()[0] = voxelSize[0];
		a.mutable_data()[1] = voxelSize[1];
		return a;
	})
	*/
	.def("voxelSize", &Domain::voxelSize )
	.def("numVoxels", &Domain::numVoxels )
	/*
	.def("numVoxels",
	[](Domain &m)
	{
		return m.numVoxels();
	})
	*/
	/*
	.def("worldToVoxel",
	[](Domain &m, py::array pWS_array)
	{
		P2d pVS = m.worldToVoxel( to_P2d(pWS_array) );
		py::array_t<double> a({ 2 });
		a.mutable_data()[0] = pVS[0];
		a.mutable_data()[1] = pVS[1];
		return a;
	})
	*/
	.def("worldToVoxel", &Domain::worldToVoxel )
	.def("voxelToWorld", &Domain::voxelToWorld )
	;


	// Field ============================================================
	py::class_<Field, Field::Ptr> class_field(m, "Field");
	class_field
	.def("__call__",[](Field &m, py::array pWS_array){return m.eval(to_P2d(pWS_array));})
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

		// Some sanity checks ...
		//std::cout << "format=" << info.format << std::endl;
		if (info.format != py::format_descriptor<std::complex<double>>::format())
			throw std::runtime_error("Incompatible format: expected a complex array!");

		if (info.ndim != 2)
			throw std::runtime_error("Incompatible buffer dimension!");

		int res_x = int(info.shape[0]);
		int res_y = int(info.shape[1]);

		auto data = static_cast<std::complex<double> *>(info.ptr);
		new (&m) VoxelGrid( data, domain, to_V2d(offset_array) );
	})
	.def("test", &VoxelGrid::test )
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
	.def("get_coeff", &SHEXP::get_coeff)
	;
}
