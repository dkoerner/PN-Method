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
#include<field/CoefficientGrid.h>


namespace py = pybind11;








PYBIND11_MODULE(pnbuilder_cpp, m)
{
	py::class_<Domain>(m, "Domain", py::buffer_protocol())
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


			//std::cout << voxel_info.format << std::endl;
			//std::cout << py::format_descriptor<int32_t>::format() << std::endl;
			//std::cout << py::format_descriptor<std::int32_t>::format() << std::endl;

			//if (voxel_info.format != py::format_descriptor<int>::format())
			//	throw std::runtime_error("Incompatible format: expected a int array!");

			//if (offset_info.format != py::format_descriptor<int>::format())
			//	throw std::runtime_error("Incompatible format: expected a int array!");

			auto size_data = static_cast<double *>(size_info.ptr);
			auto resolution_data = static_cast<int *>(resolution_info.ptr);
			auto offset_data = static_cast<double *>(offset_info.ptr);

			//std::cout << voxel_data[0] << std::endl;

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

	py::class_<GridLocation>(m, "GridLocation", py::buffer_protocol())
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


			//std::cout << voxel_info.format << std::endl;
			//std::cout << py::format_descriptor<int32_t>::format() << std::endl;
			//std::cout << py::format_descriptor<std::int32_t>::format() << std::endl;

			//if (voxel_info.format != py::format_descriptor<int>::format())
			//	throw std::runtime_error("Incompatible format: expected a int array!");

			//if (offset_info.format != py::format_descriptor<int>::format())
			//	throw std::runtime_error("Incompatible format: expected a int array!");

			auto voxel_data = static_cast<int *>(voxel_info.ptr);
			auto offset_data = static_cast<int *>(offset_info.ptr);

			//std::cout << voxel_data[0] << std::endl;

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
		.def("getShiftedLocation",
		[](GridLocation &m, py::array offset_array)
		{
			auto offset_data = static_cast<int *>(offset_array.request().ptr);
			//std::cout << "offset_data: " << offset_data[0] << " " << offset_data[1]  << std::endl;
			return m.getShiftedLocation(V2i(offset_data[0], offset_data[1]));
		})
		.def("test",
		[](GridLocation &m, int a, int b)
		{
			//auto result = std::div(a, b);
			//int quot = py_div(a,b);
			//int rem = py_mod(a,b);
			//int rem = ((a % b) + b) % b;

			//std::cout << "div(" << a << ", " << b  << ") = " << quot << " " << rem << std::endl;

			//std::cout << a/b << std::endl;
			//std::cout << a%b << std::endl;
			//std::cout << c << std::endl;
		});

	py::class_<CoefficientGrid>(m, "CoefficientGrid", py::buffer_protocol())
		.def("__init__",
		[](CoefficientGrid &m, py::array b, Domain& domain, py::array offset_array)
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

			new (&m) CoefficientGrid( data, domain, V2i(offset_data[0], offset_data[1]) );
		})
		.def("__call__",[](CoefficientGrid &m, GridLocation &l)
		{
			std::cout << "-----------------------------\n";
			auto a =m.eval(l);
			std::cout << "eval=" << a << std::endl;
			return a;
		})
		.def("dx", [](CoefficientGrid &m, GridLocation &l){return m.dx(l);})
		.def("dxdx", [](CoefficientGrid &m, GridLocation &l){return m.dxdx(l);})
		.def("dxdy", [](CoefficientGrid &m, GridLocation &l){return m.dxdy(l);})
		.def("dy", [](CoefficientGrid &m, GridLocation &l){return m.dy(l);})
		.def("dydy", [](CoefficientGrid &m, GridLocation &l){return m.dydy(l);})
		.def("dydx", [](CoefficientGrid &m, GridLocation &l){return m.dydx(l);})
		.def("dz", [](CoefficientGrid &m, GridLocation &l){return m.dz(l);})
		.def("test", [](CoefficientGrid &m){m.test();})
		;
}
