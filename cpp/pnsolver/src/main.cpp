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
#include<field/VoxelGridField.h>
#include<field/Constant.h>
#include<field/SHEXP.h>
#include <PNSystem.h>

namespace py = pybind11;




PYBIND11_MODULE(pnsolver, m)
{

	// PNSystem ==============================
	py::class_<PNSystem> class_pnsystem(m, "PNSystem");
	class_pnsystem
	.def("__init__",
	[](PNSystem &m, const std::string& stencil_name, Domain &domain, bool neumannBC)
	{
		PNSystem::Stencil stencil;

		if( stencil_name == "noop" )
			stencil = PNSystem::Stencil::noop();
		else
			stencil = PNSystem::findStencil(stencil_name);
		new (&m) PNSystem(stencil, domain, neumannBC);
	})
	.def("getNumCoefficients", &PNSystem::getNumCoefficients )
	.def("getNumVoxels", &PNSystem::getNumVoxels )
	.def("getOrder", &PNSystem::getOrder )
	.def("build", &PNSystem::build )
	.def("solve", &PNSystem::solve )
	.def("solve_cg", &PNSystem::solve_cg )
	.def("setField", &PNSystem::setField )
	.def("get_A_real", &PNSystem::get_A_real )
	.def("get_b_real", &PNSystem::get_b_real )
	.def("get_A_real_test", &PNSystem::get_A_real_test )
	.def("get_b_real_test", &PNSystem::get_b_real_test )
	.def("get_solve_convergence", &PNSystem::get_solve_convergence )

	.def("setDebugVoxel", &PNSystem::setDebugVoxel )
	//.def("computeGroundtruth", &PNSystem::computeGroundtruth )
	.def("getVoxelInfo", &PNSystem::getVoxelInfo )
	;

	// Domain ============================================================
	py::class_<Domain> class_domain(m, "Domain");
	class_domain
	.def("__init__",
	[](Domain &m, const Eigen::Matrix<double, 3, 1>& size, const Eigen::Matrix<int, 3, 1>& resolution, const Eigen::Matrix<double, 3, 1>& offset)
	{
		new (&m) Domain( size,
						 resolution,
						 offset);
	})
	.def("getResolution",
		[](Domain &m)
		{
			return m.getResolution();
		})
	.def("setResolution",
		[](Domain &m, const Eigen::Matrix<int, 3, 1>& resolution)
		{
			m.setResolution(resolution);
		})
	.def("getVoxelSize",
		[](Domain &m)
		{
			return m.getVoxelSize();
		})
	.def("numVoxels",
		[](Domain &m)
		{
			return m.numVoxels();
		})
	.def("worldToVoxel",
		[](Domain &m, const Eigen::Matrix<double, 3, 1>& pWS)
		{
			return m.worldToVoxel(pWS);
		})
	.def("voxelToWorld",
		[](Domain &m, const Eigen::Matrix<double, 3, 1>& pVS)
		{
			return m.voxelToWorld(pVS);
		})
	.def("getBoundMax",
		[](Domain &m)
		{
			return m.getBoundMax();
		})
	.def("getBoundMin",
		[](Domain &m)
		{
			return m.getBoundMin();
		})
	;


	// Field ============================================================
	py::class_<Field, Field::Ptr> class_field(m, "Field");
	class_field
	.def("__call__",
	[](Field &m, const Eigen::Matrix<double, 3, 1>& pWS)
	{
		return m.eval(pWS);
	})
	;

	// VoxelGrid ============================================================
	py::class_<VoxelGridField, VoxelGridField::Ptr> class_VoxelGridField(m, "VoxelGridField", class_field);
	class_VoxelGridField
	.def("__init__",
	[](VoxelGridField &m, py::array b, Domain& domain, const Eigen::Matrix<double, 3, 1>& offset)
	{
		typedef Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic> Strides;

		// Request a buffer descriptor from Python
		py::buffer_info info = b.request();

		// Some sanity checks ...
		if (info.format != py::format_descriptor<std::complex<double>>::format())
			throw std::runtime_error("Incompatible format: expected a complex array!");

		if (info.ndim != 3)
			throw std::runtime_error("Incompatible buffer dimension!");

		int res_x = int(info.shape[0]);
		int res_y = int(info.shape[1]);
		int res_z = int(info.shape[2]);

		auto data = static_cast<std::complex<double> *>(info.ptr);
		new (&m) VoxelGridField( data, domain, offset );
	})
	.def("test", &VoxelGridField::test)
	.def("getData", &VoxelGridField::getData)
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
	.def("__call__",
	[](RadianceField &m, const Eigen::Matrix<double, 3, 1>& pWS, const Eigen::Matrix<double, 3, 1>& omega)
	{
		return m.eval(pWS, omega);
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
	.def("setCoefficientField", &SHEXP::setCoefficientField)
	.def("getCoefficientField", &SHEXP::getCoefficientField)
	;
}
