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




PYBIND11_MODULE(pnsolver, m)
{
	// PNSystem ==============================
	py::class_<PNSystem> class_pnsystem(m, "PNSystem");
	class_pnsystem
	.def("__init__",
	[](PNSystem &m, const std::string& stencil_name, Domain &domain)
	{
		PNSystem::Stencil stencil;

		if( stencil_name == "noop" )
			stencil = PNSystem::Stencil(0, [](PNSystem::VoxelSystem&, PNSystem::Fields&){});
		else
			stencil = PNSystem::findStencil(stencil_name);
		new (&m) PNSystem(stencil, domain);
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
	.def("getLM",
	[](PNSystem &sys, int sh_index)
	{
		int l,m;
		sys.getLM(sh_index, l, m);
		return py::make_tuple(l, m);
	})
	.def("getNumCoefficients", &PNSystem::getNumCoefficients )
	.def("getNumVoxels", &PNSystem::getNumVoxels )
	.def("getOrder", &PNSystem::getOrder )
	.def("build", &PNSystem::build )
	.def("solve", &PNSystem::solve )
	.def("solveWithGuess", &PNSystem::solveWithGuess )
	.def("solve_boundary", &PNSystem::solve_boundary )
	.def("setField", &PNSystem::setField )
	//.def("solve2", &PNSystem::solve2 )
	//.def("get_A_complex", &PNSystem::get_A_complex )
	//.def("get_b_complex", &PNSystem::get_b_complex )
	.def("get_A_real", &PNSystem::get_A_real )
	.def("get_b_real", &PNSystem::get_b_real )
	.def("setDebugVoxel", &PNSystem::setDebugVoxel )
	.def("get_boundary_A_real", &PNSystem::get_boundary_A_real )
	.def("get_boundary_b_real", &PNSystem::get_boundary_b_real )
	.def("computeGroundtruth", &PNSystem::computeGroundtruth )
	;

	// Domain ============================================================
	py::class_<Domain> class_domain(m, "Domain");
	class_domain
	.def("__init__",
	[](Domain &m, const V2d& size, const V2i& resolution, const V2d& offset)
	{
		new (&m) Domain( size,
						 resolution,
						 offset);
	})
	.def("resolution", &Domain::resolution )
	.def("setResolution", &Domain::setResolution )
	.def("voxelSize", &Domain::voxelSize )
	.def("numVoxels", &Domain::numVoxels )
	.def("worldToVoxel", &Domain::worldToVoxel )
	.def("voxelToWorld", &Domain::voxelToWorld )
	.def("getBoundMax", &Domain::getBoundMax )
	.def("getBoundMin", &Domain::getBoundMin )
	;


	// Field ============================================================
	py::class_<Field, Field::Ptr> class_field(m, "Field");
	class_field
	.def("__call__", &Field::eval)
	;

	// VoxelGrid ============================================================
	py::class_<VoxelGrid, VoxelGrid::Ptr> class_VoxelGrid(m, "VoxelGrid", class_field);
	class_VoxelGrid
	.def("__init__",
	[](VoxelGrid &m, py::array b, Domain& domain, const V2d& offset)
	{
		typedef Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic> Strides;

		// Request a buffer descriptor from Python
		py::buffer_info info = b.request();

		// Some sanity checks ...
		if (info.format != py::format_descriptor<std::complex<double>>::format())
			throw std::runtime_error("Incompatible format: expected a complex array!");

		if (info.ndim != 2)
			throw std::runtime_error("Incompatible buffer dimension!");

		int res_x = int(info.shape[0]);
		int res_y = int(info.shape[1]);

		auto data = static_cast<std::complex<double> *>(info.ptr);
		new (&m) VoxelGrid( data, domain, offset );
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
	//TODO: this produces an compiler error with Eigen under MSVS2017: FLOATING_POINT_ARGUMENT_PASSED__INTEGER_WAS_EXPECTED
	//.def("__call__", &RadianceField::eval)
	.def("eval",
	 [](RadianceField &m, const P2d& pWS, const V2d& omega)
	 {
		V3d o(omega[0], omega[1], 0.0);
		return m.eval(pWS, o);
	 }
	);

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
	.def("eval2",
		[](SHEXP &m, const P2d& pWS)
		{
			return m.eval2(pWS);
		}
	)
	;
}
