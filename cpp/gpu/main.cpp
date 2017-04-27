#include <iostream>
#include <string>

#include <math/frame.h>

#include <util/string.h>
#include <util/envmap.h>
#include <util/timer.h>
#include <util/field.h>
#include <util/moexp.h>
#include <util/moexp2d.h>
#include <util/wedge.h>

#include <volume.h>
#include <shcache.h>

#include <cuda.h>
#include <builtin_types.h>
//#include <drvapi_error_string.h>

#include <cuda_runtime.h>
#include <cuda/CudaData.h>
#include <cuda/CudaEnvMap.h>
#include <cuda/CudaPixelBuffer.h>
#include <cuda/CudaVolume.h>
#include <cuda/CudaLight.h>


// Forward declare the function in the .cu file
void vectorAddition(const float* a, const float* b, float* c, int n);
void pathtrace_envmap( cumath::V3f pWS, CudaPtr<CudaEnvMap>* envmap, CudaPtr<CudaVolume>* volume, CudaPtr<CudaLight>* light, int numSamples );
void compute_sh_coefficients( CudaPtr<CudaPixelBuffer>* coords_ptr,
							  CudaPtr<CudaPixelBuffer>* f_ptr,
							  cumath::V3f* coeffs_host,
							  int order,
							  int numCoeffs);


void getCoordFromIndex( const V3i& resolution, int index, int& x, int& y, int& z )
{
	div_t divresult;

	divresult = div( index, resolution.x()*resolution.y() );

	z = divresult.quot;
	divresult = div( divresult.rem, resolution.x() );
	y = divresult.quot;
	x = divresult.rem;
}


void pathtrace_fluencefield( CudaPtr<CudaVoxelGrid<float>>* fluencegrid, CudaPtr<CudaVolume>* volume, CudaPtr<CudaLight>* light, int numSamples );
void generate_fluencefield( const std::string& filename, Transformd volume_localToWorld, CudaPtr<CudaVolume> cuvolume_ptr, CudaPtr<CudaLight> culight_ptr)
{
	int numSamples = 1000;
	V3i resolution(64, 64, 64);

	CudaVoxelGrid<float> cufluencegrid;
	CudaPtr<CudaVoxelGrid<float>> cufluencegrid_ptr(&cufluencegrid);
	cufluencegrid.initialize( cumath::V3i(resolution.x(), resolution.y(), resolution.z()));


	pathtrace_fluencefield(&cufluencegrid_ptr, &cuvolume_ptr, &culight_ptr, numSamples);

	std::vector<float> data( resolution.x()*resolution.y()*resolution.z() );
	cufluencegrid_ptr.host()->download( (unsigned char*)data.data() );
	field::VoxelGridField<float>::Ptr fluencegrid = std::make_shared<field::VoxelGridField<float>>(resolution, data.data());

	Fieldf::Ptr fluencefield = field::xform<float>(fluencegrid, volume_localToWorld);
	field::write( filename, fluencefield, resolution, volume_localToWorld );
}

void pathtrace_fluencefield_2d( CudaPtr<CudaPixelBuffer>* fluencegrid, CudaPtr<CudaVolume2D>* volume, int numSamples );


void experiment_pt_2d_anisotropic();
void experiment_nebulae();
void experiment_moexp2d();
void experiment_sh_method();
int main()
{
	int deviceCount = 0;
	int cudaDevice = 0;
	char cudaDeviceName [100];

	cuInit(0);
	cuDeviceGetCount(&deviceCount);
	cuDeviceGet(&cudaDevice, 0);
	cuDeviceGetName(cudaDeviceName, 100, cudaDevice);
	std::cout << "Number of devices: " << deviceCount << std::endl;
	std::cout << "Device name:" << cudaDeviceName << std::endl;





	//experiment_pt_2d_anisotropic();
	//experiment_nebulae();
	//experiment_moexp2d();
	experiment_sh_method();






	return 0;
}


// resolution.x -> resolution along theta angle
// resolution.y -> resolution along phi angle
template<typename T>
T integrate_sphere( moexp::SphericalFunction<T> f, V2i resolution )
{
	T result = 0.0;

	double min_theta = 0;
	double max_theta = M_PI;
	double min_phi = 0;
	double max_phi = 2.0*M_PI;

	int resolution_theta=resolution.x(); // resolution along theta angle
	int resolution_phi=resolution.y(); // resolution along phi angle
	double dtheta = (max_theta-min_theta)/double(resolution_theta);
	double dphi = (max_phi-min_phi)/double(resolution_phi);

	double pixel_area = dtheta*dphi;

	for( int t=0;t<resolution_theta;++t )
		for( int p=0;p<resolution_phi;++p )
		{
			double phi = dphi*p;
			double theta = dtheta*t;
			result+=f(theta, phi)*pixel_area*std::sin(theta);
		}

	return result;
}



// res gives resolution in theta and phi
void compute_sh_coeffs( moexp::SphericalFunction<double> f, V2i res, int order, std::vector<moexp::complex>& coeffs )
{
	int numCoeffs = moexp::numSHCoefficients(order);
	coeffs.resize(numCoeffs, 0.0);

	double min_theta = 0;
	double max_theta = M_PI;
	double min_phi = 0;
	double max_phi = 2.0*M_PI;

	//int resolution_theta=map.bitmap().rows(); // resolution along theta angle
	//int resolution_phi=map.bitmap().cols(); // resolution along phi angle
	int resolution_theta=res.x(); // resolution along theta angle
	int resolution_phi=res.y(); // resolution along phi angle
	double dtheta = (max_theta-min_theta)/double(resolution_theta);
	double dphi = (max_phi-min_phi)/double(resolution_phi);


	// cpu version for computing moments ---
	double pixel_area = dtheta*dphi;

	for( int t=0;t<resolution_theta;++t )
		for( int p=0;p<resolution_phi;++p )
		{
			double phi = dphi*p;
			double theta = dtheta*t;

			moexp::complex f_sample = f(theta, phi);

			for( int l=0;l<=order;++l )
				for( int m=-l;m<=l;++m )
				{
					coeffs[moexp::shIndex(l, m)]+=f_sample*moexp::Y(l, m, theta, phi)*pixel_area*std::sin(theta);
				}
		}
}


void functionToBgeo( const std::string& filename, moexp::_2d::Function<double> fun, bool polar = false )
{
	houio::Geometry::Ptr geo = houio::Geometry::createLineGeometry();
	houio::Attribute::Ptr pattr = geo->getAttr("P");

	int numSamples = 1000;
	double period = 2.0*M_PI;

	double dt = period/double(numSamples);
	for( int i=0;i<numSamples;++i )
	{
		double t = i*dt;
		double f = fun(t);

		houio::math::V3f p(0.0f, 0.0f, 0.0f);
		if(polar)
			p = houio::math::V3f(std::cos(t), std::sin(t), 0.0f)*float(f);
		else
			p = houio::math::V3f(t, f, 0.0f);


		pattr->appendElement<houio::math::V3f>(p);
		if(i>0)
			geo->addLine(i-1, i);
	}

	houio::HouGeoIO::xport(filename, geo);
}


void experiment_sh_method()
{
	//double exposure = 0.0;
	V2i res( 512, 1024 );


	// radiance field ---
	EnvMap L_baked("test_0_filtered.exr");
	auto L = [&]( double theta, double phi )->Color3f
	{
		return L_baked.eval(theta, phi);
	};
	rasterizeSphericalFunctionSphere("shm/fun_L.bgeo", L, 8.0);

	// phase function ---
	PhaseFunction::Ptr f = volumes::hg(0.9);
	/*
	functionToBgeo("shm/fun_f.bgeo", [&](double theta)->double
	{
		V3d wo = sphericalDirection(0.0, 0.0);
		V3d wi = sphericalDirection(theta-M_PI, 0.0);
		return f->eval(wi, wo).r();
	}, true);
	*/
	rasterizeSphericalFunctionSphere("shm/fun_f.bgeo", [&](double theta, double phi)->double
	{
		V3d wo(0.0, 0.0, 1.0);
		V3d wi = sphericalDirection(theta, phi);
		return f->eval(wi, wo).r();
	}, 0.0, 20);


	// scattering operator as convolution ---
	auto S = [&]( double theta, double phi )->Color3f
	{
		V3d wo = sphericalDirection(theta, phi);
		auto L_times_phase = [&](double theta, double phi)->Color3f
		{
			V3d wi = sphericalDirection(theta, phi);
			return L( theta, phi )*f->eval(wi, wo);
		};

		return integrate_sphere<Color3f>(L_times_phase, res);
	};
	//rasterizeSphericalFunctionSphere("shm/fun_S.bgeo", S, 8.0, 20);

	// scattering operator as inner product integral of rotated phase function ---
	auto S_ip = [&]( double theta, double phi )->Color3f
	{
		V3d d = sphericalDirection(theta, phi);
		Eigen::Quaterniond rho = Eigen::Quaterniond::FromTwoVectors( V3d(0.0, 0.0, 1.0), d );
		//Eigen::Quaterniond rho_inverse = Eigen::Quaterniond::FromTwoVectors( d, V3d(0.0, 0.0, 1.0) );
		Eigen::Quaterniond rho_inverse = Eigen::Quaterniond::FromTwoVectors( d, d );

		auto ip = [&](double theta, double phi)->Color3f
		{
			V3d wi = sphericalDirection(theta, phi);
			V3d wo = d;

			// apply inverse rotation to both arguments of the phase function
			// this is identical to forward rotation of the phase function
			V3d wi_inverse_rotation = rho_inverse*wi;
			V3d wo_inverse_rotation = rho_inverse*wo;
			return L( theta, phi )*f->eval(wi_inverse_rotation, wo_inverse_rotation);
		};

		return integrate_sphere<Color3f>(ip, res);
	};
	//rasterizeSphericalFunctionSphere("shm/fun_S_ip.bgeo", S_ip, 8.0, 20);

	// eigenfunctions of scattering operator ---
	// we want to see that the inner product integral of Ylm with a rotationally symmetric function, rotated to \omega
	// is identical to a constant factor times Ylm(\omega)
	for( int l = 2;l<0;++l )
	{
		//for( int m = -l;m<=l;++m )
		for( int m = 0;m<=l;++m )
		{
			auto Ylm = [&](double theta, double phi)->Color3f
			{
				moexp::complex v = moexp::Y( l, m, theta, phi );
				return Color3f( std::abs(v.real()), std::abs(v.imag()), 0.0f );
			};

			auto lambdal_times_Ylm = [&](double theta, double phi)->Color3f
			{
				double lambda_l = std::sqrt(4.0*M_PI / (2.0*l+1.0));
				moexp::complex v = lambda_l*moexp::Y( l, m, theta, phi );
				return Color3f( std::abs(v.real()), std::abs(v.imag()), 0.0f );
			};

			auto Cl = [&]( double theta, double phi )->Color3f
			{
				V3d d = sphericalDirection(theta, phi);

				// build rotation, which rotates the pole axis of Yl0 (0,0,1) to \omega
				Eigen::Quaterniond rho = Eigen::Quaterniond::FromTwoVectors( V3d(0.0, 0.0, 1.0), d );
				Eigen::Quaterniond rho_inverse = Eigen::Quaterniond::FromTwoVectors( d, V3d(0.0, 0.0, 1.0) );

				auto ip_integrand = [&]( double theta, double phi )->moexp::complex
				{
					moexp::complex ylm = moexp::Y( l, m, theta, phi );

					V3d d2 = sphericalDirection(theta, phi);
					V3d d2_inverse_rotation = rho_inverse*d2;
					P2d coords_rotated_theta_phi = sphericalCoordinates( d2_inverse_rotation );

					// apply rotation
					moexp::complex yl0 = moexp::Y( l, 0, coords_rotated_theta_phi.x(), coords_rotated_theta_phi.y() );

					return ylm*yl0;
				};

				moexp::complex v = integrate_sphere<moexp::complex>( ip_integrand, res);
				return Color3f( std::abs(v.real()), std::abs(v.imag()), 0.0f );
			};
			rasterizeSphericalFunctionSphere("shm/fun_Y_" + toString<int>(l) + "_" + toString<int>(m) + ".bgeo", Ylm, 0.0, 120);
			rasterizeSphericalFunctionSphere("shm/fun_lY_" + toString<int>(l) + "_" + toString<int>(m) + ".bgeo", lambdal_times_Ylm, 0.0, 120);
			rasterizeSphericalFunctionSphere("shm/fun_Cl_Y_" + toString<int>(l) + "_" + toString<int>(m) + ".bgeo", Cl, 0.0, 20);
		}//m
	}//l

	// convolution as product of sh coefficients ---
	int maxOrder = 10;
	for( int order = 0;order<=maxOrder;++order )
	{
		// we want to see, that the scattering operator can be expressed as a
		// product of the sh coefficients of radiance field and phase function

		// sh coefficients of the phase function
		std::vector<moexp::complex> f_lm;
		auto f_local = [&]( double theta, double phi )->double
		{
			V3d wi(0.0, 0.0, 1.0);
			V3d wo = sphericalDirection(theta, phi);
			return f->eval(wi, wo).r();
		};
		compute_sh_coeffs(f_local, V2i(256, 512), order, f_lm);

		// sh coefficients of the radiance field
		std::vector<moexp::complex> L_lm;
		auto L_red = [&]( double theta, double phi )->double
		{
			// there is some problem with matching the coordinates between envmap and spherical harmonics projection
			// the problem seems to something related to a mismatch between the phi coordinates
			// we fix it by simply flipping the y axis
			V3d w = sphericalDirection(theta, phi);
			V3d w2(w.x(), -w.y(), w.z());
			P2d theta_phi = sphericalCoordinates(w2);
			return L_baked.eval(theta_phi.x(), theta_phi.y()).r();

			//return L_baked.eval(theta, phi).r();
		};
		compute_sh_coeffs(L_red, V2i(256, 512), order, L_lm);

		// coefficient product
		auto rhs = [&](double theta, double phi)->Color3f
		{
			moexp::complex v = 0.0;
			for( int l=0;l<=order;++l )
			{
				for( int m = -l;m<=l;++m )
				{
					moexp::complex ylm = moexp::Y(l, m, theta, phi);
					double lambda_l = std::sqrt(4.0*M_PI / (2.0*l+1.0));
					v+= f_lm[moexp::shIndex(l, 0)] * L_lm[moexp::shIndex(l, m)] * lambda_l * ylm;
					//v+= L_lm[moexp::shIndex(l, m)]* ylm;
					//v+= f_lm[moexp::shIndex(l, m)]* ylm;
				}
			}

			//return Color3f( std::abs(v.real()), std::abs(v.imag()), 0.0f );
			return Color3f( v.real() );
		};
		rasterizeSphericalFunctionSphere("shm/fun_rhs_" + toString<int>(order) + ".bgeo", rhs, 8.0, 20);
	}











	return;
}

// ==============================================================



M22d build_SGGX2D_matrix( const V2d& w1, double pd_x, double pd_y  );

double sggx3d_sigma(V3d wi, double S_xx, double S_yy, double S_zz, double S_xy, double S_xz, double S_yz)
{
	const float sigma_squared = wi.x()*wi.x()*S_xx + wi.y()*wi.y()*S_yy + wi.z()*wi.z()*S_zz + 2.0f * (wi.x()*wi.y()*S_xy + wi.x()*wi.z()*S_xz + wi.y()*wi.z()*S_yz);
	return (sigma_squared > 0.0f) ? sqrtf(sigma_squared) : 0.0f; // conditional to avoid numerical errors
}

double interpolate( double val, double y0, double x0, double y1, double x1 ) {
	return (val-x0)*(y1-y0)/(x1-x0) + y0;
}

double base( double val ) {
	if ( val <= -0.75 ) return 0;
	else if ( val <= -0.25 ) return interpolate( val, 0.0, -0.75, 1.0, -0.25 );
	else if ( val <= 0.25 ) return 1.0;
	else if ( val <= 0.75 ) return interpolate( val, 1.0, 0.25, 0.0, 0.75 );
	else return 0.0;
}

double red( double gray ) {
	return base( gray - 0.5 );
}
double green( double gray ) {
	return base( gray );
}
double blue( double gray ) {
	return base( gray + 0.5 );
}



struct SGGX3D
{
	SGGX3D( Framed frame, V3d projectedAreas )
	{
		V3d w1 = frame.s;
		V3d w2 = frame.t;
		V3d w3 = frame.n;

		double S_11 = projectedAreas.x()*projectedAreas.x();
		double S_22 = projectedAreas.y()*projectedAreas.y();
		double S_33 = projectedAreas.z()*projectedAreas.z();


		M33d m;
		m << w1.x(), w2.x(), w3.x(),
			 w1.y(), w2.y(), w3.y(),
			 w1.z(), w2.z(), w3.z();

		M33d s = V3d(S_11, S_22, S_33).asDiagonal();

		m_S = m*s*m.transpose();


		// build transform ---
		Eigen::Affine3d unitSphereToEllipsoid;

		double sx = std::pow(S_22*S_33/S_11, 1.0/4)/std::sqrt(M_PI);
		double sy = std::pow(S_11*S_33/S_22, 1.0/4)/std::sqrt(M_PI);
		double sz = std::pow(S_11*S_22/S_33, 1.0/4)/std::sqrt(M_PI);
		unitSphereToEllipsoid = Eigen::Affine3d(m)*Eigen::Scaling(V3d(sx, sy, sz));
		m_unitSphereToEllipsoid = Transformd(unitSphereToEllipsoid);
	}

	double projectedArea(const V3d& d)const
	{
		//return std::sqrt( d.dot(m_S*d) );
		double S_xx = m_S.coeffRef(0,0);
		double S_xy = m_S.coeffRef(0,1);
		double S_xz = m_S.coeffRef(0,2);
		double S_yy = m_S.coeffRef(1,1);
		double S_yz = m_S.coeffRef(1,2);
		double S_zz = m_S.coeffRef(2,2);
		const float sigma_squared = d.x()*d.x()*S_xx + d.y()*d.y()*S_yy + d.z()*d.z()*S_zz + 2.0f * (d.x()*d.y()*S_xy + d.x()*d.z()*S_xz + d.y()*d.z()*S_yz);
		return (sigma_squared > 0.0f) ? sqrtf(sigma_squared) : 0.0f; // conditional to avoid numerical errors
	}


	M33d m_S;
	Transformd m_unitSphereToEllipsoid;
};



struct SGGX2D
{
	SGGX2D( const V2d& frame_s, const V2d& projectedDistances )
	{
		V2d w1 = frame_s;
		V2d w2 = V2d(-frame_s.y(), frame_s.x());

		double S_11 = projectedDistances.x()*projectedDistances.x();
		double S_22 = projectedDistances.y()*projectedDistances.y();


		M22d m;
		m << w1.x(), w2.x(),
			 w1.y(), w2.y();

		M22d s = V2d(S_11, S_22).asDiagonal();

		m_S = m*s*m.transpose();


		// build transform ---
		Eigen::Affine2d unitCircleToEllipsoid;

		//double sx = std::pow(S_22*S_33/S_11, 1.0/4)/std::sqrt(M_PI);
		//double sy = std::pow(S_11*S_33/S_22, 1.0/4)/std::sqrt(M_PI);
		double sx = S_11;
		double sy = S_22;
		unitCircleToEllipsoid = Eigen::Affine2d(m)*Eigen::Scaling(V2d(sx, sy));
		m_unitCircleToEllipsoid = Transform2Dd(unitCircleToEllipsoid);
	}

	double projectedDistance(const V2d& d)const
	{
		//return std::sqrt(d.dot(m_S*d));

		///*
		double S_xx = m_S.coeffRef(0,0);
		double S_xy = m_S.coeffRef(0,1);
		double S_yy = m_S.coeffRef(1,1);
		float sigma_squared = S_xx*d.x()*d.x() + S_yy*d.y()*d.y() + 2.0*S_xy*d.x()*d.y();
		if(sigma_squared>0.0)
			return sqrt(sigma_squared);
		// conditional to avoid numerical errors
		return 0.0;
		//*/
	}


	M22d m_S;
	Transform2Dd m_unitCircleToEllipsoid;
};


void experiment_moexp2d()
{

	{
		Wedge wedge;
		wedge.addParm( "d", linearSamples<float>(0.0f, M_PI*2.0f, 20) );
		wedge.addParm( "pd", linearSamples<float>(0.0f, 1.0f, 10) );
		std::vector<int> moments;
		for( int i=0;i<=10;++i )
		//for( int i=0;i<=4;++i )
			moments.push_back(i);
		wedge.addParm( "o", moments);

		for( auto it = wedge.begin(), end = wedge.end(); it!=end;++it )
		{
			double d = it.getFloat("d");
			double pd = it.getFloat("pd");
			int m = it.getInt("o");
			//std::cout << "d=" << d << " pd=" << pd << " m=" << m << std::endl;
			//std::cout << it.expand_index("$0 $1 $2") << std::endl;

			bool polar = true;

			V2d frame_s = V2d(std::cos(d), std::sin(d)).normalized();
			V2d projectedDistances(pd, 1.0);
			//SGGX2D sggx(frame_s, projectedDistances);

			CudaSGGX2D cusggx2d(cumath::V2f(frame_s.x(), frame_s.y()),
								cumath::V2f(projectedDistances.x(), projectedDistances.y()));

			moexp::_2d::Function<double> fun = [&]( double t )
			{
				// 2d sggx
				cumath::V2f d( std::cos(t), std::sin(t) );
				return cusggx2d.projectedDistance(d);

				//V2d d( std::cos(t), std::sin(t) );
				//return sggx.projectedDistance(d);
			};

			std::unique_ptr<std::vector<double>> coeffs = moexp::_2d::projectFourier(m, fun);

			moexp::_2d::Function<double> fun_rec = [&]( double t )
			{
				return moexp::_2d::fourierSum(m, (*coeffs).data(), t);
			};

			// write groundtruth function
			functionToBgeo( it.expand_index("moexp2d/sggx2d/groundtruth_$0_$1_$2.bgeo"), fun, polar );
			functionToBgeo( it.expand_index("moexp2d/sggx2d/reconstruction_$0_$1_$2.bgeo"), fun_rec, polar );
		}

		return;
	}

	/*
	{
		V3d n = V3d(1.0, 1.0, 1.0).normalized();
		//V3d n = V3d(0.0, 0.0, 1.0).normalized();
		V3d projectedArea(1.0, 1.0, 0.001);
		//V3d projectedArea(0.001, 0.001, 1.0);
		double vmin = std::min(std::min(projectedArea.x(), projectedArea.y()), projectedArea.z());
		double vmax = std::max(std::max(projectedArea.x(), projectedArea.y()), projectedArea.z());

		SGGX3D sggx( Framed(n), projectedArea);
		moexp::SphericalFunction<Color3f> func = [&]( double theta, double phi ) -> Color3f
		{
			V3d d  = sphericalDirection(theta, phi);
			double value = sggx.projectedArea(d);
			double t = ((value-vmin)/(vmax-vmin))*2.0-1.0;
			return Color3f(red(t), green(t), blue(t));
		};
		moexp::SphericalFunction<double> func_pa = [&]( double theta, double phi ) -> double
		{
			V3d d  = sphericalDirection(theta, phi);
			return sggx.projectedArea(d);
		};

		std::string filename = "moexp2d/sggx3d.bgeo";
		//filename = replace(filename, "$0", toString<int>(i));
		//rasterizeSphericalFunctionSphere( filename, func, sggx.m_unitSphereToEllipsoid);
		//rasterizeSphericalFunctionSphere( filename, func);
		displaceSphere( filename, func_pa);

		std::cout << "range=" << vmin << " " << vmax << std::endl;

		return;
	}
	*/


	/*
	int order = 30;
	bool polar = true;

	//cumath::V2f frame_s = cumath::normalize(cumath::V2f(1.0f, 1.0f));

	V2d frame_s = V2d(1.0, 0.0).normalized();
	V2d projectedDistances(0.5, 1.0);
	SGGX2D sggx(frame_s, projectedDistances);

	CudaSGGX2D cusggx2d(cumath::V2f(frame_s.x(), frame_s.y()),
						cumath::V2f(projectedDistances.x(), projectedDistances.y()));

	moexp::_2d::Function<double> fun = [&]( double t )
	{
		// 2d sggx
		cumath::V2f d( std::cos(t), std::sin(t) );
		return cusggx2d.projectedDistance(d);

		//V2d d( std::cos(t), std::sin(t) );
		//return sggx.projectedDistance(d);
	};

	std::unique_ptr<std::vector<double>> coeffs = moexp::_2d::projectFourier(order, fun);

	for( int i=0;i<=order;++i )
	{
		moexp::_2d::Function<double> fun_rec = [&]( double t )
		{
			return moexp::_2d::fourierSum(i, (*coeffs).data(), t);
		};
		std::string filename = "moexp2d/reconstruction_$0.bgeo";
		filename = replace(filename, "$0", toString<int>(i));
		functionToBgeo( filename, fun_rec, polar );
	}


	// write groundtruth function
	functionToBgeo( "moexp2d/groundtruth.bgeo", fun, polar );
	*/

	/*
	// moments ------
	std::vector<moexp::_2d::Tensor<double>> moments =  moexp::_2d::convertFourierCoefficientsToMoments((*coeffs).data());

	// validate the moments by comparing against the fourier basis functions
	for( int i=0;i<3;++i )
	{
		moexp::_2d::Function<double> fun_basis = [&]( double t )
		{
			return moexp::_2d::fourierBasis(i, (*coeffs).data(), t);
		};
		moexp::_2d::Function<double> fun_moments = [&]( double t )
		{
			V2d d( std::cos(t), std::sin(t) );
			return moexp::_2d::contract(moments[i], d );
		};

		{
			std::string filename = "moexp2d/basis_fourier_$0.bgeo";
			filename = replace(filename, "$0", toString<int>(i));
			functionToBgeo( filename, fun_basis );
		}
		{
			std::string filename = "moexp2d/basis_moments_$0.bgeo";
			filename = replace(filename, "$0", toString<int>(i));
			functionToBgeo( filename, fun_moments );
		}
	}
	*/
}

// ==============================================================


void experiment_pt_2d_anisotropic()
{

	std::string basepath = "c:/projects/epfl/experiments/python/anisotropic_absorption_2d";
	//std::string basepath = "c:/projects/epfl/experiments/python/isotropic_absorption_2d";

	///*
	{
		//Bitmap map(basepath + "/test.exr");
		Bitmap map;
		/*
		map.loadTXT(basepath + "/solver.txt");
		map.saveEXR(basepath + "/solver.exr");

		map.loadTXT(basepath + "/phi_groundtruth.txt");
		map.saveEXR(basepath + "/phi_groundtruth.exr");

		map.fill(Color3f(0.0f));
		map.loadTXT(basepath + "/gradphi_x_groundtruth.txt", 0);
		map.loadTXT(basepath + "/gradphi_y_groundtruth.txt", 1);
		map.saveEXR(basepath + "/gradphi_groundtruth.exr");

		map.loadTXT(basepath + "/discretization_phi.txt");
		map.saveEXR(basepath + "/discretization_phi.exr");

		map.fill(Color3f(0.0f));
		map.loadTXT(basepath + "/discretization_gradphi_x.txt", 0);
		map.loadTXT(basepath + "/discretization_gradphi_y.txt", 1);
		map.saveEXR(basepath + "/discretization_gradphi.exr");
		*/

		map.loadTXT(basepath + "/discretization_residuum.txt");
		map.saveEXR(basepath + "/discretization_residuum.exr");

	}

	return;
	//*/




	V2i resolution(512, 512);
	int numSamples = 1000;

	// the following variables define the anisotropic medium
	float density = 5.0;
	V2d frame_s = V2d(1.0, 0.0).normalized();
	V2d projectedDistances(0.5, 1.0);
	//V2d projectedDistances(1.0, 1.0);


	Eigen::Affine2d localToWorld;
	localToWorld = Eigen::Translation2d(V2d(-0.5, -0.5));
	Transform2Dd volume_localToWorld(localToWorld);


	// resolution

	CudaPixelBuffer cufluencegrid;
	cufluencegrid.initialize(resolution.x(), resolution.y());

	CudaVolume2D cuvolume;
	cuvolume.m_density = cumath::V3f(density);
	cuvolume.m_sigma_t_max = cuvolume.m_density;
	M33d l2w = volume_localToWorld.getMatrix();
	M33d l2w_inv = volume_localToWorld.getInverseMatrix();
	for( int j=0;j<3;++j )
		for( int i=0;i<3;++i )
		{
			// NB: we transpose the matrix here because cumath uses row vectors
			cuvolume.m_localToWorld.matrix.m[j][i] = float(l2w.coeff(i,j));
			cuvolume.m_localToWorld.inverse.m[j][i] = float(l2w_inv.coeff(i,j));
		}
	cuvolume.m_sggx = CudaSGGX2D( cumath::V2f(frame_s.x(), frame_s.y()),
								  cumath::V2f(projectedDistances.x(), projectedDistances.y()) );

	// kick of render ---
	CudaPtr<CudaPixelBuffer> fluencegrid_ptr(&cufluencegrid);
	CudaPtr<CudaVolume2D> volume_ptr(&cuvolume);
	pathtrace_fluencefield_2d( &fluencegrid_ptr, &volume_ptr, numSamples );

	// write result to disk ---
	Bitmap map(resolution);
	cufluencegrid.download((unsigned char*)map.data());
	//map.saveEXR("pt2da/test.exr");
	map.saveEXR(basepath+"/mc.exr");
	map.saveTXT(basepath+"/mc.txt");
}


void experiment_nebulae()
{
	/*
	{
		EnvMap map("test_0.exr");
		rasterizeSphericalFunctionSphere("test");
		return;
	}
	*/


	// setup volume ---------------

	//
	CudaVoxelGrid<float> cuvoxelgrid;
	CudaPtr<CudaVoxelGrid<float>> cuvoxelgrid_ptr(&cuvoxelgrid);

	//Volume::Ptr volume;
	Box3d volume_boundWS;
	float volume_sigma_t_max;
	Transformd volume_localToWorld;
	{
		std::string basePath = "c:/projects/visus/data";

		// volume ----
		double stepSize;
		field::VoxelGridField<float, float>::Ptr density = field::bgeo<float>(basePath + "/datasets/nebulae200.bgeo", &volume_localToWorld, &stepSize);

		volume_sigma_t_max = field_maximum<float>(density->grid);

		// setup localspace bounding box
		Box3d bboxLS;
		bboxLS.reset();
		bboxLS.min = P3d(0.0f,0.0f,0.0f);
		bboxLS.max = P3d(1.0f,1.0f,1.0f);

		// compute worldspace boundingbox
		volume_boundWS.reset();
		for( int i=0;i<8;++i )
			volume_boundWS.expandBy( volume_localToWorld*bboxLS.getCorner(i) );

		// initialize cudavoxelgrid
		V3i resolution = density->grid.getResolution();

		cuvoxelgrid.initialize( cumath::V3i(resolution.x(), resolution.y(), resolution.z()),
								density->grid.getRawPointer());
	}


	// ...
	CudaVolume cuvolume;
	cuvolume.m_boundWS = cumath::Box3f( cumath::V3f(volume_boundWS.min.x(), volume_boundWS.min.y(), volume_boundWS.min.z()),
										cumath::V3f(volume_boundWS.max.x(), volume_boundWS.max.y(), volume_boundWS.max.z()) );

	cuvolume.m_sigma_t = cuvoxelgrid_ptr.device();
	cuvolume.m_sigma_t_max = cumath::V3f(volume_sigma_t_max);
	M44d l2w = volume_localToWorld.getMatrix();
	M44d l2w_inv = volume_localToWorld.getInverseMatrix();
	for( int j=0;j<4;++j )
		for( int i=0;i<4;++i )
		{
			// NB: we transpose the matrix here because cumath uses row vectors
			cuvolume.m_localToWorld.matrix.m[j][i] = float(l2w.coeff(i,j));
			cuvolume.m_localToWorld.inverse.m[j][i] = float(l2w_inv.coeff(i,j));
		}

	CudaPtr<CudaVolume> cuvolume_ptr(&cuvolume);


	// setup light ----------------
	CudaLight culight;
	CudaPtr<CudaLight> culight_ptr(&culight);
	culight.m_direction = cumath::V3f(0.0, -1.0, 0.0);
	culight.m_radiance = cumath::V3f(1.0, 1.0, 1.0);
	culight.m_distance = volume_boundWS.getExtents().norm();





	generate_fluencefield("nebulae_fluence.bgeo", volume_localToWorld, cuvolume_ptr, culight_ptr);
	//generate_shcache(volume_localToWorld, cuvolume_ptr, culight_ptr);
}

void generate_shcache( Transformd volume_localToWorld, CudaPtr<CudaVolume> cuvolume_ptr, CudaPtr<CudaLight> culight_ptr )
{
	V2i envmap_resolution(128, 64);
	V3i sampling_resolution(64,64,64);
	//V3i sampling_resolution(1);
	int order = 30;
	int numVoxels = sampling_resolution.x()*sampling_resolution.y()*sampling_resolution.z();


	// setup envmap which is being sampled -----------
	CudaEnvMap cumap;
	CudaPtr<CudaEnvMap> cumap_ptr(&cumap);
	cumap.initialize(envmap_resolution.x(), envmap_resolution.y());


	// Buffer for sh coefficients
	CudaPixelBuffer cushcoeffs;
	CudaPtr<CudaPixelBuffer> cushcoeffs_ptr(&cushcoeffs);


	// environment map of the current voxel...
	EnvMap map(envmap_resolution);

	// shcache which holds the final sh coefficients for all voxels...
	SHCache shcache(order, sampling_resolution, volume_localToWorld);


	// go for each voxel...
	Timer timer;
	timer.start();
	for( int voxel=0;voxel<numVoxels;++voxel )
	{
		// clear envmap on device to have a clean plate for the next voxel
		cumap_ptr.host()->clear();

		// find voxel center in world space
		int i, j, k;
		getCoordFromIndex(sampling_resolution, voxel, i, j, k);

		P3d pLS((i+0.5)/sampling_resolution.x(),
				(j+0.5)/sampling_resolution.y(),
				(k+0.5)/sampling_resolution.z());
		P3d pWS = volume_localToWorld*pLS;
		cumath::V3f p(pWS.x(), pWS.y(), pWS.z());

		// kick of rendering
		if(0)
		{
			int numSamples = 80;
			pathtrace_envmap( p, &cumap_ptr, &cuvolume_ptr, &culight_ptr, numSamples );
		}

		// download result and save to disk
		if(0)
		{
			cumap.download((unsigned char*)map.bitmap().data());
			//std::string filename = "nebulae_envmaps/envmap_$0.exr";
			std::string filename = "test_$0.exr";
			filename = replace(filename, "$0", toString(voxel));
			map.bitmap().saveEXR(filename);
		}

		// load precomputed envmap per voxel
		if(0)
		{
			//std::string filename = "test_$0.exr";
			std::string filename = "nebulae_envmaps/envmap_$0.exr";
			filename = replace(filename, "$0", toString(voxel));
			// load rendered images
			map.bitmap() = Bitmap(filename);
		}

		// apply filter
		if(0)
		{
			double filter_width = 3.0;
			map.filter(gaussFilter(filter_width));
		}

		// compute sh coefficients
		if(0)
		{
			int numCoeffs = moexp::numSHCoefficients(order);


			moexp::SphericalFunction<Color3f> fun = [&]( double theta, double phi ) -> Color3f
			{
				return map.eval(theta, phi);
			};

			//std::unique_ptr<std::vector<Color3f>> coeffs_gt = moexp::project_Y_real( order, fun, 150000 );


			double min_theta = 0;
			double max_theta = M_PI;
			double min_phi = 0;
			double max_phi = 2.0*M_PI;

			int resolution_phi=map.bitmap().cols(); // resolution along phi angle
			int resolution_theta=map.bitmap().rows(); // resolution along theta angle
			double dtheta = (max_theta-min_theta)/double(resolution_theta);
			double dphi = (max_phi-min_phi)/double(resolution_phi);

			//std::vector<Color3f> coeffs(numCoeffs, Color3f(0.0f, 0.0f, 0.0f));
			//Color3f* coeffs_ptr = coeffs.data();
			Color3f* coeffs_ptr = shcache.getCoefficients(voxel);

			/*
			// cpu version for computing moments ---
			double pixel_area = dtheta*dphi;

			for( int t=0;t<resolution_theta;++t )
				for( int p=0;p<resolution_phi;++p )
				{
					double phi = dphi*p;
					double theta = dtheta*t;

					Color3f f = map.eval(theta, phi);

					for( int l=0;l<=order;++l )
						for( int m=-l;m<=l;++m )
						{
							double sh = moexp::Y_real(l, m, theta, phi);
							Color3f sample = f*sh;
							coeffs[moexp::shIndex(l, m)]+=sample*pixel_area*std::sin(theta);
						}
				}
			*/

			// gpu version for computing sh coefficients ---
			EnvMap input_coords(V2i(resolution_theta, resolution_phi));
			EnvMap input_f(V2i(resolution_theta, resolution_phi));
			for( int t=0;t<resolution_theta;++t )
				for( int p=0;p<resolution_phi;++p )
				{
					double phi = dphi*p;
					double theta = dtheta*t;

					Color3f f = map.eval(theta, phi);
					input_coords.bitmap().coeffRef( p, t ) = Color3f(theta, phi, 0.0f);
					input_f.bitmap().coeffRef( p, t ) = f;
				}


			CudaPixelBuffer gpu_input_coords;
			gpu_input_coords.initialize(resolution_theta, resolution_phi, (unsigned char*)input_coords.bitmap().data());
			CudaPixelBuffer gpu_input_f;
			gpu_input_f.initialize(resolution_theta, resolution_phi, (unsigned char*)input_f.bitmap().data());

			CudaPtr<CudaPixelBuffer> gpu_input_coords_ptr(&gpu_input_coords);
			CudaPtr<CudaPixelBuffer> gpu_input_f_ptr(&gpu_input_f);
			compute_sh_coefficients( &gpu_input_coords_ptr, &gpu_input_f_ptr, (cumath::V3f*)coeffs_ptr, order, numCoeffs );

			/*
			for( int l=0;l<=order;++l )
			{
				std::string filename = "sh_reconstruction_alt_$0.bgeo";
				filename = replace(filename, "$0", toString(l));

				rasterizeSphericalFunctionSphere(filename, [&](double theta, double phi)->Color3f
				{
					//return moexp::Y_real_sum<Color3f>(order, coeffs.data(), theta, phi);
					return moexp::Y_real_sum<Color3f>(l, coeffs.data(), theta, phi);
				}, 8.0);
			}
			*/

		}
	}
	timer.stop();

	std::cout << "path tracing took " << timer.elapsedSeconds()/double(numVoxels) << "s/voxel in average" << std::endl;

	//shcache.save("nebulae200.shcache");
}
