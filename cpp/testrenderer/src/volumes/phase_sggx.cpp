#include <volume.h>


namespace volumes
{

	struct SGGXGrid
	{
		// store internally using 6 bytes (as explained in the papeer)
		// sample(i,j,k) -> gives 6 S_xx, parameters...
	};

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

			M33d m_S = m*s*m.transpose();
			S_xx = m_S.coeffRef(0,0);
			S_xy = m_S.coeffRef(0,1);
			S_xz = m_S.coeffRef(0,2);
			S_yy = m_S.coeffRef(1,1);
			S_yz = m_S.coeffRef(1,2);
			S_zz = m_S.coeffRef(2,2);


			//std::cout << "S=" << m_S << std::endl;


			/*
			// build transform ---
			Eigen::Affine3d unitSphereToEllipsoid;

			double sx = std::pow(S_22*S_33/S_11, 1.0/4)/std::sqrt(M_PI);
			double sy = std::pow(S_11*S_33/S_22, 1.0/4)/std::sqrt(M_PI);
			double sz = std::pow(S_11*S_22/S_33, 1.0/4)/std::sqrt(M_PI);
			unitSphereToEllipsoid = Eigen::Affine3d(m)*Eigen::Scaling(V3d(sx, sy, sz));
			m_unitSphereToEllipsoid = Transformd(unitSphereToEllipsoid);
			*/
		}

		double D( const V3d& wm )const
		{
			const double detS = S_xx*S_yy*S_zz - S_xx*S_yz*S_yz - S_yy*S_xz*S_xz - S_zz*S_xy*S_xy + 2.0*S_xy*S_xz*S_yz;
			const double den = wm.x()*wm.x()*(S_yy*S_zz - S_yz*S_yz) + wm.y()*wm.y()*(S_xx*S_zz - S_xz*S_xz) + wm.z()*wm.z()*(S_xx*S_yy - S_xy*S_xy)
				+ 2.0*(wm.x()*wm.y()*(S_xz*S_yz - S_zz*S_xy) + wm.x()*wm.z()*(S_xy*S_yz - S_yy*S_xz) + wm.y()*wm.z()*(S_xy*S_xz - S_xx*S_yz));
			const double D = std::pow(std::abs(detS), 1.50) / (MATH_PI*den*den);
			return D;
		}

		double sigma(const V3d& d)const
		{
			const float sigma_squared = d.x()*d.x()*S_xx + d.y()*d.y()*S_yy + d.z()*d.z()*S_zz + 2.0 * (d.x()*d.y()*S_xy + d.x()*d.z()*S_xz + d.y()*d.z()*S_yz);
			return (sigma_squared > 0.0) ? std::sqrt(sigma_squared) : 0.0; // conditional to avoid numerical errors
		}

		// build orthonormal basis (Building an Orthonormal Basis from a 3D Unit Vector Without Normalization, [Frisvad2012])
		void buildOrthonormalBasis(V3d& omega_1, V3d& omega_2, const V3d& omega_3) const
		{
			if(omega_3.z() < -0.9999999f)
			{
				omega_1 = V3d( 0.0 , -1.0 , 0.0 );
				omega_2 = V3d( -1.0 , 0.0 , 0.0 );
			}else
			{
				const double a = 1.0 /(1.0 + omega_3.z() );
				const double b = -omega_3.x()*omega_3 .y()*a ;
				omega_1 = V3d (1.0 - omega_3.x()*omega_3.x()*a , b , -omega_3.x() );
				omega_2 = V3d (b , 1.0 - omega_3.y()*omega_3.y()*a , -omega_3.y() );
			}
		}

		V3d sample_VNDF( const V3d wi, const float U1, const float U2)const
		{
			// generate sample (u, v, w)
			const double r = std::sqrt(U1);
			const double phi = 2.0*M_PI*U2;
			const double u = r*std::cos(phi);
			const double v= r*std::sin(phi);
			const double w = std::sqrt(1.0f - u*u - v*v);
			// build orthonormal basis
			V3d wk, wj;
			buildOrthonormalBasis(wk, wj, wi);
			// project S in this basis
			const double S_kk = wk.x()*wk.x()*S_xx + wk.y()*wk.y()*S_yy + wk.z()*wk.z()*S_zz
			+ 2.0 * (wk.x()*wk.y()*S_xy + wk.x()*wk.z()*S_xz + wk.y()*wk.z()*S_yz);
			const double S_jj = wj.x()*wj.x()*S_xx + wj.y()*wj.y()*S_yy + wj.z()*wj.z()*S_zz
			+ 2.0 * (wj.x()*wj.y()*S_xy + wj.x()*wj.z()*S_xz + wj.y()*wj.z()*S_yz);
			const double S_ii = wi.x()*wi.x()*S_xx + wi.y()*wi.y()*S_yy + wi.z()*wi.z()*S_zz
			+ 2.0 * (wi.x()*wi.y()*S_xy + wi.x()*wi.z()*S_xz + wi.y()*wi.z()*S_yz);
			const double S_kj = wk.x()*wj.x()*S_xx + wk.y()*wj.y()*S_yy + wk.z()*wj.z()*S_zz
			+ (wk.x()*wj.y() + wk.y()*wj.x())*S_xy
			+ (wk.x()*wj.z() + wk.z()*wj.x())*S_xz
			+ (wk.y()*wj.z() + wk.z()*wj.y())*S_yz;
			const double S_ki = wk.x()*wi.x()*S_xx + wk.y()*wi.y()*S_yy + wk.z()*wi.z()*S_zz
			+ (wk.x()*wi.y() + wk.y()*wi.x())*S_xy + (wk.x()*wi.z() + wk.z()*wi.x())*S_xz + (wk.y()*wi.z() + wk.z()*wi.y())*S_yz;
			const double S_ji = wj.x()*wi.x()*S_xx + wj.y()*wi.y()*S_yy + wj.z()*wi.z()*S_zz
			+ (wj.x()*wi.y() + wj.y()*wi.x())*S_xy
			+ (wj.x()*wi.z() + wj.z()*wi.x())*S_xz
			+ (wj.y()*wi.z() + wj.z()*wi.y())*S_yz;
			// compute normal
			double sqrtDetSkji = std::sqrt(std::abs(S_kk*S_jj*S_ii - S_kj*S_kj*S_ii - S_ki*S_ki*S_jj - S_ji*S_ji*S_kk + 2.0*S_kj*S_ki*S_ji));
			double inv_sqrtS_ii = 1.0 / std::sqrt(S_ii);
			double tmp = std::sqrt(S_jj*S_ii-S_ji*S_ji);
			V3d Mk(sqrtDetSkji/tmp, 0.0, 0.0);
			V3d Mj(-inv_sqrtS_ii*(S_ki*S_ji-S_kj*S_ii)/tmp , inv_sqrtS_ii*tmp, 0);
			V3d Mi(inv_sqrtS_ii*S_ki, inv_sqrtS_ii*S_ji, inv_sqrtS_ii*S_ii);
			V3d wm_kji = normalize(u*Mk+v*Mj+w*Mi);
			// rotate back to world basis
			return wm_kji.x() * wk + wm_kji.y() * wj + wm_kji.z() * wi;
		}


		double S_xx;
		double S_xy;
		double S_xz;
		double S_yy;
		double S_yz;
		double S_zz;
		//Transformd m_unitSphereToEllipsoid;
	};

	struct SGGXSpecular : public PhaseFunction
	{
		typedef std::shared_ptr<SGGXSpecular> Ptr;

		SGGXSpecular( SGGX3D sggx ):
			PhaseFunction(),
			m_sggx(sggx)
		{
		}

		virtual double sample( const V3d& wi, V3d& wo, double& pdf, RNGd& rng ) const override
		{
			// sample VNDF
			const V3d wm = m_sggx.sample_VNDF(wi, rng.next1D(), rng.next1D());
			// specular reflection
			wo = -wi + 2.0f * wm * dot(wm, wi);
			pdf = eval(wi, wo);
			return 1.0; // eval==pdf
		}

		virtual double eval( const V3d& wi, const V3d& wo )const override
		{
			V3d wh = wi+wo;

			// normalize wh and detect zero length
			double length = wh.norm();

			if( length != 0.0 )
				wh = V3d(wh.x()/length, wh.y()/length, wh.z()/length);
			else
				return 0.0;
			return 0.25*m_sggx.D(wh)/m_sggx.sigma(wi);
		}

		SGGX3D m_sggx;
	};


	PhaseFunction::Ptr phase_sggx_specular( Framed frame, V3d projectedAreas )
	{
		SGGX3D sggx(frame, projectedAreas);
		return std::make_shared<SGGXSpecular>(sggx);
	}


} // namespace volumes
