#include <volume.h>


namespace volumes
{



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

		SGGX3D( double S_xx, double S_xy, double S_xz, double S_yy, double S_yz, double S_zz ):
			S_xx(S_xx),
			S_xy(S_xy),
			S_xz(S_xz),
			S_yy(S_yy),
			S_yz(S_yz),
			S_zz(S_zz)
		{
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

		double sample_spec( const V3d& wi, V3d& wo, double& pdf, RNGd& rng ) const
		{
			// sample VNDF
			const V3d wm = sample_VNDF(wi, rng.next1D(), rng.next1D());
			// specular reflection
			wo = -wi + 2.0f * wm * dot(wm, wi);
			pdf = eval_spec(wi, wo);
			return 1.0; // eval==pdf
		}

		double eval_spec( const V3d& wi, const V3d& wo )const
		{
			V3d wh = wi+wo;

			// normalize wh and detect zero length
			double length = wh.norm();

			if( length != 0.0 )
				wh = V3d(wh.x()/length, wh.y()/length, wh.z()/length);
			else
				return 0.0;
			return 0.25*D(wh)/sigma(wi);
		}

		SGGX3D operator*( const double & scalar )const
		{
			return SGGX3D(S_xx*scalar, S_xy*scalar, S_xz*scalar, S_yy*scalar, S_yz*scalar, S_zz*scalar );
		}

		SGGX3D operator+( const SGGX3D& other)const
		{
			return SGGX3D(S_xx+other.S_xx,
						  S_xy+other.S_xy,
						  S_xz+other.S_xz,
						  S_yy+other.S_yy,
						  S_yz+other.S_yz,
						  S_zz+other.S_zz );
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



	struct SGGXGrid : public PhaseFunctionField
	{
		typedef std::shared_ptr<SGGXGrid> Ptr;

		// store internally using 6 bytes (as explained in the papeer)
		struct PackedSGGXParms
		{
			PackedSGGXParms()
			{
				setSGGX(SGGX3D( Framed(V3d(0.0, 0.0, 1.0)), V3d(1.0, 1.0, 1.0) ));
			}


			void setSGGX( const SGGX3D& sggx )
			{
				// see paper Eq.19
				quantized_sigma_x = clamp<unsigned char>(unsigned char(std::sqrt(sggx.S_xx)*255.0), 0, 255);
				quantized_sigma_y = clamp<unsigned char>(unsigned char(std::sqrt(sggx.S_yy)*255.0), 0, 255);
				quantized_sigma_z = clamp<unsigned char>(unsigned char(std::sqrt(sggx.S_zz)*255.0), 0, 255);

				quantized_r_xy = clamp<unsigned char>(unsigned char(((sggx.S_xy/std::sqrt(sggx.S_xx*sggx.S_yy)+1.0)*0.5)*255.0), 0, 255);
				quantized_r_xz = clamp<unsigned char>(unsigned char(((sggx.S_xz/std::sqrt(sggx.S_xx*sggx.S_zz)+1.0)*0.5)*255.0), 0, 255);
				quantized_r_yz = clamp<unsigned char>(unsigned char(((sggx.S_yz/std::sqrt(sggx.S_yy*sggx.S_zz)+1.0)*0.5)*255.0), 0, 255);
			}

			SGGX3D getSGGX()const
			{
				// see paper Eq.20
				double sigma_x = quantized_sigma_x/255.0;
				double sigma_y = quantized_sigma_y/255.0;
				double sigma_z = quantized_sigma_z/255.0;
				double r_xy = quantized_r_xy/255.0*2.0-1.0;
				double r_xz = quantized_r_xz/255.0*2.0-1.0;
				double r_yz = quantized_r_yz/255.0*2.0-1.0;

				return SGGX3D(sigma_x*sigma_x,      // S_xx
							  r_xy*sigma_x*sigma_y, // S_xy
							  r_xz*sigma_x*sigma_z, // S_xz
							  sigma_y*sigma_y,      // S_yy
							  r_yz*sigma_y*sigma_z, // S_yz
							  sigma_z*sigma_z       // S_zz
							  );
			}

			unsigned char quantized_sigma_x;
			unsigned char quantized_sigma_y;
			unsigned char quantized_sigma_z;
			unsigned char quantized_r_xy;
			unsigned char quantized_r_xz;
			unsigned char quantized_r_yz;
		};

		SGGXGrid():
			PhaseFunctionField(),
			m_localToWorld()
		{
			V3i resolution = V3i(32, 32, 32);
			resize(resolution);

			SGGX3D sggx_iso(Framed(V3d(0.0, 0.0, 1.0)), V3d(1.0, 1.0, 1.0));
			SGGX3D sggx_dir(Framed(V3d(0.0, 0.0, 1.0)), V3d(1.0, 1.0, 0.1));

			RNGd rng(123);
			for( int k=0;k<resolution.z();++k )
				for( int j=0;j<resolution.y();++j )
					for( int i=0;i<resolution.x();++i )
					{
						P3d pVS(i+0.5, j+0.5, k+0.5);
						P3d pLS = voxelToLocal(pVS);

						// distance from center axis
						Ray3d center_axisLS( P3d(0.5, 1.0, 0.1), V3d(0.0, -1.0, 0.0) );
						double distance = (closestPoint(center_axisLS, pLS)-pLS).norm();
						if(distance < 0.1)
						{
							//Framed frame = Framed(V3d(0.0, 0.0, 1.0));
							//V3d projectedAreas=V3d(0.1, 0.1, 1.0);
							//Framed frame = Framed(sampleSphere(rng));
							//V3d projectedAreas=V3d(rng.next1D(), rng.next1D(), rng.next1D());
							//SGGX3D sggx = sggx_dir;
							SGGX3D sggx = lerp(sggx_iso, sggx_dir, 1.0);

							lvalue(i, j, k).setSGGX(sggx);
						}
					}




			// temp ----
			//Framed frame = Framed(V3d(0.0, 0.0, 1.0));
			//V3d projectedAreas=V3d(1.0, 1.0, 1.0);
		}

		// pWS is the world position at which the phase function is to be evaluated
		// wi is the incident light direction (pointing from the scattering point towards the light)
		// wo is the outgoing light direction (pointing from the scattering point towards the camera)
		virtual double eval( const P3d& pWS, const V3d& wi, const V3d& wo )const override
		{
			SGGX3D sggx = evaluate(worldToVoxel(pWS));
			return sggx.eval_spec(wi, wo);
		}


		// pWS is the world position at which the phase function is to be sampled
		// wi point outwards, away from the scattering event
		// wo points outwards, away from the scattering event
		// this is inline with the convention for bsdfs and mitsuba
		// pdf is given in solid angle measure
		virtual double sample( const P3d& pWS, const V3d& wi, V3d& wo, double& pdf, RNGd& rng )const override
		{
			SGGX3D sggx = evaluate(worldToVoxel(pWS));
			return sggx.sample_spec(wi, wo, pdf, rng);
		}

		virtual void setLocalToWorld( const Transformd& localToWorld ) override
		{
			m_localToWorld = localToWorld;
		}

	private:

		const PackedSGGXParms& sample( int i, int j, int k )const
		{
			return m_data[k*m_resolution.x()*m_resolution.y() + j*m_resolution.x() + i];
		}

		PackedSGGXParms& lvalue( int i, int j, int k )
		{
			return m_data[k*m_resolution.x()*m_resolution.y() + j*m_resolution.x() + i];
		}

		SGGX3D evaluate( const P3d& pVS )const
		{
			// take sample location within voxel into account
			V3d vs = pVS - V3d(0.5, 0.5, 0.5);
			double tx = vs.x() - floor(vs.x());
			double ty = vs.y() - floor(vs.y());
			double tz = vs.z() - floor(vs.z());

			// lower left corner
			V3i c1;
			c1[0] = (int)floor(vs.x());
			c1[1] = (int)floor(vs.y());
			c1[2] = (int)floor(vs.z());

			// upper right corner
			V3i c2 = c1+V3i(1);
			V3i res = m_resolution;

			// clamp the indexing coordinates
			c1[0] = std::max(0, std::min(c1.x(), res.x()-1));
			c2[0] = std::max(0, std::min(c2.x(), res.x()-1));
			c1[1] = std::max(0, std::min(c1.y(), res.y()-1));
			c2[1] = std::max(0, std::min(c2.y(), res.y()-1));
			c1[2] = std::max(0, std::min(c1.z(), res.z()-1));
			c2[2] = std::max(0, std::min(c2.z(), res.z()-1));

			//return sample( c1.x(), c1.y(), c1.z() ).getSGGX();
			return lerp( lerp( lerp( sample( c1.x(), c1.y(), c1.z() ).getSGGX(),
									 sample( c2.x(), c1.y(), c1.z() ).getSGGX(), (double)tx ),
							   lerp( sample( c1.x(), c2.y(), c1.z() ).getSGGX(),
									 sample( c2.x(), c2.y(), c1.z() ).getSGGX(), (double)tx ), (double)ty ),
						 lerp( lerp( sample( c1.x(), c1.y(), c2.z() ).getSGGX(),
									 sample( c2.x(), c1.y(), c2.z() ).getSGGX(), (double)tx ),
							   lerp( sample( c1.x(), c2.y(), c2.z() ).getSGGX(),
									 sample( c2.x(), c2.y(), c2.z() ).getSGGX(), (double)tx ), (double)ty ), (double)tz );
		}

		void resize( const V3i& resolution )
		{
			m_resolution = resolution;
			m_data.resize(m_resolution.x()*m_resolution.y()*m_resolution.z());
		}

		P3d localToVoxel( const P3d& pLS )const
		{
			return P3d( pLS.x()*m_resolution.x(),
						pLS.y()*m_resolution.y(),
						pLS.z()*m_resolution.z());
		}

		P3d voxelToLocal( const P3d& pVS )const
		{
			return P3d( pVS.x()/m_resolution.x(),
						pVS.y()/m_resolution.y(),
						pVS.z()/m_resolution.z());
		}

		P3d worldToLocal( const P3d& pWS )const
		{
			return m_localToWorld.inverse()*pWS;
		}

		P3d worldToVoxel( const P3d& pWS )const
		{
			return localToVoxel( worldToLocal(pWS) );
		}

		V3i                          m_resolution;
		std::vector<PackedSGGXParms> m_data;
		Transformd                   m_localToWorld;
	};

	PhaseFunction::Ptr phase_sggx_specular( Framed frame, V3d projectedAreas )
	{
		SGGX3D sggx(frame, projectedAreas);
		return std::make_shared<SGGXSpecular>(sggx);
	}

	PhaseFunctionField::Ptr phasefield_sggx()
	{
		SGGXGrid::Ptr sggxfield = std::make_shared<SGGXGrid>();


		return sggxfield;
	}


} // namespace volumes
