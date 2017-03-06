#include <util/moexp.h>


#include <map>
#include <math/transform.h>
#include <util/cas.h>
#include <util/bitmap.h>
#include <houio/Geometry.h>
#include <houio/HouGeoIO.h>






namespace moexp
{

	int ipow(int base, int exp)
	{
		int result = 1;
		while (exp)
		{
			if (exp & 1)
				result *= base;
			exp >>= 1;
			base *= base;
		}

		return result;
	}

	//const double M_PI = 3.14159265358979323846;

	// Number of precomputed factorials and double-factorials that can be
	// returned in constant time.
	const int kCacheSize = 16;

	// Compute the factorial for an integer @x. It is assumed x is at least 0.
	// This implementation precomputes the results for low values of x, in which
	// case this is a constant time lookup.
	//
	// The vast majority of SH evaluations will hit these precomputed values.
	double factorial(int x)
	{
		const double factorial_cache[kCacheSize] = {1, 1, 2, 6, 24, 120, 720, 5040,
												  40320, 362880, 3628800, 39916800,
												  479001600, 6227020800,
												  87178291200, 1307674368000};

		if (x < kCacheSize)
		{
			return factorial_cache[x];
		}else
		{
			double s = 1.0;
			for (int n = 2; n <= x; n++)
			{
				s *= n;
			}
			return s;
		}
	}

	// Compute the double factorial for an integer @x. This assumes x is at least
	// 0.  This implementation precomputes the results for low values of x, in
	// which case this is a constant time lookup.
	//
	// The vast majority of SH evaluations will hit these precomputed values.
	// See http://mathworld.wolfram.com/DoubleFactorial.html
	double doubleFactorial(int x)
	{
	  const double dbl_factorial_cache[kCacheSize] = {1, 1, 2, 3, 8, 15, 48, 105,
													  384, 945, 3840, 10395, 46080,
													  135135, 645120, 2027025};

		if (x < kCacheSize)
		{
			return dbl_factorial_cache[x];
		}else
		{
			double s = 1.0;
			double n = x;
			while (n > 1.0)
			{
				s *= n;
				n -= 2.0;
			}
			return s;
		}
	}

	int numTensorComponents( int order, int dimension )
	{
		return ipow(dimension, order);
	}

	// condon-shortley phase
	double csp( int m )
	{
		return (m % 2 == 0 ? 1.0 : -1.0);
	}

	double C(int l, int m)
	{
		double a = csp(m);
		double b1 = (2*l+1)*INV_FOURPI;
		double b2 = factorial(l-m)/factorial(l+m);
		return a*std::sqrt(b1*b2);
		return a;
	}

	double a( int l, int m, int j )
	{
		double a = csp(j);
		double b = ipow(2, l)*factorial(j)*factorial(l-j);
		double c = factorial(2*l-2*j);
		double frac1 = a/b;
		double frac2;
		double d_fac = l-m-2*j;
		// it appears that fractions which contain negative factorials are considered zero by convention
		// see http://mathoverflow.net/questions/10124/the-factorial-of-1-2-3#comment14723_10129
		if( d_fac < 0)
			frac2 = 0.0;
		else
			frac2 = c/factorial(d_fac);
		return frac1*frac2;
	}


	// Evaluate the associated Legendre polynomial of degree @l and order @m at
	// coordinate @x. The inputs must satisfy:
	// 1. l >= 0
	// 2. 0 <= m <= l
	// 3. -1 <= x <= 1
	// See http://en.wikipedia.org/wiki/Associated_Legendre_polynomials
	//
	// This implementation is based off the approach described in [1],
	// instead of computing Pml(x) directly, Pmm(x) is computed. Pmm can be
	// lifted to Pmm+1 recursively until Pml is found
	//
	// note that the Condon-Shorteley Phase is included...
	double P(int l, int m, double x)
	{
		// Compute Pmm(x) = (-1)^m(2m - 1)!!(1 - x^2)^(m/2), where !! is the double
		// factorial.
		double pmm = 1.0;
		// P00 is defined as 1.0, do don't evaluate Pmm unless we know m > 0
		if (m > 0)
		{
			//double sign = (m % 2 == 0 ? 1 : -1);
			double sign = 1.0;
			pmm = sign * doubleFactorial(2 * m - 1) * pow(1 - x * x, m / 2.0);
		}

		if (l == m)
		{
			// Pml is the same as Pmm so there's no lifting to higher bands needed
			return pmm;
		}

		// Compute Pmm+1(x) = x(2m + 1)Pmm(x)
		double pmm1 = x * (2 * m + 1) * pmm;
		if (l == m + 1)
		{
			// Pml is the same as Pmm+1 so we are done as well
			return pmm1;
		}

		// Use the last two computed bands to lift up to the next band until l is
		// reached, using the recurrence relationship:
		// Pml(x) = (x(2l - 1)Pml-1 - (l + m - 1)Pml-2) / (l - m)
		for (int n = m + 2; n <= l; n++)
		{
			double pmn = (x * (2 * n - 1) * pmm1 - (n + m - 1) * pmm) / (n - m);
			pmm = pmm1;
			pmm1 = pmn;
		}
		// Pmm1 at the end of the above loop is equal to Pml
		return pmm1;
	}





	moexp::complex Y_gg( int l, int m, double theta, double phi )
	{
		return C(l,m)*P(l, m, std::cos(theta))*complex(std::cos(m*phi), std::sin(m*phi));
	}

	moexp::complex Y_cc( int l, int m, double theta, double phi )
	{
		return C(l,m)*P(l, m, std::cos(theta))*complex(std::cos(m*phi), -std::sin(m*phi));
	}

	complex Y( int l, int m, double theta, double phi )
	{
		if(m>=0)
		{
			return Y_gg(l, m, theta, phi);
		}else
		{
			return csp(m)*Y_cc(l, std::abs(m), theta, phi);
		}
	}

	// real spherical harmonics...
	// l must be >= 0
	// m must be >= -l and <= l
	double Y_real(int l, int m, double theta, double phi)
	{
		//CHECK(l >= 0, "l must be at least 0.");
		//CHECK(-l <= m && m <= l, "m must be between -l and l.");
		double kml = sqrt((2.0 * l + 1) * factorial(l - abs(m)) / (4.0 * M_PI * factorial(l + abs(m))));
		if (m > 0)
		{
			return sqrt(2.0)*kml*cos(m * phi)*csp(m)*P(l, m, cos(theta));
		}else
		if(m < 0)
		{
			return sqrt(2.0)*kml*sin(-m * phi)*csp(-m)*P(l, -m, cos(theta));
		}else
		{
			return kml*P(l, 0, cos(theta));
		}
	}

	std::unique_ptr<std::vector<complex>> project_Y(int order, const SphericalFunction<double>& func, int sample_count)
	{
		std::unique_ptr<std::vector<complex>> coeffs(new std::vector<complex>());
		coeffs->assign(numSHCoefficients(order), 0.0);

		RNGd rng;
		// integrate func times sh basis function to get the coefficient associated with that basis
		for( int i=0;i<sample_count;++i )
		{
			V3d d = sampleSphere<double>(rng);
			P2d theta_phi = sphericalCoordinates(d);
			double theta = theta_phi.x();
			double phi = theta_phi.y();
			double d_pdf = sampleSpherePDF();

			double f = func( theta_phi.x(), theta_phi.y() );

			// integrate and accumulate for each sh basis function
			for( int l=0;l<=order;++l )
				for( int m=-l;m<=l;++m )
				{
					complex sh = Y(l, m, theta, phi);
					complex sample = f*sh/d_pdf;
					complex& sh_coeff = (*coeffs)[shIndex(l, m)];
					sh_coeff += (sample - sh_coeff)/double(i+1);
				}
		}

		return coeffs;
	}

	complex Y_sum(int order, const std::vector<complex>& coeffs, double theta, double phi)
	{
		//CHECK(GetCoefficientCount(order) == coeffs.size(), "Incorrect number of coefficients provided.");
		complex sum = complex(0.0);//Zero<T>();
		for (int l = 0; l <= order; l++)
		{
			for (int m = -l; m <= l; m++)
			{
				sum += Y(l, m, theta, phi) * coeffs[shIndex(l, m)];
			}
		}
		return sum;
	}






	struct YTensorBuilder
	{

		std::vector<Tensor<complex>> compute_Y_tensors( int l, int m )
		{
			buildComponentCounts(l);

			// this will hold all the tensors by their rank
			std::map<int, Tensor<complex>> tensors;

			cas::Expression::Ptr expanded;
			cas::Scope scope;
			scope.m_variables["l"] = cas::num(l);

			// here we instantiate the variables l and m and expand the equation 2.7 in the Thorne paper
			if( m>= 0 )
			{
				scope.m_variables["m"] = cas::num(m);
				expanded = expand(expr_Y->deep_copy(), scope);
				if(l==3 && m==0)
					std::cout << expanded->toLatex() << std::endl;
			}else
			{
				scope.m_variables["m"] = cas::num(std::abs(m));
				expanded = expand(expr_Y_negative_m->deep_copy(), scope);
			}

			// now we extract the coefficients by analysing the resulting terms of the expanded equation
			cas::Addition::Ptr add = std::dynamic_pointer_cast<cas::Addition>(expanded);
			if(!add)
			{
				std::cout << "YTensorBuilder::compute_Y_tensors: error: equation expansion failed\n";
				throw std::runtime_error("hmghmhgm");
			}

			int numTerms = add->getNumOperands();
			for( int i=0;i<numTerms;++i )
			{
				cas::Multiplication::Ptr mul = std::dynamic_pointer_cast<cas::Multiplication>(add->getOperand(i));
				if(!mul)
				{
					std::cout << "YTensorBuilder::compute_Y_tensors: error: equation expansion failed\n";
					throw std::runtime_error("hmghmhgm");
				}
				//std::cout << "term " << i << "= " << mul->toLatex() << std::endl;

				// each term contributes to one tensor component
				double component_contribution = 1.0;

				// in each term, we expect l number of variables to occur, from
				// which we can figure out the component...this is called the component code
				V3i code(0,0,0);

				bool flag = false;


				int numFactors = mul->getNumOperands();
				for( int j=0;j<numFactors;++j )
				{
					cas::Expression::Ptr factor = mul->getOperand(j);

					cas::Variable::Ptr var = std::dynamic_pointer_cast<cas::Variable>(factor);
					cas::Power::Ptr pow = std::dynamic_pointer_cast<cas::Power>(factor);

					//std::cout << "\tfactor " << j << "= " << factor->toLatex() << std::endl;

					if( std::dynamic_pointer_cast<cas::Number>(factor) )
					{
						cas::Number::Ptr n = std::dynamic_pointer_cast<cas::Number>(factor);
						double value = 1.0;
						switch(n->getType())
						{
							case cas::Number::EInteger:value = n->get_int();break;
							case cas::Number::EReal:value = n->get_real();break;
							default:
							{
								std::cout << "YTensorBuilder::compute_Y_tensors: error: unable to handle number type\n";
								throw std::runtime_error("hrhhrhth");
							}break;
						};

						component_contribution *= value;
						//std::cout << "\tgot number  value=" << value << std::endl;
					}else
					// check for Clm or almj
					if( std::dynamic_pointer_cast<cas::Index>(factor) )
					{
						double index_contribution = 1.0;
						cas::Index::Ptr index = std::dynamic_pointer_cast<cas::Index>(factor);
						if(index->getBaseName() == "C")
							index_contribution *= C(index->getExponent(0), index->getExponent(1));
						else
						if(index->getBaseName() == "a")
							index_contribution *= a(index->getExponent(0), index->getExponent(1), index->getExponent(2));

						component_contribution *= index_contribution;
						//std::cout << "\tgot index " << index->getBaseName() << " value=" << index_contribution << std::endl;
					}else
					if(var)
					{
						if(var->getName() == "n_x")
							code[0] += 1;
						else
						if(var->getName() == "in_y")
							code[1] += 1;
						else
						if(var->getName() == "n_z")
							code[2] += 1;
					}else
					if(!var && pow)
					{
						var = std::dynamic_pointer_cast<cas::Variable>(pow->getOperand(0));
						cas::Number::Ptr exp = std::dynamic_pointer_cast<cas::Number>(pow->getOperand(1));
						if(!var)
						{
							std::cout << "YTensorBuilder::compute_Y_tensors: error: power of non variable encountered\n";
							throw std::runtime_error("dsfhdfhdfh");
						}
						if(!exp)
						{
							std::cout << "YTensorBuilder::compute_Y_tensors: error: power with non number exponent encountered\n";
							throw std::runtime_error("dsfhdfhdfh");
						}
						if(exp->getType() != cas::Number::EInteger)
						{
							std::cout << "YTensorBuilder::compute_Y_tensors: error: power with non integer exponent encountered\n";
							throw std::runtime_error("dsfhdfhdfh");
						}

						int exp_number = exp->get_int();
						if(exp_number >=0)
						{
							if(var->getName() == "n_x")
								code[0] += exp_number;
							else
							if(var->getName() == "in_y")
								code[1] += exp_number;
							else
							if(var->getName() == "n_z")
								code[2] += exp_number;
						}else
							flag = true;
					}else
					{
						std::cout << "YTensorBuilder::compute_Y_tensors: error: unable to handle factor type\n";std::flush(std::cout);
						std::cout << factor->toLatex() << std::endl;
						throw std::runtime_error("hrhhrhth");
					}

				} // for each factor of current term


				if(flag)
				{
					//std::cout << "YTensorBuilder::compute_Y_tensors: warning: variable with negative exponent found contribution=" << component_contribution << std::endl;
					continue;
				}

				// complete the final contribution associated with the current term
				complex final_contribution = std::pow( complex(0.0, 1.0), code[1] )*component_contribution;

				// now we need to find out which tensor component this coefficient is associated with
				// find the rank of the tensor to which the current component belongs
				int component_l = code[0] + code[1] + code[2];

				// look for the tensor in our map and create one if needed
				if( tensors.find(component_l) == tensors.end() )
					tensors[component_l] = Tensor<complex>(component_l);


				// we use the tensor associated with the rank which is implied by the current component code
				Tensor<complex>* tensor = &tensors[component_l];
				for( auto it = tensor->begin(), end = tensor->end();it!=end;++it )
				{
					if( it.code() == code )
					{
						it.value() += final_contribution/double(m_component_count[code]);
					}
				}


			} // for each term in the addition


			std::vector<Tensor<complex>> tensor_list;
			for( auto& it:tensors )
				tensor_list.push_back(it.second);
			return tensor_list;
		}

		static YTensorBuilder* instance()
		{
			if(!g_instance)
				g_instance = new YTensorBuilder();
			return g_instance;
		}


	private:
		YTensorBuilder():
			m_built_level(-1)
		{
			{
				using namespace cas;
				expr_Y = mul( index(var("C"), var("l"), var("m")), pow(add(var("n_x"), var("in_y")), var("m")),sum( "j", largest_integer(mul(pow(num(2), num(-1)) , add(var("l"), mul(num(-1), var("m"))))), mul( index( var("a"), var("l"), var("m"), var("j")), pow( var("n_z"), add(var("l"), mul( num(-1), var("m")), mul( num(-2), var("j")))))));
				expr_Y_negative_m = mul( pow(num(-1), var("m")), index(var("C"), var("l"), var("m")), pow(add(var("n_x"), mul(num(-1), var("in_y"))), var("m")),sum( "j", largest_integer(mul(pow(num(2), num(-1)) , add(var("l"), mul(num(-1), var("m"))))), mul( index( var("a"), var("l"), var("m"), var("j")), pow( var("n_z"), add(var("l"), mul( num(-1), var("m")), mul( num(-2), var("j")))))));
			}
		}

		struct compare_V3i
		{
			bool operator()(const V3i& a, const V3i& b) const
			{
				return std::make_tuple(a.x(), a.y(), a.z()) < std::make_tuple(b.x(), b.y(), b.z());
			}
		};

		void buildComponentCounts( int order )
		{
			while( m_built_level < order )
			{
				++m_built_level;
				int l = m_built_level;
				{
					Tensor<double> t(l);
					for( auto it = t.begin(), end = t.end(); it!=end;++it )
					{
						V3i code = it.code();
						if( m_component_count.find(code) == m_component_count.end() )
							m_component_count[code] = 0;
						++m_component_count[code];
					}
				}
			}
		}


		// this is used for extracting the tensor components from an expanded expression for Y
		cas::Expression::Ptr expr_Y;
		cas::Expression::Ptr expr_Y_negative_m;

		// this contains for each component the number of times it appears in its tensor
		// this is used for our way of dealing with symmetries for tensors of arbitrary rank
		std::map<V3i, int, compare_V3i> m_component_count;
		int m_built_level; // just a way to remember which component counts we already have built

		static YTensorBuilder* g_instance;
	};

	YTensorBuilder* YTensorBuilder::g_instance = 0;


	std::vector<Tensor<complex>> Y_tensors( int l, int m )
	{
		YTensorBuilder* builder = YTensorBuilder::instance();
		return builder->compute_Y_tensors( l, m );
	}
} // namespace moexp




struct EnvMap
{
	typedef std::shared_ptr<EnvMap> Ptr;

	EnvMap( const std::string& filename )
	{
		m_bitmap = Bitmap(filename);
		m_transform = Transformd();
	}

	EnvMap()
	{
		m_bitmap = Bitmap(V2i(512, 256));
		m_transform = Transformd();
	}

	// evaluate environment map
	Color3f eval( double theta, double phi )const
	{
		V3d d = sphericalDirection<double>(theta, phi);
		P2d uv = directionToUV(d);
		return m_bitmap.eval(uv);
	}


	P2d directionToUV( const V3d& d )const
	{
		// using formulas given in http://gl.ict.usc.edu/Data/HighResProbes/
		// with the difference that u=[0,1] (instead of [0,2]) and we negate z
		P2d uv( (1+std::atan2(d.x(), d.z())/M_PI)/2,
				 safe_acos(d.y())/M_PI );
		return uv;
	}
	V3d uvToDirection( const P2d& uv )const
	{
		// using formulas given in http://gl.ict.usc.edu/Data/HighResProbes/
		// with the difference that u=[0,1] (instead of [0,2]) and we negate z
		// azimuthal angle
		double theta = M_PI*(uv.x()*2.0-1.0);
		// elevation angle
		double phi = M_PI*uv.y();
		return V3d( std::sin(phi)*std::sin(theta), std::cos(phi), std::sin(phi)*cos(theta) );
	}
	P2d uvToXY(const P2d& uv)const
	{
		P2d xy(
			(uv.x()*(m_bitmap.cols()-1)),
			(uv.y()*(m_bitmap.rows()-1))
			);
		return xy;
	}

	P2d xyToUV(const P2d& xy)const
	{
		return P2d(
			(xy.x())/double(m_bitmap.cols()-1),
			(xy.y())/double(m_bitmap.rows()-1)
			);
	}

	V3d xyToDirection( const P2d& xy )const
	{
		return uvToDirection( xyToUV(xy) );
	}
	P2d directionToXY( const V3d& d )const
	{
		return uvToXY(directionToUV(d));
	}

	Bitmap& bitmap()
	{
		return m_bitmap;
	}

private:
	Transformd m_transform;
	Bitmap m_bitmap;
};


void rasterizeSphericalFunctionSphere(const std::string& filename, moexp::SphericalFunction<Color3f> func )
{
	houio::Geometry::Ptr geo = houio::Geometry::createSphere(120, 120, 1.0);
	houio::Attribute::Ptr pAttr = geo->getAttr("P");
	houio::Attribute::Ptr cdAttr = houio::Attribute::createV3f(pAttr->numElements());
	for( int i=0;i<pAttr->numElements();++i )
	{
		houio::math::V3f p = pAttr->get<houio::math::V3f>(i);
		P2d theta_phi = sphericalCoordinates<double>(V3d(p.x, p.y, p.z));
		double theta = theta_phi.x();
		double phi = theta_phi.y();
		Color3f col = func(theta, phi);
		cdAttr->set<houio::math::V3f>( i, houio::math::V3f(col.r(), col.g(), col.b()) );
	}
	geo->setAttr("Cd", cdAttr);
	houio::HouGeoIO::xport( filename, geo);
}

void rasterizeSphericalFunctionMap( EnvMap& envmap, moexp::SphericalFunction<Color3f> func )
{
	//std::ofstream f( "test_pixels_coords.txt", std::ios::binary | std::ios::trunc );
	int xres = envmap.bitmap().cols();
	int yres = envmap.bitmap().rows();
	for( int j=0;j<yres;++j )
		for( int i=0;i<xres;++i )
		{
			P2d xy( i+0.5f, j+0.5f );
			V3d d = envmap.xyToDirection(xy);
			P2d theta_phi = sphericalCoordinates(d);
			double theta = theta_phi.x();
			double phi = theta_phi.y();
			envmap.bitmap().coeffRef(j, i) = func(theta, phi);
			//f << theta << " "  << phi << std::endl;
		}
}



// compute sh coefficients
void validate_moexp( int order )
{
	EnvMap map("envmap.exr");


	moexp::SphericalFunction<double> fun = [&](double theta, double phi) -> double
	{
		return map.eval(theta, phi).getLuminance();
	};

	rasterizeSphericalFunctionSphere( "groundtruth.bgeo", fun );

	// complex sh ---
	{
		// project
		int order = 30;
		int numSamples = 50000;
		std::unique_ptr<std::vector<moexp::complex>> sh_coeffs;
		sh_coeffs = moexp::project_Y(order, fun, numSamples);

		// reconstruct
		for( int l=0;l<=order;++l )
		{
			moexp::SphericalFunction<double> reconstruction = [&](double theta, double phi) -> double
			{
				// for some reason our complexSH reconstruction is flipped in y
				// when compared to real SH reconstruction
				V3d d = sphericalDirection(theta, phi);
				d.y() = -d.y();
				double theta2, phi2;
				P2d theta_phi2 = sphericalCoordinates(d);
				theta2 = theta_phi2.x();
				phi2 = theta_phi2.y();

				moexp::complex value = moexp::Y_sum(l, *sh_coeffs.get(), theta2, phi2);
				//return std::abs(value);
				return value.real();
			};

			std::string filename("testfit_complexsh_reconstruction_$0.bgeo");
			filename = replace(filename, "$0", toString(l));

			rasterizeSphericalFunctionSphere( filename, reconstruction );
		}
	}

	// real sh
	{
		// project
		int order = 30;
		int numSamples = 50000;
		std::unique_ptr<std::vector<double>> sh_coeffs;
		sh_coeffs = moexp::project_Y_real(order, fun, numSamples);

		// reconstruct
		for( int l=0;l<=order;++l )
		{
			moexp::SphericalFunction<double> reconstruction = [&](double theta, double phi) -> double
			{
				return moexp::Y_real_sum<double>(l, *sh_coeffs.get(), theta, phi);
			};
			std::string filename("testfit_realsh_reconstruction_$0.bgeo");
			filename = replace(filename, "$0", toString(l));
			rasterizeSphericalFunctionSphere( filename, reconstruction );
		}
	}

	// moment expansion ---
	{
		// produce tensors...
		std::vector<std::vector<moexp::Tensor<moexp::complex>>> Y_tensors_all;
		for( int l=0;l<order;++l )
		{
			std::cout << "building Y_tensor for order l=" << l << std::endl;
			for( int m=-l;m<=l;++m )
			{
				std::cout << "\tm=" << m << std::endl;
				std::vector<moexp::Tensor<moexp::complex>> tensors = moexp::Y_tensors(l, m);
				Y_tensors_all.push_back(tensors);
			}
		}

		// validate that comples SH basis functions are identical to our ylm tensors ---
		/*
		std::vector<Bitmap::Ptr> error_bitmaps;
		int Ylm_index = 0;
		for( int l=0;l<order;++l )
		{
			for( int m=-l;m<=l;++m, ++Ylm_index )
			{
				//if( (l != 1) && (m!=1))
				//	continue;
				std::vector<moexp::Tensor<moexp::complex>>& Y_tensors = Y_tensors_all[Ylm_index];

				EnvMap map_error;


				moexp::SphericalFunction<Color3f> error = [&](double theta, double phi) -> Color3f
				{
					V3d n = sphericalDirection(theta, phi);
					moexp::complex csh = moexp::Y(l, m, theta, phi);
					moexp::complex csh_from_tensor_contraction;

					for( auto& tensor:Y_tensors )
						csh_from_tensor_contraction += moexp::contract(tensor, n);

					double error = std::abs(csh-csh_from_tensor_contraction);

					return Color3f(error);
				};

				rasterizeSphericalFunction(map_error, error);

				{
					std::string filename("error_$0_$1.exr");
					filename = replace(filename, "$0", toString(l));
					filename = replace(filename, "$1", toString(m));
					map_error.bitmap().saveEXR(filename);
					error_bitmaps.push_back(std::make_shared<Bitmap>(filename));
				}
			}
		}
		compositeImages(error_bitmaps, "error_all.exr");
		*/


		// project
		int numSamples = 50000;
		std::unique_ptr<std::vector<moexp::complex>> sh_coeffs;
		sh_coeffs = moexp::project_Y(order, fun, numSamples);

		// we wedge this for every order
		for( int k=1;k<order;++k )
		{
			// compute F_kl tensors according to eq. 2.13a in [thorn80]
			std::vector<moexp::Tensor<moexp::complex>> F_list;

			for( int l=0;l<k;++l )
				F_list.push_back(moexp::Tensor<moexp::complex>(l));

			for( int l=0;l<k;++l )
				for( int m=-l;m<=l;++m )
				{
					int shindex = moexp::shIndex(l,m);
					std::vector<moexp::Tensor<moexp::complex>> Y_tensors = Y_tensors_all[shindex];
					for( auto& tensor : Y_tensors )
						F_list[tensor.getOrder()].multiplyAdd((*sh_coeffs)[shindex], tensor);
				}

			// reconstruct
			moexp::SphericalFunction<double> reconstruction = [&](double theta, double phi) -> double
			{
				V3d n = sphericalDirection(theta, phi);

				// for some reason our complexSH reconstruction is flipped in y
				// when compared to real SH reconstruction
				n.y() = -n.y();

				moexp::complex value(0.0, 0.0);
				for( auto& F:F_list )
					value += moexp::contract(F, n);

				return value.real();
			};

			std::string filename("testfit_moexp_reconstruction_$0.bgeo");
			filename = replace(filename, "$0", toString(k-1));
			rasterizeSphericalFunctionSphere( filename, reconstruction );
		}
	}
}


