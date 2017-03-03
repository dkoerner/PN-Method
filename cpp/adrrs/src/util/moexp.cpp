#include <util/moexp.h>


#include <map>
#include <util/cas.h>







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


