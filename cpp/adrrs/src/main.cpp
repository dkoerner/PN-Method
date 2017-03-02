#include <iostream>
#include <complex>

#include <scene.h>
#include <integrator.h>
#include <util/threadpool.h>
#include <util/field.h>
#include <util/wedge.h>
#include <util/sh.h>

#include <houio/Geometry.h>

#include <pncache.h>





void render_volume( RenderTaskInfo& ti )
{
	int width = ti.g.scene->camera->getResolutionX();
	int height = ti.g.scene->camera->getResolutionY();


	// note that y goes from bottom=0 to top=max

	for (int scanline = ti.taskid; scanline < ti.g.crop_window.max.y(); scanline += ti.numTasks)
	{
		if( scanline < ti.g.crop_window.min.y() )
			continue;

		for (int x = ti.g.crop_window.min.x(); x < ti.g.crop_window.max.x(); ++x)
		{
			int y = height-1-scanline;
			int index = y*width + x;
			Color3f f(0.0f);
			Color3f T(1.0f);


			bool debug = false;

			if( ti.g.debug_pixel.x()>0 )
			{
				debug = (x == ti.g.debug_pixel.x()) &&
						(y == ti.g.debug_pixel.y());

				if(!debug)
					continue;
			}

			Ray3d rayWS;
			ti.g.scene->camera->sampleRay( P2d(x+0.5, y+0.5), rayWS, debug );
			//ti.g.scene->camera->sampleRay( P2d(x+ti.rng.next1D(), y+ti.rng.next1D()), rayWS );

			// do raycast ---
			try
			{
				RadianceQuery rq;
				rq.ray = rayWS;
				rq.pixel = V2i(x, y);
				rq.debug = debug;
				f = ti.g.scene->integrator->Li(ti.g.scene, rq, ti.rng);
				T = rq.transmittance;
			}
			catch (std::exception& e)
			{
				std::cout << "render_volume: caught exception at task=" << ti.taskid << " x=" << x << " y=" << y << " " << " index=" <<  index << " sample=" << ti.samples << std::endl;
				std::cout << e.what() << '\n';
				std::flush(std::cout);
				throw e;
			}


			if( std::isnan(f.getLuminance()) )
				std::cout << "PathTracingTask::run: got NaN value @ index=" << index << " sample=" << ti.samples << std::endl;

			// update pixel color
			Color3f& c = ti.g.image->coeffRef(index);
			c += (f - c)/float(ti.samples+1);

			// update transmittance
			Color3f& c_transmittance = ti.g.image_transmittance->coeffRef(index);
			c_transmittance += (T - c_transmittance)/float(ti.samples+1);

			if(debug)
			{
				c = Color3f(1.0f, 0.0f, 0.0f);
				c_transmittance = Color3f(1.0f, 0.0f, 0.0f);
			}
		} // pixel x
	} // pixel y


	++ti.samples;
}




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


void writeSphericalFunction(const std::string& filename, sh::SphericalFunction<Color3f> func )
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

void rasterizeSphericalFunction( EnvMap& envmap, sh::SphericalFunction<Color3f> func )
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

double test( int n, int l )
{
	int n2 = 2*n;
	int l2 = 2*l;
	return tensor::ipow(-1, n)*double(sh::factorial(l)*sh::doubleFactorial(l2-n2-1))/double(sh::factorial(l-n2)*sh::doubleFactorial(l2-1)*sh::doubleFactorial(n2));
}

// condon-shortley phase
double csp( int m )
{
	return (m % 2 == 0 ? 1.0 : -1.0);
}

double Clm( int l, int m )
{
	//double a = csp(m);
	double a = 1.0;
	double b1 = (2*l+1)*INV_FOURPI;
	double b2 = sh::factorial(l-m)/sh::factorial(l+m);
	return a*std::sqrt(b1*b2);
	return a;
}

double Clm_including_csp( int l, int m )
{
	double a = csp(m);
	double b1 = (2*l+1)*INV_FOURPI;
	double b2 = sh::factorial(l-m)/sh::factorial(l+m);
	return a*std::sqrt(b1*b2);
	return a;
}

using complex = std::complex<double>;
double sh_eval2( int l, int m, double theta, double phi )
{
	std::complex<double> test2(0.0, 1.0);
	std::complex<double> test = std::exp(test2*phi*double(m));
	return test.real()*sh::P(l, m, std::cos(theta));
}

complex Y( int l, int m, double theta, double phi )
{
	return Clm(l,m)*sh::P(l, m, std::cos(theta))*complex(std::cos(m*phi), std::sin(m*phi));
}

complex Y_cc( int l, int m, double theta, double phi )
{
	return Clm(l,m)*sh::P(l, m, std::cos(theta))*complex(std::cos(m*phi), -std::sin(m*phi));
}

complex complex_sh( int l, int m, double theta, double phi )
{
	if(m>=0)
		return Y(l, m, theta, phi);
	else
		return csp(m)*Y_cc(l, std::abs(m), theta, phi);
}



double almj( int l, int m, int j, bool debug = false )
{
	double a = csp(j);
	double b = tensor::ipow(2, l)*sh::factorial(j)*sh::factorial(l-j);
	double c = sh::factorial(2*l-2*j);
	double frac1 = a/b;
	double frac2;
	double d_fac = l-m-2*j;
	// it appears that fractions which contain negative factorials are considered zero by convention
	// see http://mathoverflow.net/questions/10124/the-factorial-of-1-2-3#comment14723_10129
	if( d_fac < 0)
		frac2 = 0.0;
	else
		frac2 = c/sh::factorial(d_fac);
	return frac1*frac2;
}

double P2( int l, int m, double theta, bool debug = false )
{
	double a = std::pow( std::sin(theta), double(m) );
	double b = 0.0;
	int j_end = int(std::ceil((l-m)/2.0));
	for( int j=0;j<=j_end;++j )
		b+= almj(l, m, j)*std::pow(std::cos(theta), double(l-m-2*j));
	return csp(m)*a*b;
}


complex Y2( int l, int m, double theta, double phi )
{
	return Clm(l,m)*P2(l, m, theta)*complex(std::cos(m*phi), std::sin(m*phi));
}

complex Y2_cc( int l, int m, double theta, double phi )
{
	return Clm(l,m)*P2(l, m, theta)*complex(std::cos(m*phi), -std::sin(m*phi));
}


complex complex_sh2( int l, int m, double theta, double phi )
{
	if(m>=0)
		return Y2(l, m, theta, phi);
	else
		return csp(m)*Y2_cc(l, std::abs(m), theta, phi);
}


double P3( int l, int m, const V3d& n )
{
	double a = 1.0;
	double b = 0.0;
	int j_end = int(std::ceil((l-m)/2.0));
	for( int j=0;j<=j_end;++j )
		b+= almj(l, m, j)*std::pow(n.z(), double(l-m-2*j));
	return csp(m)*a*b;
}

complex Y3( int l, int m, const V3d& n )
{
	return Clm(l,m)*P3(l, m, n)*std::pow(complex(n.x(), n.y()), double(m));
}

complex Y3_cc( int l, int m, const V3d& n )
{
	return Clm(l,m)*P3(l, m, n)*std::pow(complex(n.x(), -n.y()), double(m));
}


complex complex_sh3( int l, int m, const V3d& n )
{
	if(m>=0)
		return Y3(l, m, n);
	else
		return csp(m)*Y3_cc(l, std::abs(m), n);
}


namespace xp
{


	struct Expression
	{
		typedef std::shared_ptr<Expression> Ptr;
		virtual std::string toLatex()const=0;
		virtual void print_hierarchy( int indent = 0 )const
		{
			for(int i=0;i<indent;++i)std::cout<<"\t";
			std::cout << "Expression" << std::endl;
		}
		virtual Ptr deep_copy()const=0;

	};

	Expression::Ptr num( int value );
	Expression::Ptr var( const std::string& name );
	Expression::Ptr add( Expression::Ptr a, Expression::Ptr b );
	Expression::Ptr add( Expression::Ptr a, Expression::Ptr b, Expression::Ptr c );
	Expression::Ptr sum( const std::string& index, Expression::Ptr end, Expression::Ptr body );
	Expression::Ptr mul( Expression::Ptr a, Expression::Ptr b );
	Expression::Ptr mul( Expression::Ptr a, Expression::Ptr b, Expression::Ptr c );
	Expression::Ptr mul( Expression::Ptr a, Expression::Ptr b, Expression::Ptr c, Expression::Ptr d );
	Expression::Ptr pow( Expression::Ptr base, Expression::Ptr exp );
	Expression::Ptr largest_integer( Expression::Ptr expr );

	Expression::Ptr index( Expression::Ptr base, Expression::Ptr a, Expression::Ptr b, Expression::Ptr c );
	Expression::Ptr index( Expression::Ptr base, Expression::Ptr a, Expression::Ptr b);

	std::string toLatex(Expression::Ptr expr);


	struct Scope
	{
		std::map<std::string, Expression::Ptr> m_variables;
	};

	namespace rewrite
	{
		Expression::Ptr expand( Expression::Ptr expr, Scope& scope );
		Expression::Ptr expand(Expression::Ptr e);
		void replace_variable( Expression::Ptr expr, const std::string& variable, Expression::Ptr replacement );
		Expression::Ptr fold_constants( Expression::Ptr expr );
	}

	struct Variable : public Expression
	{
		typedef std::shared_ptr<Variable> Ptr;
		Variable( const std::string& name )
			:Expression(),
			 m_name(name)
		{
		}

		Variable( const Variable& other )
			:Expression(),
			 m_name(other.m_name)
		{

		}


		virtual Expression::Ptr deep_copy()const override
		{
			return std::make_shared<Variable>(*this);
		}

		const std::string& getName()const
		{
			return m_name;
		}

		virtual std::string toLatex()const override
		{
			return m_name;
		}

		virtual void print_hierarchy( int indent = 0 )const override
		{
			for(int i=0;i<indent;++i)std::cout<<"\t";
			std::cout << "Variable " << m_name << std::endl;
		}

	private:
		std::string m_name;
	};

	struct Number : public Expression
	{
		typedef std::shared_ptr<Number> Ptr;
		Number(int value) :
			Expression()
		{
			set_int(value);
		}

		Number( const Number& other )
			:Expression(),
			 m_type(other.m_type)
		{
			switch(m_type)
			{
			case EInteger:m_i = other.m_i;break;
			case EReal:m_d = other.m_d;break;
			};
		}


		virtual Expression::Ptr deep_copy()const override
		{
			return std::make_shared<Number>(*this);
		}

		enum EType
		{
			EInteger,
			EReal//,
			//EComplex
		};

		virtual std::string toLatex()const override
		{
			switch(m_type)
			{
			case EInteger:return toString(m_i);break;
			case EReal:return toString(m_d);break;
			};

			return "number:n/a";
		}

		virtual void print_hierarchy( int indent = 0 )const override
		{
			for(int i=0;i<indent;++i)std::cout<<"\t";
			std::cout << "Number ";
			switch(m_type)
			{
			case EInteger:std::cout << "integer " << m_i;break;
			case EReal:std::cout << "real " << m_d;break;
			};
			std::cout << std::endl;
		}

		Number& operator*=(const Number& rhs)
		{
			if( (m_type == EInteger) && (rhs.m_type == EInteger) )
				set_int( m_i*rhs.m_i );
			else
			if( (m_type == EInteger) && (rhs.m_type == EReal) )
				set_real( m_i*rhs.m_d );
			else
			if( (m_type == EReal) && (rhs.m_type == EInteger) )
				set_real( m_d*rhs.m_i );
			else
			if( (m_type == EReal) && (rhs.m_type == EReal) )
				set_real( m_d*rhs.m_d );
			else
			{
				std::cout << "Number::*= error: unknown type combination\n";
				throw std::runtime_error("safasfsaf");
			}
			return *this;
		}
		Number& operator+=(const Number& rhs)
		{
			if( (m_type == EInteger) && (rhs.m_type == EInteger) )
				set_int( m_i+rhs.m_i );
			else
			if( (m_type == EInteger) && (rhs.m_type == EReal) )
				set_real( m_i+rhs.m_d );
			else
			if( (m_type == EReal) && (rhs.m_type == EInteger) )
				set_real( m_d+rhs.m_i );
			else
			if( (m_type == EReal) && (rhs.m_type == EReal) )
				set_real( m_d+rhs.m_d );
			else
			{
				std::cout << "Number::*= error: unknown type combination\n";
				throw std::runtime_error("safasfsaf");
			}
			return *this;
		}
		Number& pow(const Number& exp)
		{
			if( (m_type == EInteger) && (exp.m_type == EInteger) )
				set_real( std::pow(m_i, exp.m_i) );
			else
			if( (m_type == EInteger) && (exp.m_type == EReal) )
				set_real( std::pow(double(m_i), exp.m_d) );
			else
			if( (m_type == EReal) && (exp.m_type == EInteger) )
				set_real( std::pow(m_d, double(exp.m_i)) );
			else
			if( (m_type == EReal) && (exp.m_type == EReal) )
				set_real( std::pow(m_d, exp.m_d) );
			else
			{
				std::cout << "Number::*= error: unknown type combination\n";
				throw std::runtime_error("safasfsaf");
			}
			return *this;
		}

		EType getType()
		{
			return m_type;
		}
		int get_int()
		{
			return m_i;
		}
		double get_real()
		{
			return m_d;
		}

	private:
		union
		{
			int m_i;
			double m_d;
		};
		EType m_type;

		void set_int( int value )
		{
			m_i = value;
			m_type = EInteger;
		}
		void set_real( double value )
		{
			m_d = value;
			m_type = EReal;
		}
	};

	struct Operator : public Expression
	{
		typedef std::shared_ptr<Operator> Ptr;
		virtual int getNumOperands()const=0;
		virtual Expression::Ptr getOperand(int index)=0;
		virtual void setOperand(int index, Expression::Ptr expr)=0;
	};

	struct LargestInteger : public Operator
	{
		typedef std::shared_ptr<LargestInteger> Ptr;
		LargestInteger( Expression::Ptr expr ):
			Operator(),
			m_expr(expr)
		{
		}

		virtual Expression::Ptr deep_copy()const override
		{
			return std::make_shared<LargestInteger>(m_expr->deep_copy());
		}

		virtual std::string toLatex()const override
		{
			return "\\lfloor " + m_expr->toLatex() + " \\rfloor";
		}

		virtual int getNumOperands()const
		{
			return 1;
		}
		virtual Expression::Ptr getOperand(int index)
		{
			if( index == 0 )
				return m_expr;
			return Expression::Ptr();
		}

		virtual void setOperand(int index, Expression::Ptr expr)
		{
			if(index==0)
				m_expr = expr;
		}


	private:
		Expression::Ptr m_expr;
	};

	struct Index : public Operator
	{
		typedef std::shared_ptr<Index> Ptr;

		Index(Expression::Ptr base):Operator(),m_base(base)
		{
		}


		virtual Expression::Ptr deep_copy()const override
		{
			Ptr idx = std::make_shared<Index>(m_base->deep_copy());
			for( auto index:m_indices )
				idx->addIndex(index->deep_copy());
			return idx;
		}

		virtual std::string toLatex()const override
		{
			std::string result = m_base->toLatex() + "^{";

			int numIndices = m_indices.size();
			for( int i=0;i<numIndices;++i )
			{
				result+=m_indices[i]->toLatex();
				if( i<numIndices-1 )
					result+=",";
			}
			return result + "}";
		}

		virtual void print_hierarchy( int indent = 0 )const
		{
			for(int i=0;i<indent;++i)std::cout<<"\t";
			std::cout << "Index" << std::endl;
		}

		void addIndex( Expression::Ptr expr )
		{
			m_indices.push_back(expr);
		}

		virtual int getNumOperands()const
		{
			return m_indices.size()+1;
		}
		virtual Expression::Ptr getOperand(int index)
		{
			if( index == 0 )
				return m_base;
			return m_indices[index-1];
		}

		virtual void setOperand(int index, Expression::Ptr expr)
		{
			if(index==0)
				m_base = expr;
			else
			{
				m_indices[index-1] = expr;
			}
		}

		std::string getBaseName()
		{
			Variable::Ptr base = std::dynamic_pointer_cast<Variable>(m_base);
			return base->getName();
		}

		int getExponent(int index)
		{
			Number::Ptr num =  std::dynamic_pointer_cast<Number>(m_indices[index]);
			return num->get_int();
		}


	private:
		Expression::Ptr m_base;
		std::vector<Expression::Ptr> m_indices;
	};

	struct Addition : public Operator
	{
		typedef std::shared_ptr<Addition> Ptr;


		virtual Expression::Ptr deep_copy()const override
		{
			Ptr add = std::make_shared<Addition>();
			for( auto operand:m_operands )
				add->addOperand(operand->deep_copy());
			return add;
		}

		virtual std::string toLatex()const override
		{
			std::string result;
			int numOperands = m_operands.size();
			for( int i=0;i<numOperands;++i )
			{
				result+=m_operands[i]->toLatex();
				if( i<numOperands-1 )
					result+="+";
			}
			return result;
		}

		virtual void print_hierarchy( int indent = 0 )const override
		{
			for(int i=0;i<indent;++i)std::cout<<"\t";
			std::cout << "Addition " << std::endl;
			for( auto it:m_operands )
				it->print_hierarchy(indent+1);
		}

		void addOperand( Expression::Ptr expr )
		{
			m_operands.push_back(expr);
		}

		virtual int getNumOperands()const
		{
			return m_operands.size();
		}
		virtual Expression::Ptr getOperand(int index)
		{
			return m_operands[index];
		}

		virtual void setOperand(int index, Expression::Ptr expr)
		{
			m_operands[index] = expr;
		}

	private:
		std::vector<Expression::Ptr> m_operands;
	};

	struct Multiplication : public Operator
	{
		typedef std::shared_ptr<Multiplication> Ptr;

		virtual Expression::Ptr deep_copy()const override
		{
			Ptr mul = std::make_shared<Multiplication>();
			for( auto operand:m_operands )
				mul->addOperand(operand->deep_copy());
			return mul;
		}

		virtual std::string toLatex()const override
		{
			std::string result;
			int numOperands = m_operands.size();
			for( int i=0;i<numOperands;++i )
			{
				bool parentheses = false;
				Expression::Ptr op = m_operands[i];
				// if operand is an addition operator, we will put it in parentheses
				if( std::dynamic_pointer_cast<Addition>(op) )
					parentheses=true;
				if(parentheses)
					result+= "\\left ( " + op->toLatex() + " \\right )";
				else
					result+= op->toLatex();
				if( i<numOperands-1 )
					result+="";
			}
			return result;
		}

		virtual void print_hierarchy( int indent = 0 )const override
		{
			for(int i=0;i<indent;++i)std::cout<<"\t";
			std::cout << "Multiplication " << std::endl;
			for( auto it:m_operands )
				it->print_hierarchy(indent+1);
		}

		void addOperand( Expression::Ptr expr )
		{
			m_operands.push_back(expr);
		}

		virtual int getNumOperands()const
		{
			return m_operands.size();
		}
		virtual Expression::Ptr getOperand(int index)
		{
			return m_operands[index];
		}

		virtual void setOperand(int index, Expression::Ptr expr)
		{
			m_operands[index] = expr;
		}

	private:
		std::vector<Expression::Ptr> m_operands;
	};

	struct Sum : public Operator
	{
		typedef std::shared_ptr<Sum> Ptr;

		Sum( const std::string& index, Expression::Ptr end, Expression::Ptr body ):
			Operator(),
			m_index(index),
			m_end(end),
			m_body(body)
		{
		}

		virtual Expression::Ptr deep_copy()const override
		{
			return std::make_shared<Sum>(m_index, m_end->deep_copy(), m_body->deep_copy());
		}

		virtual std::string toLatex()const override
		{
			std::string result = "\\sum_"+m_index+"^{" + m_end->toLatex() + "}{" + m_body->toLatex() + "}";
			return result;
		}

		virtual int getNumOperands()const
		{
			return 2;
		}
		virtual Expression::Ptr getOperand(int index)
		{
			if( index == 0 )
				return m_end;
			else
			if( index == 1 )
				return m_body;
			return Expression::Ptr();
		}

		virtual void setOperand(int index, Expression::Ptr expr)
		{
			if(index==0)
				m_end = expr;
			else
			if( index == 1 )
				m_body = expr;
		}

		const std::string& getIndexName()const
		{
			return m_index;
		}


	private:
		std::string m_index;
		Expression::Ptr m_end;
		Expression::Ptr m_body;
	};

	struct Power : public Operator
	{
		typedef std::shared_ptr<Power> Ptr;
		Power(Expression::Ptr base, Expression::Ptr exp):
			Operator(),
			m_base(base),
			m_exp(exp)
		{

		}

		virtual Expression::Ptr deep_copy()const override
		{
			return std::make_shared<Power>( m_base->deep_copy(), m_exp->deep_copy());
		}

		virtual std::string toLatex()const override
		{
			bool doparentheses = false;
			if( std::dynamic_pointer_cast<Operator>(m_base) )
				doparentheses = true;

			std::string result;

			if(doparentheses)
				result = "\\left ( " + m_base->toLatex() + " \\right )";
			else
				result = m_base->toLatex();
			result += "^{" + m_exp->toLatex() + "}";
			return result;
		}

		virtual void print_hierarchy( int indent = 0 )const override
		{
			for(int i=0;i<indent;++i)std::cout<<"\t";
			std::cout << "Power" << std::endl;
			m_base->print_hierarchy(indent+1);
			m_exp->print_hierarchy(indent+1);
		}

		virtual int getNumOperands()const
		{
			return 2;
		}
		virtual Expression::Ptr getOperand(int index)
		{
			if( index == 0 )
				return m_base;
			else
			if( index == 1 )
				return m_exp;
			return Expression::Ptr();
		}

		virtual void setOperand(int index, Expression::Ptr expr)
		{
			if(index==0)
				m_base = expr;
			else
			if( index == 1 )
				m_exp = expr;
		}
	private:
		Expression::Ptr m_base;
		Expression::Ptr m_exp;
	};





}


complex contract_moment( tensor::Tensor<complex>& a, const V3d& d  )
{
	complex result = 0.0;

	// iterate over all components of tensor a
	for( auto it = a.begin(), end = a.end(); it!=end;++it  )
	{
		result += it.weight(d)*it.value();
	}

	return result;
}

struct ylm_t
{
	struct TensorData
	{
		typedef std::shared_ptr<TensorData> Ptr;
		tensor::Tensor<complex> tensor;
		std::vector<complex> tensor_components;
	};

	static std::vector<std::vector<V3i>> g_component_codes;

	struct compare_V3i
	{
		bool operator()(const V3i& a, const V3i& b) const
		{
			return std::make_tuple(a.x(), a.y(), a.z()) < std::make_tuple(b.x(), b.y(), b.z());
		}
	};
	static std::map<V3i, int, compare_V3i> g_component_count;

	static void build(int order)
	{
		g_component_codes.resize(order);
		for( int l=0;l<order;++l)
		{
			std::vector<V3i>& component_codes = g_component_codes[l];

			tensor::Tensor<double> t( 0, l );
			for( auto it = t.begin(), end = t.end(); it!=end;++it )
			{
				V3i code(0,0,0);
				for( int j=0;j<l;++j )
					++code[it.index(j)];
				component_codes.push_back(code);

				if( g_component_count.find(code) == g_component_count.end() )
					g_component_count[code] = 0;
				++g_component_count[code];
			}
		}
	}

	ylm_t( int l, int m ):
		m_l(l),
		m_m(m)
	{
		/*
		m_components_l.resize(tensor::numComponents(l), 0.0);
		m_tensor_l = tensor::Tensor<complex>(m_components_l.data(), m_l);
		if( l>= 2 )
		{
			m_components_lm2.resize(tensor::numComponents(l-2), 0.0);
			m_tensor_lm2 = tensor::Tensor<complex>(m_components_lm2.data(), m_l-2);
		}
		*/
	}

	ylm_t( const ylm_t& other ):
		m_l(other.m_l),
		m_m(other.m_m),
		m_tensordata(other.m_tensordata.begin(), other.m_tensordata.end())
	{
	}


	void print()
	{
		std::cout << "ylm l=" << m_l << " m=" << m_m << std::endl;
		for( auto& it:m_tensordata )
		{
			int l = it.first;
			TensorData::Ptr td = it.second;
			std::cout << "\tl=" << l << std::endl;
			int numComponents = td->tensor_components.size();
			for( int i=0;i<numComponents;++i )
			{
				std::cout << "\t" << g_component_codes[l][i].toString() << "=" << td->tensor_components[i] << std::endl;
			}
		}
	}

	void add_contribution_to_component( const V3i& code, complex contribution )
	{
		// find the rank of the tensor to which the current component belongs
		int component_l = code[0] + code[1] + code[2];

		TensorData::Ptr td;
		if( m_tensordata.find(component_l) == m_tensordata.end() )
		{
			td = std::make_shared<TensorData>();
			m_tensordata[component_l] = td;
			td->tensor_components.resize(tensor::numComponents(component_l));
			td->tensor = tensor::Tensor<complex>( td->tensor_components.data(), component_l );
		}else
			td = m_tensordata[component_l];

		if( (component_l != m_l) && (component_l != m_l-2))
		{
		//if( (component_l != m_l) )
			std::cout << "error: component rank is not of type l-2=" << m_l-2 << " or l=" << m_l << " code=" << code.toString() << std::endl;
			//return;
		}


		std::vector<V3i>& component_codes = g_component_codes[component_l];

		int numComponents = td->tensor_components.size();
		for( int c=0;c<numComponents;++c )
		{
			if( component_codes[c] == code )
			{
				td->tensor_components[c] += contribution/double(g_component_count[code]);
			}
		}
	}



	std::map<int, TensorData::Ptr> m_tensordata;

	/*
	tensor::Tensor<complex> m_tensor_l;
	tensor::Tensor<complex> m_tensor_lm2;
	std::vector<complex> m_components_l; // tensor components of rank l tensor
	std::vector<complex> m_components_lm2; // tensor components of rank l-2 tensor
	*/
	int m_l;
	int m_m;
};

std::vector<std::vector<V3i>> ylm_t::g_component_codes;
std::map<V3i, int, ylm_t::compare_V3i> ylm_t::g_component_count;

ylm_t get_components( int l, int m, xp::Expression::Ptr ylm_expr )
{
	ylm_t ylm(l,m);

	xp::Addition::Ptr add = std::dynamic_pointer_cast<xp::Addition>(ylm_expr);
	if(!add)
		return ylm;
	//std::cout << "extracting components\n";



	/*
	// print component codes for debugging
	for( auto& code:component_codes )
	{
		for( int i=0;i<3;++i )
			std::cout << code[i];
		std::cout << std::endl;
	}
	*/

	int numTerms = add->getNumOperands();
	for( int i=0;i<numTerms;++i )
	{
		xp::Multiplication::Ptr mul = std::dynamic_pointer_cast<xp::Multiplication>(add->getOperand(i));
		if(!mul)
			std::cout << "error\n";
		//std::cout << "term " << i << "= " << mul->toLatex() << std::endl;

		// each term contributes to one tensor component
		double component_contribution = 1.0;

		// in each term, we expect l number of variables to occur, from
		// which we can figure out the component
		V3i code(0,0,0);

		bool flag = false;

		int numFactors = mul->getNumOperands();
		for( int j=0;j<numFactors;++j )
		{
			xp::Expression::Ptr factor = mul->getOperand(j);

			//std::cout << "\tfactor " << j << "= " << factor->toLatex() << std::endl;

			if( std::dynamic_pointer_cast<xp::Number>(factor) )
			{
				xp::Number::Ptr n = std::dynamic_pointer_cast<xp::Number>(factor);
				double value = 1.0;
				switch(n->getType())
				{
					case xp::Number::EInteger:value = n->get_int();break;
					case xp::Number::EReal:value = n->get_real();break;
					default:
					{
						std::cout << "unable to handle number type\n";
						throw std::runtime_error("hrhhrhth");
					}break;
				};

				component_contribution *= value;
				//std::cout << "\tgot number  value=" << value << std::endl;
			}

			// check for Clm or almj
			if( std::dynamic_pointer_cast<xp::Index>(factor) )
			{
				double index_contribution = 1.0;
				xp::Index::Ptr index = std::dynamic_pointer_cast<xp::Index>(factor);
				if(index->getBaseName() == "C")
					index_contribution *= Clm_including_csp(index->getExponent(0), index->getExponent(1));
				else
				if(index->getBaseName() == "a")
					index_contribution *= almj(index->getExponent(0), index->getExponent(1), index->getExponent(2));

				component_contribution *= index_contribution;
				//std::cout << "\tgot index " << index->getBaseName() << " value=" << index_contribution << std::endl;
			}

			xp::Variable::Ptr var = std::dynamic_pointer_cast<xp::Variable>(factor);
			xp::Power::Ptr pow = std::dynamic_pointer_cast<xp::Power>(factor);
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
				var = std::dynamic_pointer_cast<xp::Variable>(pow->getOperand(0));
				xp::Number::Ptr exp = std::dynamic_pointer_cast<xp::Number>(pow->getOperand(1));
				if(!var)
					std::cout << "error: power of non variable encountered\n";
				if(!exp)
					std::cout << "error: power with non number exponent encountered\n";
				if(exp->getType() != xp::Number::EInteger)
					std::cout << "error: power with non integer exponent encountered\n";

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
				//...
			}

		} // for each factor of current term

		if(flag)
		{
			std::cout << "warning: variable with negative exponent found contribution=" << component_contribution << std::endl;
			continue;
		}

		complex final_contribution = std::pow( complex(0.0, 1.0), code[1] )*component_contribution;


		ylm.add_contribution_to_component(code, final_contribution);
	}

	return ylm;
}

int main()
{
	int order = 5;

	ylm_t::build(order);

	std::vector<ylm_t> ylm_all;
	///*
	{
		using namespace xp;
		Expression::Ptr Ylm = mul( index(var("C"), var("l"), var("m")), pow(add(var("n_x"), var("in_y")), var("m")),sum( "j", largest_integer(mul(pow(num(2), num(-1)) , add(var("l"), mul(num(-1), var("m"))))), mul( index( var("a"), var("l"), var("m"), var("j")),
																		pow( var("n_z"), add(var("l"), mul( num(-1), var("m")), mul( num(-2), var("j")))))));
		Expression::Ptr Ylm2 = mul( pow(num(-1), var("m")), index(var("C"), var("l"), var("m")), pow(add(var("n_x"), mul(num(-1), var("in_y"))), var("m")),sum( "j", largest_integer(mul(pow(num(2), num(-1)) , add(var("l"), mul(num(-1), var("m"))))), mul( index( var("a"), var("l"), var("m"), var("j")),
																		pow( var("n_z"), add(var("l"), mul( num(-1), var("m")), mul( num(-2), var("j")))))));


		//std::cout << "$$" << toLatex(Ylm)  << "$$" << std::endl;
		//std::cout << "$$" << toLatex(Ylm2)  << "$$" << std::endl;




		for( int l=0;l<order;++l )
			for( int m=-l;m<=l;++m )
			{
				//if( !(l==3 && m==1) )
				//	continue;

				ylm_t ylm(l,m);

				// m>=0
				if( m>= 0 )
				{
					Scope scope;
					scope.m_variables["l"] = num(l);
					scope.m_variables["m"] = num(m);

					Expression::Ptr ylm_expr = rewrite::expand(Ylm->deep_copy(), scope);
					//std::cout << "l=" << l << " m=" << m << std::endl;
					//std::cout << "$$" << toLatex(ylm) << "$$" << std::endl;
					ylm = get_components(l, m, ylm_expr);
				}
				else
				{
					// m<0
					Scope scope;
					scope.m_variables["l"] = num(l);
					scope.m_variables["m"] = num(std::abs(m));

					Expression::Ptr ylm2_expr= rewrite::expand(Ylm2->deep_copy(), scope);
					//std::cout << "l=" << l << " m=" << -m << std::endl;
					//std::cout << "$$" << toLatex(ylm2) << "$$" << std::endl;
					ylm = get_components(l, m, ylm2_expr);
				}

				//std::cout << "y" << l << m << "=\n";
				//for( auto& it:ylm_components )
				//{
				//	std::cout << "\t" << it << std::endl;
				//}
				ylm_all.push_back(ylm);
			}
		//e->print_hierarchy();
		//return 0;
	}


	for(auto it:ylm_all)
		it.print();

	//*/
	/*
	{
		int order = 2;
		for( int l=0;l<order;++l )
		{
			std::cout << "l=" << l <<std::endl;
			for( int m=-l;m<=l;++m )
			{
				std::cout << "\tm=" << m <<std::endl;
				int j_end = int(std::ceil((l-m)/2.0));
				std::cout << "\tjend=" << j_end <<std::endl;
				for( int j=0;j<=j_end;++j )
				{
					std::cout << "\tj=" << j << " l-m-2*j=" << l-m-2*j << std::endl;
					std::cout << "\talmj=" << almj(l,m,j) << std::endl;
				}

				std::cout << std::endl;
			}
		}
		return 0;
	}
	*/

	/*
	// gg
	{
		std::vector<double> x_list = linearSamples(-1.0, 1.0, 50);

		int order = 6;
		for( int l=0;l<order;++l )
			for( int m=0;m<=l;++m )
			{
				//bool debug = (l==2)&&(m==1);
				//if(!debug)
				//	continue;
				bool debug = false;

				std::vector<double> y_list;
				std::vector<double> y2_list;
				for( double x:x_list )
				{
					double theta = safe_acos(x);
					double y = sh::P(l, m, x);

					double y2 = P2(l, m, theta, debug);
					y_list.push_back(y);
					y2_list.push_back(y2);
				}
				{
					std::string filename("test_P_$0_$1.txt");
					filename = replace(filename, "$0", toString(l));
					filename = replace(filename, "$1", toString(m));
					writeSamples( filename, x_list, y_list );
				}
				{
					std::string filename("test_P2_$0_$1.txt");
					filename = replace(filename, "$0", toString(l));
					filename = replace(filename, "$1", toString(m));
					writeSamples( filename, x_list, y2_list );
				}
			}
	}
	*/



	// convert groudtruth results to exr
	/*
	{
		int resx = 512;
		int resy = 256;
		int order = 2;
		for( int l=0;l<order;++l )
		{
			for( int m=-l;m<=l;++m )
			{
				EnvMap map;

				std::string prefix = "c:/projects/epfl/epfl17/python";
				std::string filename("test_$0_$1.txt");
				filename = replace(filename, "$0", toString(l));
				filename = replace(filename, "$1", toString(m));

				std::ifstream in( prefix + "/" + filename );

				for( int j=0;j<resy;++j )
					for( int i=0;i<resx;++i )
					{
						std::string line;
						std::getline(in, line);
						std::vector<std::string> parts;
						splitString(line, parts);
						complex csh( fromString<double>(parts[0]), fromString<double>(parts[1]) );
						map.bitmap().coeffRef(j, i) = Color3f(std::abs(csh.real()), std::abs(csh.imag()), 0.0f);
					}

				{
					std::string filename("test2_cshnp_$0_$1.exr");
					filename = replace(filename, "$0", toString(l));
					filename = replace(filename, "$1", toString(m));
					map.bitmap().saveEXR(filename);
				}
			}
		}
	}
	*/


	// investigate SH functions ---
	///*
	int ylm_index = 0;
	for( int l=0;l<order;++l )
	{
		for( int m=-l;m<=l;++m, ++ylm_index )
		{
			//if( (l != 2) && (m!=0))
			//	continue;
			ylm_t ylm = ylm_all[ylm_index];


			EnvMap map_csh_groundtruth;
			EnvMap map_csh;
			EnvMap map_error;
			//EnvMap map_tensor;


			sh::SphericalFunction<Color3f> ylm_complex = [&](double theta, double phi) -> Color3f
			{
				complex csh = complex_sh(l,m,theta, phi);
				return Color3f(std::abs(csh.real()), std::abs(csh.imag()), 0.0f);
			};

			sh::SphericalFunction<Color3f> ylm_complex_ours = [&](double theta, double phi) -> Color3f
			{
				V3d n = sphericalDirection(theta, phi);
				complex csh = complex_sh(l,m,theta, phi);
				if( (l==2) && (m==0) )
				{
					csh = Clm(2, 0)*almj(2, 0, 0)*n.z()*n.z();// + Clm(2, 0)*almj(2, 0, 1);
				}else
				if( (l==2) && (m==1) )
				{
					csh = Clm(2, 1)*almj(2, 1, 0)*n.x()*n.z() + Clm(2, 1)*almj(2, 1, 0)*n.z()*complex(0.0, n.y());
				}

				return Color3f(std::abs(csh.real()), std::abs(csh.imag()), 0.0f);
			};


			sh::SphericalFunction<Color3f> error = [&](double theta, double phi) -> Color3f
			{
				V3d n = sphericalDirection(theta, phi);
				complex csh = complex_sh(l, m, theta, phi);
				//complex csh3 = complex_sh3(l, m, n);
				//complex csh3 = contract_moment(ylm.m_tensor_l, n);
				complex csh3;

				//if(l>=2)
				//	csh3 += contract_moment(ylm.m_tensor_lm2, n);
				for( auto& it:ylm.m_tensordata )
				{
					csh3 += contract_moment(it.second->tensor, n);
				}

				double error = std::abs(csh-csh3);
				//double error = std::abs(csh-csh4);

				return Color3f(error);
			};

			rasterizeSphericalFunction(map_error, error);
			rasterizeSphericalFunction(map_csh_groundtruth, ylm_complex);
			rasterizeSphericalFunction(map_csh, ylm_complex_ours);

			{
				std::string filename("error_$0_$1.exr");
				filename = replace(filename, "$0", toString(l));
				filename = replace(filename, "$1", toString(m));
				map_error.bitmap().saveEXR(filename);
			}
			{
				std::string filename("csh_$0_$1.exr");
				filename = replace(filename, "$0", toString(l));
				filename = replace(filename, "$1", toString(m));
				map_csh_groundtruth.bitmap().saveEXR(filename);
			}
			{
				std::string filename("csh_ours_$0_$1.exr");
				filename = replace(filename, "$0", toString(l));
				filename = replace(filename, "$1", toString(m));
				map_csh.bitmap().saveEXR(filename);
			}
		}
	}
	//*/







	/*
	EnvMap map("envmap.exr");

	sh::SphericalFunction<Color3f> fun = [&](double theta, double phi) -> Color3f
	{
		return map.eval(theta, phi);
	};

	writeSphericalFunction( "groundtruth.bgeo", fun );

	// project
	int order = 30;
	int numSamples = 50000;
	std::unique_ptr<std::vector<Color3f>> sh_coeffs;
	sh_coeffs = sh::project<Color3f>(order, fun, numSamples);

	// reconstruction
	for( int l=0;l<=order;++l )
	{
		sh::SphericalFunction<Color3f> reconstruction = [&](double theta, double phi) -> Color3f
		{
			return sh::evalSum(l, *sh_coeffs.get(), theta, phi);
		};

		std::string filename("testfit_reconstruction_$0.bgeo");
		filename = replace(filename, "$0", toString(l));

		writeSphericalFunction( filename, reconstruction );
	}
	*/


	/*
	int L = 6;
	for( int l=0;l<L;++l )
	{
		for( int m=-l;m<=l;++m )
		{

			// visualization for houdini
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
					//Color3f col = map.eval(theta_phi.x(), theta_phi.y());

					float sh = EvalSHSlow(l, m, theta, phi);

					//Color3f col(sh);
					//cdAttr->set<houio::math::V3f>( i, houio::math::V3f(col.r(), col.g(), col.b()) );

					pAttr->set<houio::math::V3f>( i, p*sh );
				}
				//geo->setAttr("Cd", cdAttr);
				std::string filename("test_$0_$1.bgeo");
				filename = replace(filename, "$0", toString(l));
				filename = replace(filename, "$1", toString(m));
				houio::HouGeoIO::xport( filename, geo);
			}
		}
	}
	*/




	return 0;


	//std::string basePath = "c:/projects/visus/data";
	std::string basePath = ".";
	//std::string outpath = basePath + "/noisereduction/fitting4";
	std::string outpath = basePath + "";




	// output image ---
	V2i res = V2i(512, 512);
	//V2i res = V2i(2048, 2048);
	Bitmap image_color(res);
	Bitmap image_transmittance(res);


	// scene elements ---
	//Volume::Ptr volume = volumes::C60();
	Volume::Ptr volume = volumes::nebulae();
	//Light::Ptr light = lights::point( volume->getBound().getCenter() + P3d(0.0, volume->getBound().getExtents().y()*0.6, 0.0) );
	Light::Ptr light = lights::directional( V3d(0.0, -1.0, 0.0), volume->getBound() );

	Scene scene;
	scene.bound = volume->getBound();
	scene.id = "nebulae";
	scene.volume = volume.get();
	scene.light = light.get();


	//Camera::Ptr camera = createView( V3d(0.0, 0.0, 1.0), scene.bound, res );



	/*
	// generate pn cache -----------------------
	{
		Integrator::Ptr integrator = integrators::volume_path_tracer();
		scene.integrator = integrator.get();
		std::string filename = outpath + "/nebulae.moments";

		PNCache cache;
		int numMoments = 4;
		int numSamples = 100;
		int res = 128;
		cache.generate( outpath + "/nebulae.moments", &scene, numMoments, numSamples, res );
		return;
	}
	*/

	// pn analysis and debugging ---
	{
		Integrator::Ptr integrator = integrators::dummy();
		scene.integrator = integrator.get();

		PNCache cache;
		cache.m_override = true;
		cache.m_doZeroScattering = false;
		int numMoments = 4;
		int numSamples = 1000000;
		int res = 1;
		cache.generate( outpath + "/analysis", &scene, numMoments, numSamples, res );


		Wedge wedge;
		std::vector<int> moment_values = {0, 1, 2, 3, 4};
		wedge.addParm("moment", moment_values);
		wedge.build();

		std::cout << "4pi=" << 4.0*M_PI << std::endl;
		std::cout << "4pi/3=" << 4.0*M_PI/3.0 << std::endl;

		std::vector<Wedge::Iteration> iterations = wedge.iterations();
		for( auto it : iterations )
		{
			int moment_index = it.getInt("moment");
			std::cout << "moment=" << moment_index << std::endl;

			tensor::Tensor<double> moment_tensor = cache.getMoment( 0, 0, 0, moment_index );

			for(auto component = moment_tensor.begin(), end=moment_tensor.end(); component!=end;++component)
			{
				std::cout << "\t" << component.index_str() << "=" << component.value() << std::endl;
				/*
				houio::Geometry::Ptr geo = houio::Geometry::createSphere(50, 50, 1.0);
				houio::Attribute::Ptr pattr = geo->getAttr("P");
				int numPoints = pattr->numElements();
				for(int i=0;i<numPoints;++i)
				{
					houio::math::V3f d = pattr->get<houio::math::V3f>(i).normalized();
					d *= std::abs(component.weight(V3d(d.x, d.y, d.z)));
					pattr->set<houio::math::V3f>(i, d);
				}
				houio::HouGeoIO::xport(it.expand_value( outpath+"/analysis_$0_" + component.index_str() + ".bgeo"), geo);
				*/
			}
		}
	}


	Eigen::Affine3d cameraToWorld;
	cameraToWorld = Eigen::Translation3d(V3d(0.0, 0.0, -2.5));

	//Camera::Ptr camera = std::make_shared<OrthographicCamera>(res.x(), res.y(), 2.0, 2.0, 1e-4f, 5000.0);
	Camera::Ptr camera = std::make_shared<PerspectiveCamera>(res.x(), res.y() );
	//camera->setCameraToWorld( Eigen::Affine3d(lookAt<double>(targetWS-viewDir*b, targetWS)) );
	camera->setCameraToWorld(Transformd(cameraToWorld));

	scene.camera = camera.get();


	/*
	// RENDERING -----------------------------------
	{
		//Integrator::Ptr integrator = integrators::raymarcher(0.005);
		PNCache::Ptr cache = std::make_shared<PNCache>(outpath + "/nebulae.moments");
		cache->eval(P3d(0.0f, 0.0, 0.0), normalize(V3d(1.0, 1.0, 1.0)), true);

		//Bitmap::Ptr pixel_estimates = readImage(outpath + "/nebulae_pixel_estimate.exr");
		Integrator::Ptr integrator = integrators::cache_raymarcher(0.005, cache);
		//Integrator::Ptr integrator = integrators::adrrs_volume_path_tracer(pixel_estimates, fluence_field);
		scene.integrator = integrator.get();

		int numSamples = 1;

		// prepare render info ---
		RenderTaskInfo::g.scene = &scene;
		RenderTaskInfo::g.image = &image_color;
		RenderTaskInfo::g.image_transmittance = &image_transmittance;
		RenderTaskInfo::g.crop_window = Box2i( V2i(0,0), res );
		//RenderTaskInfo::g.debug_pixel = P2i(318, 209);
		//RenderTaskInfo::g.debug_pixel = P2i(256, 256);
		//RenderTaskInfo::g.debug_pixel = P2i(340, 340);


		// execute multithreaded render ---
		Terminator terminator(numSamples);
		std::cout << "rendering image..."<< std::endl;std::flush(std::cout);
		runGenericTasks<RenderTaskInfo>( render_volume, terminator, ThreadPool::getNumSystemCores() );


		// save results ---
		flip(image_color);
		flip(image_transmittance);
		image_color.saveEXR( outpath + "/" + scene.id + "_color.exr" );
		std::string transmittance_exr = outpath + "/test_transmittance.exr";
		image_transmittance.saveEXR(transmittance_exr);
	}
	*/




	return 0;
}


namespace xp
{
	Expression::Ptr num( int value )
	{
		return std::make_shared<Number>(value);
	}


	Expression::Ptr var( const std::string& name )
	{
		return std::make_shared<Variable>(name);
	}

	Expression::Ptr add( Expression::Ptr a, Expression::Ptr b )
	{
		Addition::Ptr add = std::make_shared<Addition>();
		add->addOperand(a);
		add->addOperand(b);
		return add;
	}

	Expression::Ptr add( Expression::Ptr a, Expression::Ptr b, Expression::Ptr c )
	{
		Addition::Ptr add = std::make_shared<Addition>();
		add->addOperand(a);
		add->addOperand(b);
		add->addOperand(c);
		return add;
	}

	Expression::Ptr sum( const std::string& index, Expression::Ptr end, Expression::Ptr body )
	{
		return std::make_shared<Sum>(index, end, body);
	}

	Expression::Ptr mul( Expression::Ptr a, Expression::Ptr b )
	{
		// if one of the two factors is a one, we will just return the other side
		Number::Ptr num_a = std::dynamic_pointer_cast<Number>(a);
		Number::Ptr num_b = std::dynamic_pointer_cast<Number>(b);
		if(num_a && (num_a->getType() == Number::EInteger)&& (num_a->get_int() == 1) )
			return b;
		if(num_b && (num_b->getType() == Number::EInteger)&& (num_b->get_int() == 1) )
			return a;


		Multiplication::Ptr mul = std::make_shared<Multiplication>();
		mul->addOperand(a);
		mul->addOperand(b);
		return mul;
	}
	Expression::Ptr mul( Expression::Ptr a, Expression::Ptr b, Expression::Ptr c )
	{
		Multiplication::Ptr mul = std::make_shared<Multiplication>();
		mul->addOperand(a);
		mul->addOperand(b);
		mul->addOperand(c);
		return mul;
	}
	Expression::Ptr mul( Expression::Ptr a, Expression::Ptr b, Expression::Ptr c, Expression::Ptr d )
	{
		Multiplication::Ptr mul = std::make_shared<Multiplication>();
		mul->addOperand(a);
		mul->addOperand(b);
		mul->addOperand(c);
		mul->addOperand(d);
		return mul;
	}

	Expression::Ptr pow( Expression::Ptr base, Expression::Ptr exp )
	{
		Number::Ptr num_exp = std::dynamic_pointer_cast<Number>(exp);
		// if exponent is 0, we return 1
		if( num_exp && (num_exp->getType() == Number::EInteger) && (num_exp->get_int() == 0) )
			return num(1);
		// if exponent is 1, we return base
		if( num_exp && (num_exp->getType() == Number::EInteger) && (num_exp->get_int() == 1) )
			return base;


		return std::make_shared<Power>(base, exp);
	}

	Expression::Ptr index( Expression::Ptr base, Expression::Ptr a, Expression::Ptr b)
	{
		Index::Ptr idx = std::make_shared<Index>(base);
		idx->addIndex(a);
		idx->addIndex(b);
		return idx;
	}

	Expression::Ptr index( Expression::Ptr base, Expression::Ptr a, Expression::Ptr b, Expression::Ptr c )
	{
		Index::Ptr idx = std::make_shared<Index>(base);
		idx->addIndex(a);
		idx->addIndex(b);
		idx->addIndex(c);
		return idx;
	}

	Expression::Ptr largest_integer( Expression::Ptr expr )
	{
		return std::make_shared<LargestInteger>(expr);
	}

	std::string toLatex(Expression::Ptr expr)
	{
		return expr->toLatex();
	}

	namespace rewrite
	{
		Expression::Ptr expand( Expression::Ptr expr, Scope& scope )
		{
			Expression::Ptr result = expr;

			// we replace all variables
			for( auto& it : scope.m_variables )
			{
				std::string name = it.first;
				Expression::Ptr replacement = it.second;
				replace_variable(result, name, replacement);
			}

			// and fold constants
			result = fold_constants(result);

			// expand sum, binomial coefficients
			result = expand(result);

			// fold constants again
			result = fold_constants(result);

			// expand sum, binomial coefficients
			result = expand(result);

			return result;
		}


		Expression::Ptr expand_sum(Sum::Ptr sum)
		{
			Addition::Ptr add = std::make_shared<Addition>();
			Number::Ptr end = std::dynamic_pointer_cast<Number>(sum->getOperand(0));
			Expression::Ptr body = sum->getOperand(1);
			if(!end && end->getType() != Number::EInteger)
			{
				std::cout << "expand_sum: error: end index is no number or not an integer\n";
				throw std::runtime_error("afasfasf");
			}
			int numIterations = end->get_int();
			for(int i=0;i<=numIterations;++i)
			{
				Expression::Ptr expr = body->deep_copy();
				rewrite::replace_variable( expr, sum->getIndexName(), num(i));
				add->addOperand(expr);
			}
			return add;
		}

		Expression::Ptr expand_binomial( Power::Ptr p )
		{
			Addition::Ptr base = std::dynamic_pointer_cast<Addition>( p->getOperand(0) );
			Number::Ptr exp = std::dynamic_pointer_cast<Number>( p->getOperand(1) );

			if(base && exp && base->getNumOperands()==2 && exp->getType()==Number::EInteger)
			{
				Expression::Ptr x = base->getOperand(0);
				Expression::Ptr y = base->getOperand(1);
				Addition::Ptr add = std::make_shared<Addition>();

				int n = exp->get_int();

				if(n==0)
					return num(1);

				for( int k=0;k<=n;++k )
				{
					int exp_x = n-k;
					int exp_y = k;

					Expression::Ptr pow_x;
					Expression::Ptr pow_y;

					if( exp_x > 0 )
					{
						if(exp_x>1)
							pow_x = pow(x->deep_copy(), num(exp_x));
						else
							pow_x = x->deep_copy();
					}
					if( exp_y > 0 )
					{
						if(exp_y>1)
							pow_y = pow(y->deep_copy(), num(exp_y));
						else
							pow_y = y->deep_copy();
					}

					int coeff = sh::factorial(n)/(sh::factorial(k)*sh::factorial(n-k));

					if( pow_x && pow_y )
					{
						if(coeff > 1)
							add->addOperand( mul(num(coeff), pow_x, pow_y) );
						else
							add->addOperand( mul(pow_x, pow_y) );
					}else
					if( pow_x )
					{
						if(coeff > 1)
							add->addOperand( mul(num(coeff), pow_x) );
						else
							add->addOperand( pow_x );
					}
					else
					if( pow_y )
					{
						if(coeff > 1)
							add->addOperand( mul(num(coeff), pow_y) );
						else
							add->addOperand( pow_y );
					}else
					{
						std::cout << "expand_binomial: error: \n";
						throw std::runtime_error("asfsfasf");
					}

				} // for all combinations

				return add;
			}

			return p;
		}

		Expression::Ptr expand_distributive( Multiplication::Ptr m )
		{
			// find addition
			Addition::Ptr add;
			Multiplication::Ptr factor = std::make_shared<Multiplication>();
			{
				int numOperands = m->getNumOperands();
				for( int i=0;i<numOperands;++i )
				{
					add = std::dynamic_pointer_cast<Addition>(m->getOperand(i));
					if(add)
						break;
				}
				if(!add)
					// no addition term
					return m;

				for( int i=0;i<numOperands;++i )
				{
					Expression::Ptr operand = m->getOperand(i);
					if(operand != add)
						factor->addOperand(operand);
				}
			}


			Addition::Ptr result = std::make_shared<Addition>();

			int numOperands = add->getNumOperands();
			for( int i=0;i<numOperands;++i )
			{
				Expression::Ptr term = mul(add->getOperand(i)->deep_copy(), factor->deep_copy());
				result->addOperand( expand(term) );
			}

			return expand(result);
		}

		Expression::Ptr expand(Expression::Ptr expr)
		{
			// first we expand all childs...
			if(std::dynamic_pointer_cast<Operator>(expr))
			{
				Operator::Ptr op = std::dynamic_pointer_cast<Operator>(expr);
				int numOperands = op->getNumOperands();
				for( int i=0;i<numOperands;++i )
					op->setOperand(i, expand(op->getOperand(i)));
			}

			// now handle concrete expansion rules..
			// expand sumation signs
			if(std::dynamic_pointer_cast<Sum>(expr))
				return expand_sum(std::dynamic_pointer_cast<Sum>(expr));
			else
			if(std::dynamic_pointer_cast<Power>(expr))
			{
				Power::Ptr p = std::dynamic_pointer_cast<Power>(expr);

				// apply binomial expansion
				{
					Addition::Ptr base = std::dynamic_pointer_cast<Addition>( p->getOperand(0) );
					Number::Ptr exp = std::dynamic_pointer_cast<Number>( p->getOperand(1) );

					if(base && exp && base->getNumOperands()==2 && exp->getType()==Number::EInteger)
						return expand_binomial(p);
				}
				// apply exponent shortcuts
				{
					Expression::Ptr base = p->getOperand(0);
					Number::Ptr exp = std::dynamic_pointer_cast<Number>( p->getOperand(1) );
					if( exp && (exp->getType()==Number::EInteger) && (exp->get_int()==0) )
						return num(1);
					if( exp && (exp->getType()==Number::EInteger) && (exp->get_int()==1) )
						return p->getOperand(0);
				}
				// distribute exponent over multiplication
				{
					Multiplication::Ptr mul = std::dynamic_pointer_cast<Multiplication>(p->getOperand(0));
					Number::Ptr exp = std::dynamic_pointer_cast<Number>( p->getOperand(1) );
					if(mul)
					{
						int numFactors = mul->getNumOperands();
						Multiplication::Ptr mul_new = std::make_shared<Multiplication>();
						for( int i=0;i<numFactors;++i )
							mul_new->addOperand(pow(mul->getOperand(i), exp->deep_copy()));
						return expand(mul_new);
					}
				}
			}else
			// apply distributive law
			if(std::dynamic_pointer_cast<Multiplication>(expr))
			{
				Multiplication::Ptr mul = std::dynamic_pointer_cast<Multiplication>(expr);

				int numOperands = mul->getNumOperands();
				for( int i=0;i<numOperands;++i )
				{
					if(std::dynamic_pointer_cast<Addition>(mul->getOperand(i)))
						return expand_distributive(mul);
				}

				// flatten nested multiplcations
				{
					Multiplication::Ptr mul_new = std::make_shared<Multiplication>();
					bool got_new = false;
					int numOperands = mul->getNumOperands();
					for( int i=0;i<numOperands;++i )
					{
						Expression::Ptr term = mul->getOperand(i);
						Multiplication::Ptr term_mul = std::dynamic_pointer_cast<Multiplication>(term);
						if(term_mul)
						{
							got_new = true;
							int num_nested = term_mul->getNumOperands();
							for(int j=0;j<num_nested;++j)
								mul_new->addOperand(term_mul->getOperand(j));
						}else
						{
							Number::Ptr num = std::dynamic_pointer_cast<Number>(term);
							// we ignore constants of one
							if( (num && (num->getType() == Number::EInteger) && (num->get_int() == 1))||
								(num && (num->getType() == Number::EReal) && (num->get_real() == 1.0)))
							{
								got_new = true;
							}else
								mul_new->addOperand(term);
						}
					}
					if(got_new)
					{
						if( mul_new->getNumOperands() == 1 )
							return mul_new->getOperand(0);
						else
							return mul_new;
					}
				}
			}else
			// flatten nested additions
			if(std::dynamic_pointer_cast<Addition>(expr))
			{
				Addition::Ptr add = std::dynamic_pointer_cast<Addition>(expr);
				Addition::Ptr add_new = std::make_shared<Addition>();
				bool got_nested = false;
				int numOperands = add->getNumOperands();
				for( int i=0;i<numOperands;++i )
				{
					Expression::Ptr term = add->getOperand(i);
					Addition::Ptr term_add = std::dynamic_pointer_cast<Addition>(term);
					if(term_add)
					{
						got_nested = true;
						int num_nested = term_add->getNumOperands();
						for(int j=0;j<num_nested;++j)
							add_new->addOperand(term_add->getOperand(j));
					}else
						add_new->addOperand(term);
				}
				if(got_nested)
					return add_new;
			}



			return expr;
		}

		void replace_variable( Expression::Ptr expr, const std::string& variable, Expression::Ptr replacement )
		{
			if( std::dynamic_pointer_cast<Operator>(expr) )
			{
				Operator::Ptr op = std::dynamic_pointer_cast<Operator>(expr);
				int numChilds = op->getNumOperands();
				for( int i=0;i<numChilds;++i )
				{
					Expression::Ptr child = op->getOperand(i);

					Variable::Ptr var = std::dynamic_pointer_cast<Variable>(child);
					if(var)
					{
						if(var->getName() == variable)
							op->setOperand(i, replacement);
					}else
						replace_variable(child, variable, replacement);
				}
			}
		}

		Expression::Ptr fold_constants( Addition::Ptr add )
		{
			std::vector<Number::Ptr> number_childs;
			std::vector<Expression::Ptr> other_childs;

			int numChilds = add->getNumOperands();
			for( int i=0;i<numChilds;++i )
			{
				Expression::Ptr child = fold_constants( add->getOperand(i) );

				if( std::dynamic_pointer_cast<Number>(child) )
					number_childs.push_back(std::dynamic_pointer_cast<Number>(child));
				else
					other_childs.push_back(child);
			}

			// fold all the numbers we have
			Number::Ptr num;
			if( number_childs.size() > 0 )
				num = std::dynamic_pointer_cast<Number>(number_childs[0]->deep_copy());

			for( int i=1;i<number_childs.size();++i )
				*num += *number_childs[i];

			if( number_childs.size() == numChilds )
			{
				// all childs are numbers, we can return a pure number
				return num;
			}else
			{
				// only sum (or even no) childs are numbers, fold all the numbers we have
				Addition::Ptr add_folded = std::make_shared<Addition>();
				if(num)
					add_folded->addOperand(num);
				for( auto expr:other_childs )
					add_folded->addOperand(expr);
				return add_folded;
			}
		}

		Expression::Ptr fold_constants( Multiplication::Ptr mul )
		{
			std::vector<Number::Ptr> number_childs;
			std::vector<Expression::Ptr> other_childs;

			int numChilds = mul->getNumOperands();
			for( int i=0;i<numChilds;++i )
			{
				Expression::Ptr child = fold_constants( mul->getOperand(i) );

				if( std::dynamic_pointer_cast<Number>(child) )
					number_childs.push_back(std::dynamic_pointer_cast<Number>(child));
				else
					other_childs.push_back(child);
			}

			// fold all the numbers we have
			Number::Ptr num;
			if( number_childs.size() > 0 )
				num = std::dynamic_pointer_cast<Number>(number_childs[0]->deep_copy());

			for( int i=1;i<number_childs.size();++i )
			{
				*num *= *number_childs[i];
			}

			if( number_childs.size() == numChilds )
			{
				// all childs are numbers, we can return a pure number
				return num;
			}else
			{
				// only sum (or even no) childs are numbers, fold all the numbers we have
				Multiplication::Ptr mul_folded = std::make_shared<Multiplication>();
				if(num)
					mul_folded->addOperand(num);
				for( auto expr:other_childs )
					mul_folded->addOperand(expr);
				return mul_folded;
			}
		}

		Expression::Ptr fold_constants( Power::Ptr pow )
		{
			pow->setOperand(0, fold_constants(pow->getOperand(0)));
			pow->setOperand(1, fold_constants(pow->getOperand(1)));

			Number::Ptr base = std::dynamic_pointer_cast<Number>(pow->getOperand(0));
			Number::Ptr exp = std::dynamic_pointer_cast<Number>(pow->getOperand(1));
			if( base && exp )
			{
				base = std::dynamic_pointer_cast<Number>(base->deep_copy());
				base->pow(*exp);
				return base;
			}
			return pow;
		}

		Expression::Ptr fold_constants( LargestInteger::Ptr li )
		{
			li->setOperand(0, fold_constants(li->getOperand(0)));

			Number::Ptr base = std::dynamic_pointer_cast<Number>(li->getOperand(0));
			if(base)
			{
				if( base->getType() == Number::EInteger )
					return num(base->get_int());
				else
				if( base->getType() == Number::EReal )
					return num(int(std::ceil( base->get_real() )));
			}
			return li;
		}

		Expression::Ptr fold_constants( Expression::Ptr expr )
		{
			if( std::dynamic_pointer_cast<Operator>(expr) )
			{
				Operator::Ptr op = std::dynamic_pointer_cast<Operator>(expr);

				if( std::dynamic_pointer_cast<Multiplication>(op) )
					return fold_constants( std::dynamic_pointer_cast<Multiplication>(op));
				else
				if( std::dynamic_pointer_cast<Addition>(op) )
					return fold_constants( std::dynamic_pointer_cast<Addition>(op));
				else
				if( std::dynamic_pointer_cast<Power>(op) )
					return fold_constants( std::dynamic_pointer_cast<Power>(op));
				else
				if( std::dynamic_pointer_cast<LargestInteger>(op) )
					return fold_constants( std::dynamic_pointer_cast<LargestInteger>(op));

				int numChilds = op->getNumOperands();
				for( int i=0;i<numChilds;++i )
				{
					Expression::Ptr child = op->getOperand(i);
					op->setOperand(i, fold_constants(child));
				}
			}

			return expr;
		}
	}
}
