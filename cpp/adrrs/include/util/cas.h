/*
  cas stands for Computer Algebra System.
  It is a very crude and simple framework, which allows to work with symbolic
  algebraic expressions.

  With this framework, you can create algebraic expressions, including addition,
  multiplication and power and run a couple of manipulations on them.

  Simple example:
  Consider the algebraic term (a+b)^m
  Now consider that we want to have a program, which is able to list all the
  coefficients of the various terms, which result from expanding this equation for
  a given m. In this example, it is easy to write down all the terms using
  the binomial coefficient. However, once additional factors, such as a sum sign
  (\Sigma) with an index j running from 0 to e.g. (l-m)/2 get added, it can become
  quite involved to come up with an algorithm for listing all the terms.
  This framework is basically an over-engineered solution to this problem. You give
  it the equation and values for some variables, and it can expand it for you.

  For the example above you can do:

  using namespace cas;
  // this is how you create expressions
  Expression::Ptr expr = pow(add(var("a"), var("b")), var("m"));

  // now we want to specify m, which could be another expression or some concrete value:
  Scope scope;
  scope.m_variables["m"] = num(7);

  // now we want to expand our expression.
  // cas has the expand function, which tries a couple of rules for flattening
  // the equation. Expanding binomial coefficients is one of them...
  Expression::Ptr expanded = expand(expr, scope);

  // every expression can be turned into latex for inspection
  std::cout << "$$" << expanded->toLatex() << "$$" << std::endl;

  // this will be printed:
  // $$a^{7}+7a^{6}b+21a^{5}b^{2}+35a^{4}b^{3}+35a^{3}b^{4}+21a^{2}b^{5}+7ab^{6}+b^{7}$$

*/
#pragma once

#include <memory>
#include <string>
#include <sstream>
#include <iostream>
#include <map>
#include <vector>



namespace cas
{
	// This is the fudamental datastructure which represents an expression in cas.
	struct Expression
	{
		typedef std::shared_ptr<Expression> Ptr;
		virtual std::string toLatex()const=0;
		virtual void print_hierarchy( int indent = 0 )const;
		virtual Ptr deep_copy()const=0;
	};

	// the following api calls allow to build complex expression trees
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

	// the expand function takes additional variable definitions as second input
	struct Scope
	{
		std::map<std::string, Expression::Ptr> m_variables;
	};

	// expand will run a set of rules for manipulating algebraic expressions,
	// such as expanding binomial coefficients, constant folding, distributive law, etc...
	Expression::Ptr expand( Expression::Ptr expr, Scope& scope );

	// what follows are different type of expressions, such as variables, operators etc.
	// the reason why this is exposed to the use is to allow him to query expressions
	// for its components. This is useful when you want to use the expanded result.
	// For example if you want to query all the coefficients for the different terms
	// of a summation...

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

	template<class T>
	inline std::string toString(const T& t)
	{
		std::ostringstream stream;
		stream << t;
		return stream.str();
	}

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
} // namespace cas
