#include <util/cas.h>


namespace cas
{
	void Expression::print_hierarchy( int indent )const
	{
		for(int i=0;i<indent;++i)std::cout<<"\t";
		std::cout << "Expression" << std::endl;
	}


	namespace rewrite
	{
		Expression::Ptr expand(Expression::Ptr e);
		void replace_variable( Expression::Ptr expr, const std::string& variable, Expression::Ptr replacement );
		Expression::Ptr fold_constants( Expression::Ptr expr );
	}


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


	Expression::Ptr expand( Expression::Ptr expr, Scope& scope )
	{
		Expression::Ptr result = expr;

		// we replace all variables
		for( auto& it : scope.m_variables )
		{
			std::string name = it.first;
			Expression::Ptr replacement = it.second;
			rewrite::replace_variable(result, name, replacement);
		}

		// and fold constants
		result = rewrite::fold_constants(result);

		// expand sum, binomial coefficients
		result = rewrite::expand(result);

		// fold constants again
		result = rewrite::fold_constants(result);

		// expand sum, binomial coefficients
		result = rewrite::expand(result);

		return result;
	}

	namespace rewrite
	{
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

		double factorial(int x)
		{
			double s = 1.0;
			for (int n = 2; n <= x; ++n)
				s *= n;
			return s;
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

					int coeff = factorial(n)/(factorial(k)*factorial(n-k));

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
} // namespace cas
