import traceback

class Expression(object):
	def __init__(self, children = []):
		if isinstance(children, tuple):
			print("!!!!!!!!!!!!!!!!!!")
			raise ValueError('A very specific bad thing happened')

		self.children = children

	def toLatex(self):
		return 

	def indent_string( self, indent ):
		return ''.join( ['\t' for i in range(indent)] )

	def toHierarchy( self, indent = 0 ):
		result = self.indent_string(indent) + "{}\n".format(self.__class__.__name__)
		for child in self.children:
			result += child.toHierarchy(indent+1)
		return result


	def numChildren(self):
		return len(self.children)

	def getChildren(self, index):
		return self.children[index]

	def getAllChildren(self):
		return self.children

	def setChildren(self, index, new_child):
		self.children[index] = new_child

	def depends_on(self, variable):
		for i in range(self.numChildren()):
			if self.getChildren(i).depends_on(variable):
				return True
		return False

	def __eq__(self, other):
		return False



class Number(Expression):
	def __init__( self, value ):
		super().__init__()
		self.value = value

	def toLatex(self):
		return str(self.value)

	def getValue(self):
		return self.value

	def toHierarchy( self, indent = 0 ):
		return self.indent_string(indent) + "Number {}\n".format(self.getValue())

	def deep_copy(self):
		return Number(self.getValue())

class Variable(Expression):
	def __init__( self, symbol ):
		super().__init__()
		self.symbol = symbol

	def getSymbol(self):
		return self.symbol

	def depends_on(self, variable):
		return self.symbol == variable.symbol

	def __eq__(self, other):
		if other.__class__ == Variable:
			if other.getSymbol() == self.getSymbol():
				return True
		return False

	def deep_copy(self):
		return Variable(self.getSymbol())

	def toLatex(self):
		return str(self.symbol)

	def toHierarchy( self, indent = 0 ):
		return self.indent_string(indent) + "Variable {}\n".format(self.symbol)


class TensorComponent(Variable):
	def __init__( self, tensor_symbol, component_symbol ):
		super().__init__(tensor_symbol)
		self.component_symbol = component_symbol
	def deep_copy(self):
		return TensorComponent(self.getSymbol(), self.component_symbol)
	def toLatex(self):
		return "{}_{{{}}}".format(self.getSymbol(), self.component_symbol)


		
class Function( Expression ):
	def __init__( self, symbol, arguments ):
		super().__init__(arguments)
		self.symbol = symbol
		self.latexArgumentPositions = [0 for arg in arguments]

		for arg in arguments:
			if isinstance(arg, str):
				print("===================================")
				raise ValueError()

	def getSymbol(self):
		return self.symbol

	def numArguments(self):
		return self.numChildren()

	def getArgument(self, index):
		return self.getChildren(index)

	def getArguments(self):
		return self.getAllChildren()

	def setLatexArgumentPosition( self, index, level ):
		# this method specifies, how function arguments are rendered in latex
		# default is, that all arguments come in parentheses after the function symbol
		# however, arguments could also be rendered as sub- or superscripts
		# level=-1: argument is a subscript
		# level=0: render as default
		# level=1: argument is a superscript
		self.latexArgumentPositions[index] = level

	def setAllSuperScripts(self):
		self.latexArgumentPositions = [1 for arg in self.getArguments()]

	def deep_copy(self):
		cp = Function(self.getSymbol(), [arg.deep_copy() for arg in self.getArguments()])
		cp.latexArgumentPositions = [level for level in self.latexArgumentPositions]
		return cp

	def toLatex(self):
		string = self.getSymbol()
		if self.numArguments() == 0:
			return string
		else:
			superscript_args = []
			normal_args = []
			subscript_args = []
			for i in range(self.numArguments()):
				if self.latexArgumentPositions[i] == -1:
					subscript_args.append(i)
				elif self.latexArgumentPositions[i] == 0:
					normal_args.append(i)
				elif self.latexArgumentPositions[i] == 1:
					superscript_args.append(i)
				else:
					raise ValueError("asdasdasdasd")

			# render superscript arguments
			if superscript_args:
				string += "^{{"
				for i in range(len(superscript_args)):
					string += self.getArgument(superscript_args[i]).toLatex()
					if i < len(superscript_args)-1:
						string+=","
				string += "}}"
			# render subscript args
			if subscript_args:
				string += "_{{"
				for i in range(len(subscript_args)):
					string += self.getArgument(subscript_args[i]).toLatex()
					if i < len(subscript_args)-1:
						string+=","
				string += "}}"
			# render normal arguments
			if normal_args:
				string += "\\left ("
				for i in range(len(normal_args)):
					string += self.getArgument(normal_args[i]).toLatex()
					if i < len(normal_args)-1:
						string+=","
				string += "\\right )"
			return string

	def __eq__(self, other):
		# two functions are considered the same, if they have the
		# same symbol and the same number of arguments
		# later we might take physical dimension into account
		if other.__class__ == Function:
			if other.getSymbol() == self.getSymbol() and other.numArguments() == self.numArguments():
				return True
		return False
		



class Kronecker(Function):
	def __init__(self, index_i, index_j):
		super().__init__( "\\delta", [index_i, index_j])

	def getFirstIndex(self):
		return self.getArgument(0)

	def getSecondIndex(self):
		return self.getArgument(1)

	def deep_copy(self):
		return Kronecker(self.getFirstIndex().deep_copy(), self.getSecondIndex().deep_copy())

	def toLatex(self):
		return "\\delta_{{{},{}}}".format(self.getArgument(0).toLatex(), self.getArgument(1).toLatex())

class SHBasis(Function):
	def __init__(self, l, m, omega, conjugate_complex=False):
		super().__init__( "Y", [l, m, omega])
		self.conjugate_complex = conjugate_complex

	def get_l(self):
		return self.getChildren(0)
	def get_m(self):
		return self.getChildren(1)
	def get_direction_argument(self):
		return self.getChildren(2)
	def is_conjugate_complex(self):
		return self.conjugate_complex

	def deep_copy(self):
		l = self.get_l()
		m = self.get_m()
		omega = self.get_direction_argument()

		if isinstance(l, str):
			print("YYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYY")

		return SHBasis( self.get_l().deep_copy(), self.get_m().deep_copy(), self.get_direction_argument().deep_copy(), self.is_conjugate_complex() )

	def toLatex(self):
		l = self.get_l()
		m = self.get_m()
		omega = self.get_direction_argument()
		#print([l,m,omega])
		if self.is_conjugate_complex():
			return "\overline{Y^{" + l.toLatex() + ", " + m.toLatex() + "}}(" + omega.toLatex() + ")"
		else:
			return "Y^{" + l.toLatex() + ", " + m.toLatex() + "}(" + omega.toLatex() + ")"


class SHCoefficient(Function):
	def __init__(self, symbol, l, m, x):
		super().__init__( symbol, [l, m, x])

		if isinstance(l, str):
			raise ValueError()

		self.setLatexArgumentPosition(0, 1)
		self.setLatexArgumentPosition(1, 1)
	def get_l(self):
		return self.getChildren(0)
	def get_m(self):
		return self.getChildren(1)
	def get_position_argument(self):
		return self.getChildren(2)
	def deep_copy(self):
		return SHCoefficient( self.getSymbol(), self.get_l().deep_copy(), self.get_m().deep_copy(), self.get_position_argument().deep_copy() )



class Operator( Expression ):
	def __init__( self, operands ):
		super().__init__( operands)
	def numOperands(self):
		return self.numChildren()
	def getOperand(self, index):
		return self.children[index]
	def getOperands(self):
		return self.getAllChildren()



class Negate( Operator ):
	def __init__(self, operand):
		super().__init__([operand])
	def deep_copy(self):
		return Negate(self.getOperand(0).deep_copy())
	def toLatex(self):
		return "-" + self.getOperand(0).toLatex()




class Addition( Operator ):
	def __init__(self, operands):
		flattened = []
		# here we flatten nested Multiplications
		for op in operands:
			if op.__class__ == Addition:
				for i in range(op.numOperands()):
					flattened.append(op.getOperand(i))
			else:
				flattened.append(op)
		super().__init__(flattened)
	def deep_copy(self):
		return Addition( [op.deep_copy() for op in self.getOperands()] )
	def toLatex(self):
		result = ""
		for i in range(self.numOperands()):
			expr = self.getOperand(i)
			if isinstance(expr,Negate):
				result+="-" + expr.getOperand(0).toLatex()
			else:
				if i> 0:
					result+="+" + expr.toLatex()
				else:
					result+=expr.toLatex()
		return result




class Multiplication( Operator ):
	def __init__(self, operands):
		flattened = []
		# here we flatten nested Multiplications
		for op in operands:
			if op.__class__ == Multiplication:
				for i in range(op.numOperands()):
					flattened.append(op.getOperand(i))
			else:
				flattened.append(op)

		super().__init__(flattened)
	def deep_copy(self):
		return Multiplication( [op.deep_copy() for op in self.getOperands()] )
	def toLatex(self):
		result = ""
		for i in range(self.numOperands()):
			op = self.getOperand(i)
			parentheses = False

			if op.__class__ == Addition or op.__class__ == Negate:
				parentheses = True

			if parentheses:
				result += "\\left(" + op.toLatex() + "\\right)"
			else:
				result += op.toLatex()
		return result

class Quotient(Operator):
	def __init__( self, numerator, denominator ):
		super().__init__([numerator, denominator])
	def getNumerator(self):
		return self.getOperand(0)
	def getDenominator(self):
		return self.getOperand(1)
	def deep_copy(self):
		return Quotient(self.getNumerator().deep_copy(), self.getDenominator().deep_copy())
	def toLatex(self):
		return "\\frac{{{}}}{{{}}}".format(self.getNumerator().toLatex(),self.getDenominator().toLatex())

class ScopedOperator( Operator ):
	# this is an operator which establishes a scope for a specific variable within the whole expression
	# an example would be integration, where the variable would be the integration variable. all child
	# expressions of the integration, which contain the integration variable, must not leave the scope
	def __init__(self, expr, variable, other_operands):
		pass

class Integration( Operator ):
	def __init__( self, integrand, variable ):
		super().__init__([integrand, variable])
	def getIntegrand(self):
		return self.children[0]
	def getVariable(self):
		return self.children[1]
	def deep_copy(self):
		return Integration(self.getIntegrand().deep_copy(), self.getVariable().deep_copy())
	def toLatex(self):
		integrand_latex = self.getIntegrand().toLatex()

		# special case: if the integrand is 1, we wont render it
		if self.getIntegrand().__class__ == Number:
			if self.getIntegrand().getValue() == 1:
				integrand_latex = ""

		return "\int{" + integrand_latex + "\\mathbf{d}" + self.getVariable().getSymbol() + "}"

class Summation( Operator ):
	def __init__(self, expr, index, start, end):
		self.index = index
		super().__init__([expr, start, end])
	def getIndex(self):
		return self.index
	def getExpr(self):
		return self.getOperand(0)
	def getStart(self):
		return self.getOperand(1)
	def getEnd(self):
		return self.getOperand(2)
	def deep_copy(self):
		return Summation( self.getExpr().deep_copy(), self.getIndex().deep_copy(), self.getStart().deep_copy(), self.getEnd().deep_copy() )
	def toLatex(self):
		result = "\\sum_{{{}={}}}^{{{}}}{{{}}}".format(self.index.toLatex(), self.getStart().toLatex(), self.getEnd().toLatex(), self.getExpr().toLatex())
		return result



def num( value ):
	if value < 0:
		return neg(num(-value))
	return Number(value)


def var( symbol ):
	return Variable(symbol)

def tensor_component( tensor_symbol, component_symbol ):
	return TensorComponent(tensor_symbol, component_symbol)

def fun( symbol, *args ):
	return Function(symbol, list(args))

def add( *args ):
	return Addition(list(args))

def sub( a, b ):
	return add( a, neg(b) )

def mul( *args ):
	return Multiplication(list(args))

def div( numerator, denominator ):
	return Quotient(numerator, denominator)

def pow( value ):
	pass

def neg( expr ):
	return Negate(expr)

def integrate( integrand, variable ):
	return Integration( integrand, variable )

def sum( expr, index, start, end ):
	return Summation(expr, index, start, end)

def infty():
	return var("\\infty")

def kronecker( index_i, index_j ):
	return Kronecker(index_i, index_j)



def latex( expr ):
	return expr.toLatex()

def hierarchy( expr ):
	return expr.toHierarchy()



class FoldConstants(object):
	def visit_Addition(self, expr):
		sum = None
		other_childs = []

		numOperands = expr.numOperands()
		for i in range(numOperands):
			child = expr.getOperand(i)
			if isinstance(child, Number):
				if not sum:
					sum = child.getValue()
				else:
					sum += child.getValue()
				continue
			if isinstance(child, Negate):
				if isinstance( child.getOperand(0), Number ):
					if not sum:
						sum = -child.getOperand(0).getValue()
					else:
						sum += -child.getOperand(0).getValue()
					continue
			other_childs.append(child)


		if sum and other_childs:
			# be aware that sum == False in case when sum is 0
			return Addition( other_childs + [num(sum)] )
		elif sum:
			return num(sum)
		else:
			if len(other_childs) == 1:
				return other_childs[0]
			return Addition( other_childs  )
				



	'''
	def visit_Multiplication(self, expr):
		number_childs = []
		other_childs = []
		numNegates = 0

		numOperands = expr.numOperands()
		for i in range(numOperands):
			child = expr.getOperand(i)
			if isinstance(child, Number):
				number_childs.append(child)
			elif isinstance(child, Negate):
				numNegates += 1
				other_childs.append(child.getOperand(0))
			else:
				other_childs.append(child)

		folded_number = None
		if number_childs:
			folded_number = number_childs[0].getValue()

		for i in range(1, len(number_childs)):
			folded_number *= number_childs[i].getValue()


		folded_result = None

		if not other_childs:
			if not folded_number:
				raise ValueError
			# all childs are numbers, we can return a pure number
			folded_result = num(folded_number)
		elif folded_number:
			# there are some expression childs and a folded number
			folded_result = Multiplication( [folded_number] + other_childs )
		else:
			# there are only other childs, no folded number
			folded_result = Multiplication(other_childs)


		if numNegates % 2 == 0:
			return folded_result
		else:
			return neg(folded_result)
	'''

	def visit_Number(self, expr):
		# turn negative numbers into negates
		# here we actually do the opposite of folding, however
		# using negate is much better for equation manipulation
		# since it allows to explicitly treat the sign
		if expr.getValue() < 0:
			return neg(num(-expr.getValue()))
		return expr
	def visit_Negate(self, expr):
		arg = expr.getOperand(0)
		if isinstance(arg, Negate):
			# double negation
			return arg.getOperand(0)
		return expr

class SHRecursiveRelation(object):
	def visit_Multiplication(self, expr):
		#print("visit_Multiplication {}".format(expr.__class__.__name__))
		children = []
		omega_index = None
		shbasis_index = None
		pairs = []

		# check if multiplication has \omega_x and shbasis
		for i in range(expr.numChildren()):
			child = expr.getChildren(i)
			children.append(child)
			if not omega_index and child.__class__ == TensorComponent:
				if child.getSymbol() == "\\omega":
					omega_index = i
			if not shbasis_index and child.__class__ == SHBasis:
				if child.getArgument(2).__class__ == Variable:
					if child.getArgument(2).getSymbol() == "\\omega":
						shbasis_index = i
			if omega_index != None and shbasis_index != None:
				pairs.append((omega_index, shbasis_index))
				omega_index = None
				shbasis_index = None

		if pairs:
			for pair in pairs:
				omega = children[pair[0]]
				shbasis = children[pair[1]]
				l = shbasis.get_l()
				m = shbasis.get_m()
				sharg = shbasis.get_direction_argument()

				a = fun( "a", sub(l, num(1)), m )
				a.setAllSuperScripts()
				a_basis = SHBasis(sub(l, num(1)), m, sharg, conjugate_complex = True)
				b = fun( "b", add(l, num(1)), m )
				b.setAllSuperScripts()
				b_basis = SHBasis(add(l, num(1)), m, sharg, conjugate_complex = True)
				c = fun( "c", sub(l, num(1)), sub(m, num(1)) )
				c.setAllSuperScripts()
				c_basis = SHBasis(sub(l, num(1)), sub(m, num(1)), sharg, conjugate_complex = True)
				d = fun("d", add(l, num(1)), sub(m, num(1)))
				d.setAllSuperScripts()
				d_basis = SHBasis(add(l, num(1)), sub(m, num(1)), sharg, conjugate_complex = True)
				e = fun("e", sub(l, num(1)), add(m, num(1)))
				e.setAllSuperScripts()
				e_basis = SHBasis(sub(l, num(1)), add(m, num(1)), sharg, conjugate_complex = True)
				f = fun("f", add(l, num(1)), add(m, num(1)))
				f.setAllSuperScripts()
				f_basis = SHBasis(add(l, num(1)), add(m, num(1)), sharg, conjugate_complex = True)

				if omega.__class__ == TensorComponent:
					if omega.component_symbol == 'x':
						# recurrence relation for w_xYlm
						term0 = neg( mul( c, c_basis ) )
						term1 = mul( d, d_basis )
						term2 = mul( e, e_basis )
						term3 = neg( mul( f, f_basis ) )
						children[pair[0]] = mul( div(num(1), num(2)), add( term0, term1, term2, term3 ) )
					elif omega.component_symbol == 'y':
						# recurrence relation for w_yYlm
						term0 = mul( c, c_basis )
						term1 = neg( mul( d, d_basis ) )
						term2 = mul( e, e_basis )
						term3 = neg( mul( f, f_basis ) )
						children[pair[0]] = mul( div(var('i'), num(2)), add( term0, term1, term2, term3 ) )
					elif omega.component_symbol == 'z':
						term0 = mul( a, a_basis )
						term1 = mul( b, b_basis )
						children[pair[0]] = add( term0, term1 )
					del children[pair[1]]

		if len(children) == 1:
			return children[0]
		return Multiplication(children)

class SHOrthogonalityProperty(object):
	def visit_Integration(self, expr):
		if expr.getChildren(0).__class__ == Multiplication:
			m = expr.getChildren(0)
			if m.numOperands() == 2:
				op0 = m.getOperand(0)
				op1 = m.getOperand(1)
				if op0.__class__ == SHBasis and op1.__class__ == SHBasis:
					# both operands are SH basis functions!
					# make sure one is conjugate complex, and the other is not
					if (op0.is_conjugate_complex() and not op1.is_conjugate_complex()) or (op1.is_conjugate_complex() and not op0.is_conjugate_complex()):
						# identity applies
						return mul(kronecker(op0.get_l(), op1.get_l()), kronecker(op0.get_m(), op1.get_m()))
		return expr
	def visit_Multiplication(self, expr):
		# flatten nested multiplications
		factors = []
		for i in range(expr.numOperands()):
			f = expr.getOperand(i)
			if f.__class__ == Multiplication:
				for j in range(f.numOperands()):
					factors.append(f.getOperand(j))
			else:
				factors.append(f)
		return Multiplication(factors)
	def visit_Addition(self, expr):
		# flatten nested additions
		factors = []
		for i in range(expr.numOperands()):
			f = expr.getOperand(i)
			if f.__class__ == Addition:
				for j in range(f.numOperands()):
					factors.append(f.getOperand(j))
			else:
				factors.append(f)
		return Addition(factors)


class DistributiveLaw(object):
	def __init__(self, debug = False):
		self.debug = debug
	def visit_Multiplication(self, expr):
		add = None
		others = []
		factors_before = 0
		factors_after = 0
		# find addition and the rest
		for i in range(expr.numOperands()):
			op = expr.getOperand(i)
			if add == None and op.__class__ == Addition:
				add = op
			else:
				if not add:
					factors_before += 1
				else:
					factors_after += 1
				others.append(op)

		#	print(others)
		if add == None:
			# no addition term...do nothing
			return expr

		terms = []
		for i in range(add.numOperands()):
			others_copies = [ op.deep_copy() for op in others]
			if factors_before >= factors_after:
				terms.append( Multiplication( others_copies + [add.getOperand(i)]) )
			else:
				terms.append( Multiplication([add.getOperand(i)] + others_copies) )
				
		return Addition(terms)
	def visit_Addition(self, expr):
		# flatten nested additions
		factors = []
		for i in range(expr.numOperands()):
			f = expr.getOperand(i)
			if f.__class__ == Addition:
				for j in range(f.numOperands()):
					factors.append(f.getOperand(j))
			else:
				factors.append(f)
		return Addition(factors)

class CleanupSigns(object):
	def visit_Integration(self, expr):
		integrand = expr.getOperand(0)
		variable = expr.getOperand(1)

		if integrand.__class__ == Negate:
			return neg(Integration(integrand.getOperand(0), variable))
		else:
			return expr


	def visit_Multiplication(self, expr):
		numNegates = 0
		factors = []
		# find negations 
		for i in range(expr.numOperands()):
			op = expr.getOperand(i)
			if op.__class__ == Negate:
				numNegates += 1
				factors.append(op.getOperand(0))
			else:
				factors.append(op)
		if numNegates % 2 == 0:
			return Multiplication(factors)
		else:
			return neg(Multiplication(factors))

class SplitIntegrals(object):
	def visit_Integration(self, expr):
		integrand = expr.getOperand(0)
		variable = expr.getOperand(1)
		terms = []
		if integrand.__class__ == Addition:
			for i in range(integrand.numOperands()):
				terms.append( Integration(integrand.getOperand(i), variable) )
			return Addition(terms)
		return expr
		
class Substitute(object):
	def __init__(self, expr, replacement):
		self.expr = expr
		self.replacement = replacement

	def visit_Expression(self, expr):
		if self.expr == expr:
			return self.replacement.deep_copy()
		return expr
'''
	def visit_Multiplication(self, expr):
		# if one of the terms is an integrand or summation
		# we will move all other terms into the integral/summation

		# TODO: now it becomes tricky in terms of domain dependencies
		# moving factors like this can cause invalid manipulations,
		# when the other factors are domain holders as well
		# currently we blindly assume that one can move all factors
		summation_or_integral = None
		other_factors = []
		for i in range(expr.numOperands()):
			factor = expr.getOperand(i)
			if summation_or_integral == None and (factor.__class__ == Summation or factor.__class__ == Integration):
				summation_or_integral = factor.deep_copy()
				continue
			other_factors.append(factor.deep_copy())
		if summation_or_integral != None and len(other_factors) > 0:
			current = summation_or_integral.getOperand(0)
			summation_or_integral.setChildren( 0, Multiplication( other_factors + [current] ) )
			return summation_or_integral
		return expr
'''

class SwitchDomains(object):
	def visit_Integration(self, expr):
		# if we come across an integration which holds a sum,
		# then we can switch the domains easily
		if expr.getOperand(0).__class__ == Summation:
			summation = expr.getOperand(0)
			summation_child = summation.getOperand(0)
			# now we set the child of the integration to the child of the summation
			expr.setChildren(0, summation_child)
			# and we set the integration to the child of the summation
			summation.setChildren(0, expr)
			# the summation is now the head
			return summation
		# if we come across an integration which holds a multiplication of which
		# one factor is a sum, then we will move all factors into the sum
		# and switch domains
		elif expr.getOperand(0).__class__ == Multiplication:
			m = expr.getOperand(0)
			outside_factors = []
			summation = None
			for i in range(m.numOperands()):
				op = m.getOperand(i)
				if not summation and op.__class__ == Summation:
					summation = op
				else:
					outside_factors.append(op)
			if summation:
				summation_child = summation.getOperand(0)
				# now we set the child of the integration to the child of the summation
				# multiplied by the outside factors
				expr.setChildren(0, Multiplication(outside_factors + [summation_child]))
				# and we set the integration to the child of the summation
				summation.setChildren(0, expr)
				# the summation is now the head
				return summation
		return expr

class Factorize(object):
	def visit_Integration(self, expr):
		# if the interal holds a multiplication operator,
		# then we will try to pull out all factors which
		# do not depend on the domain variable
		if expr.getOperand(0).__class__ == Multiplication:
			multiplication = expr.getOperand(0)
			# iterate all factors 
			non_dependent_factors = []
			dependent_factors = []
			for i in range(multiplication.numOperands()):
				op = multiplication.getOperand(i)
				if not op.depends_on( expr.getVariable() ):
					non_dependent_factors.append(op)
				else:
					dependent_factors.append(op)
			# now we extract all non dependent factors and set only the remaining factors
			if len(dependent_factors) ==0:
				expr.setChildren(0, num(1))
			elif len(dependent_factors) ==1:
				expr.setChildren(0, dependent_factors[0])
			else:
				expr.setChildren(0, Multiplication(dependent_factors))

			return Multiplication( non_dependent_factors + [expr] )

		# do nothing
		return expr



class SummationOverKronecker(object):
	def visit_Summation(self, expr):
		var = expr.getIndex()
		#print(expr.getExpr().__class__.__name__)
		if expr.getExpr().__class__ == Multiplication:
			m = expr.getExpr()
			factors = []
			kroneckers = []
			for i in range(m.numOperands()):
				f = m.getOperand(i)
				if f.__class__ == Kronecker:
					if f.getFirstIndex() == var or f.getSecondIndex() == var:
						kroneckers.append(f)
						continue
				factors.append(f)
			# TODO: account for multiple kroneckers, selecting different terms
			# TODO: deep_copy()
			terms = []
			for k in kroneckers:
				term = None
				if len(factors) == 0:
					term = num(1)
				elif len(factors) == 1:
					term = factors[0]
				else:
					term = Multiplication(factors)
				if k.getFirstIndex() == var:
					terms.append( apply_recursive(term, Substitute(var, f.getSecondIndex())) )
				elif k.getSecondIndex() == var:
					terms.append( apply_recursive(term, Substitute(var, f.getFirstIndex())) )
			if len(terms) == 1:
				return terms[0]
			elif len(terms) > 1:
				return Addition(terms)
		return expr




def apply( expr, cls, visitor ):
	#indent_str = expr.indent_string(indent)
	#print("{}{}".format(indent_str, cls.__name__))
	try:
		#print("{}trying {}".format(indent_str, cls.__name__))
		visitor_function = getattr(visitor.__class__, "visit_{}".format(cls.__name__))
		#print("{}yepp!".format(indent_str))
		# there is a visitor function:call it
		return visitor_function(visitor, expr)
	except AttributeError:
		pass
		#print("{}nope".format(indent_str))

	# Iterates the class hierarchy of the Expression class (and its subclasses).
	# Tries to call visitor function for each subclass name (starting from the leave).
	bases = cls.__bases__
	for basecls in bases:
	#for i in range(len(bases)):
		#basecls = bases[i]
		# recursve up within the class hierarchy
		result = apply( expr, basecls, visitor )
		# we are done as soon as this function returns something valid
		if result != None:
			return result

	# our visitor didnt have a visit function for the given expression type
	# or any of its base functions
	return expr



def apply_recursive( expr, visitor ):
	for i in range(expr.numChildren()):
		expr.setChildren(i, apply_recursive(expr.getChildren(i), visitor))
	return apply(expr, expr.__class__, visitor)


if __name__ == "__main__":

	# setup equation

	omega = var("\\omega")
	omega_x = tensor_component("\\omega", "x")
	omega_y = tensor_component("\\omega", "y")
	omega_z = tensor_component("\\omega", "z")
	dx = var("\\partial_x")
	dy = var("\\partial_y")
	dz = var("\\partial_z")
	L = fun( "L", var("\\vec{x}"), omega)
	Ylm = SHBasis(sub(var("l'"), num(1)), var("m'"), omega, conjugate_complex=True)

	a = fun("a", sub(var("l'"), num(1)), var("m'"))
	a.setAllSuperScripts()
	b = fun("b", add(var("l'"), num(1)), var("m'"))
	b.setAllSuperScripts()

	# expression for the expanded radiance field
	L_expanded = sum( sum( mul( SHCoefficient( "L", var("l"), var("m"), var("\\vec{x}") ), SHBasis(var("l"), var("m"), omega, conjugate_complex=False) ), var('m'), neg(var('l')), var('l') ), var('l'), num(0), infty() )

	#expr = mul( a, var("\\partial_x"), integrate( mul(omega_x, Ylm, L), omega ) )
	expr_a = add( mul(a, SHBasis(sub(var("l'"), num(1)), var("m'"), omega, conjugate_complex=True)), mul(b, SHBasis(add(var("l'"), num(1)), var("m'"), omega, conjugate_complex=True)) )
	expr_b = add( mul(omega_x, dx, L), mul(omega_y, dy, L), mul(omega_z, dz, L))
	expr = integrate( mul(expr_a, expr_b), omega )

	#'''
	#print("$$\n" + latex(expr) + "\n$$")
	expr = apply_recursive(expr, DistributiveLaw())
	expr = apply_recursive(expr, DistributiveLaw())
	expr = apply_recursive(expr, SplitIntegrals())
	expr = apply_recursive(expr, Factorize())
	expr = apply_recursive(expr, SHRecursiveRelation())
	expr = apply_recursive(expr, FoldConstants())
	expr = apply_recursive(expr, Factorize())
	expr = apply_recursive(expr, DistributiveLaw())
	expr = apply_recursive(expr, SplitIntegrals())
	expr = apply_recursive(expr, CleanupSigns())
	expr = apply_recursive(expr, Factorize())
	expr = apply_recursive(expr, DistributiveLaw(True))
	expr = apply_recursive(expr, CleanupSigns())
	expr = apply_recursive(expr, Substitute(L, L_expanded))
	expr = apply_recursive(expr, SwitchDomains())
	expr = apply_recursive(expr, SwitchDomains())
	expr = apply_recursive(expr, Factorize())
	expr = apply_recursive(expr, SHOrthogonalityProperty())
	#expr = expr.getOperand(0)
	#print(hierarchy(expr))
	print("\n----------------------------\n")
	print("$$\n" + latex(expr) + "\n$$")
	expr = apply_recursive(expr, SummationOverKronecker())
	print("\n----------------------------\n")
	print("$$\n" + latex(expr) + "\n$$")


	#expr = add( expr.getOperand(0), expr.getOperand(1) )
	#'''



	'''
	expr = integrate( mul(var('a'), sum(mul(var('b'), var('c')), var('i'), num(0), infty())), var('\\omega') )
	print("\n----------------------------\n")
	print("$$\n" + latex(expr) + "\n$$")
	expr = apply_recursive(expr, SwitchDomains())

	print("\n----------------------------\n")
	print("$$\n" + latex(expr) + "\n$$")
	'''


	#expr = apply_recursive(expr, MoveFactors())
	#print("\n----------------------------\n")
	#print("$$\n" + latex(expr) + "\n$$")

	'''
	# run manipulations
	print("$$\n" + latex(expr) + "\n$$")
	expr = apply_recursive(expr, SplitIntegrals())
	expr = apply_recursive(expr, CleanupSigns())
	
	expr = apply_recursive(expr, MoveFactors())
	expr = apply_recursive(expr, MoveFactors())
	expr = apply_recursive(expr, MoveFactors())
	expr = apply_recursive(expr, FoldConstants())
	expr = apply_recursive(expr, DistributiveLaw())
	expr = apply_recursive(expr, CleanupSigns())
	print("$$\n" + latex(expr) + "\n$$")
	'''




	#print(hierarchy(expr))

	#print("$$\n" + latex(expr) + "\n$$")
	#print("$$\n" + latex(expr) + "\n$$")
	'''
	#

	#expr = kronecker(var('i'), var('j'))
	#print("$$\n" + latex(expr) + "\n$$")

	# output result
	#print(hierarchy(expr))
	#print(latex(expr))
	'''

