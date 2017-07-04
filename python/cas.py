

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

	def setChildren(self, index, new_child):
		self.children[index] = new_child

	def fold_constants(self):
		return self


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

class Variable(Expression):
	def __init__( self, symbol ):
		super().__init__()
		self.symbol = symbol

	def getName(self):
		return self.symbol
	def toLatex(self):
		return str(self.symbol)
	def toHierarchy( self, indent = 0 ):
		return self.indent_string(indent) + "Variable {}\n".format(self.symbol)
		
class Function( Expression ):
	def __init__( self, symbol, arguments ):
		super().__init__(arguments)
		self.symbol = symbol

	def getArgument(self, index):
		return self.getChildren(index)

	def toLatex(self):
		string = str(self.symbol)
		if not self.children:
			return string
		else:
			numArgs = self.numChildren()
			string += "("
			for i in range(numArgs):
				string += self.children[i].toLatex()
				if i < numArgs-1:
					string += ", "
			string += ")"
			return string


class SHBasis(Function):
	def __init__(self, l, m, omega):
		super().__init__( "Y", [l, m, omega])

	def get_l(self):
		return self.getChildren(0)
	def get_m(self):
		return self.getChildren(1)
	def get_direction_argument(self):
		return self.getChildren(2)

	def toLatex(self):
		l = self.children[0]
		m = self.children[1]
		omega = self.children[2]
		return "\overline{Y^{" + l.toLatex() + ", " + m.toLatex() + "}}(" + omega.toLatex() + ")"
	def fold_constants(self):
		return SHBasis( self.l.fold_constants(), self.m.fold_constants(), self.omega.fold_constants() )


class Operator( Expression ):
	def __init__( self, operands ):
		super().__init__( operands)
	def numOperands(self):
		return self.numChildren()
	def getOperand(self, index):
		return self.children[index]



class Negate( Operator ):
	def __init__(self, operand):
		super().__init__([operand])
	def toLatex(self):
		return "-" + self.getOperand(0).toLatex()
	def fold_constants(self):
		arg = self.getOperand(0)
		if isinstance(arg, Number):
			# turn negative numbers into negates
			# here we actually do the opposite of folding, however
			# using negate is much better for equation manipulation
			# since it allows to explicitly treat the sign
			if arg.getValue() < 0:
				return neg(-arg.getValue())
		if isinstance(arg, Negate):
			# double negation
			return arg.getOperand(0)
		return self



class Addition( Operator ):
	def __init__(self, operands):
		super().__init__(operands)
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
	def fold_constants(self):
		return Addition( [operand.fold_constants() for operand in self.operands] )



class Multiplication( Operator ):
	def __init__(self, operands):
		super().__init__(operands)
	def toLatex(self):
		result = ""
		for i in range(self.numOperands()):
			op = self.getOperand(i)
			parentheses = False

			if op.__class__ == Addition:
				parentheses = True

			if parentheses:
				result += "\\left(" + op.toLatex() + "\\right)"
			else:
				result += op.toLatex()
		return result
	def fold_constants(self):
		number_childs = []
		other_childs = []
		numNegates = 0

		numOperands = self.numOperands()
		for i in range(numOperands):
			child = self.getOperand(i).fold_constants()
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

		if folded_number:
			if folded_number == -1:
				numNegates += 1
				folded_number = None


		folded_result = None

		if not other_childs:
			# all childs are numbers, we can return a pure number
			folded_result = num(folded_number)
		elif folded_number:
			# there are some expression childs and a folded number
			folded_result = Multiplication( [folded_number] + other_childs )
		else:
			folded_result = Multiplication(other_childs)


		if numNegates % 2 == 0:
			return folded_result
		else:
			return neg(folded_result)


class Integration( Operator ):
	def __init__( self, integrand, variable ):
		super().__init__([integrand, variable])
	def getIntegrand(self):
		return self.children[0]
	def getVariable(self):
		return self.children[1]
	def toLatex(self):
		return "\int{" + self.getIntegrand().toLatex() + "\\mathbf{d}" + self.getVariable().getName() + "}"


def num( value ):
	return Number(value)


def var( symbol ):
	return Variable(symbol)

def fun( symbol, *args ):
	return Function(symbol, list(args))

def add( *args ):
	return Addition(list(args))

def sub( a, b ):
	return add( a, neg(b) )

def sum( value ):
	pass

def mul( *args ):
	return Multiplication(list(args))

def pow( value ):
	pass

def neg( expr ):
	return Negate(expr)

def integrate( integrand, variable ):
	return Integration( integrand, variable )



def latex( expr ):
	return expr.toLatex()

def hierarchy( expr ):
	return expr.toHierarchy()


#def fold_constants( expr ):
#	return expr.fold_constants()




class FoldConstants(object):
	pass

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
			if not omega_index and child.__class__ == Variable:
				if child.getName() == "\omega_x":
					omega_index = i
			if not shbasis_index and child.__class__ == SHBasis:
				if child.getArgument(2).__class__ == Variable:
					if child.getArgument(2).getName() == "\omega":
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

				# recurrence relation for w_xYlm
				term0 = neg( SHBasis(sub(l, num(1)), sub(m, num(1)), sharg) )
				term1 = SHBasis(add(l, num(1)), sub(m, num(1)), sharg)
				term2 = SHBasis(sub(l, num(1)), add(m, num(1)), sharg)
				term3 = neg( SHBasis(add(l, num(1)), add(m, num(1)), sharg) )
				children[pair[0]] = add( term0, term1, term2, term3 )
				del children[pair[1]]

		if len(children) == 1:
			return children[0]
		return Multiplication(children)


class DistributiveLaw(object):
	def visit_Multiplication(self, expr):
		add = None
		others = []
		# find addition and the rest
		for i in range(expr.numOperands()):
			op = expr.getOperand(i)
			if add == None and op.__class__ == Addition:
				add = op
			else:
				others.append(op)

		if add == None:
			# no addition term...do nothing
			return expr

		terms = []
		for i in range(add.numOperands()):
			terms.append( Multiplication([add.getOperand(i)] + others) )
		return Addition(terms)


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
		


def apply( expr, cls, visitor, indent = 0 ):
	indent_str = expr.indent_string(indent)
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
		result = apply( expr, basecls, visitor, indent+1 )
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
	omega = var("\omega")
	wx = var("\omega_x")
	L = fun( "L", var("\\vec{x}"), omega)
	Ylm = SHBasis(sub(var("l'"), num(1)), var("m'"), omega)
	#expr = mul(Ylm, L)
	expr = integrate( mul(wx, Ylm, L), omega )

	#print(latex(expr))
	#print(hierarchy(expr))

	# run manipulations
	#expr = fold_constants( expr )
	#expr = apply( expr, FoldConstants() )
	expr = apply_recursive(expr, SHRecursiveRelation())
	expr = apply_recursive(expr, DistributiveLaw())
	expr = apply_recursive(expr, SplitIntegrals())
	print(hierarchy(expr))
	expr = apply_recursive(expr, CleanupSigns())
	#expr = fold_constants( expr )



	# output result
	#print(hierarchy(expr))
	print(latex(expr))

