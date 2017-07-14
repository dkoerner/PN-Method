import traceback
import math
import shtools

# only needed for exporting
import numpy as np
import scipy.io

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

	def depends_on(self, variable, debug = False):
		for i in range(self.numChildren()):
			if self.getChildren(i).depends_on(variable, debug):
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
	def __init__( self, symbol, children = [] ):
		super().__init__(children)
		self.symbol = symbol

	def getSymbol(self):
		return self.symbol

	def depends_on(self, variable, debug = False):
		for i in range(self.numChildren()):
			if self.getChildren(i).depends_on(variable, debug):
				return True

		result = self.symbol == variable.symbol
		#if debug:
		#	print("Variable depends_on: self.symbol={} variable.symbol={} result={}".format(self.symbol, variable.symbol, result))
		return result
		

	def __eq__(self, other):
		if other.__class__ == Variable:
			if other.getSymbol() == self.getSymbol():
				return True
		return False

	def deep_copy(self):
		return Variable(self.getSymbol())

	def toLatex(self):
		return str(self.getSymbol())

	def toHierarchy( self, indent = 0 ):
		return self.indent_string(indent) + "Variable {}\n".format(self.symbol)

class ImaginaryUnit(Variable):
	def __init__(self):
		super().__init__("i")
	def depends_on(self, variable, debug = False):
		return False
	def deep_copy(self):
		return ImaginaryUnit()

class Tensor(Variable):
	def __init__( self, symbol, rank, dimension ):
		self.dimension = dimension
		self.rank = rank
		self.has_components = False
		self.collapsed = True # causes just to print the symbol in latex

		dimension_symbols=['x','y','z']
		self.component_symbol_to_index = {}
		for s in range(len(dimension_symbols)):
			self.component_symbol_to_index[dimension_symbols[s]] = s

		numComponents = dimension**rank
		components = []
		indices = [0 for i in range(rank)]
		for c in range(numComponents):
			components.append(TensorComponent(symbol, "".join([ dimension_symbols[k] for k in indices ])))

			for i in range(rank):
				indices[i] += 1
				if indices[i] < dimension:
					break
				indices[i] = 0

		super().__init__(symbol, components)

	def getChildIndex(self, indices):
		child_index = 0
		for i in range(self.rank):
			if type(indices) is list:
				if isinstance(indices[i], str):
					index = self.component_symbol_to_index[indices[i]]
				else:
					index = indices[i]
			else:
				if isinstance(indices, str):
					index = self.component_symbol_to_index[indices]
				else:
					index = indices
			child_index += index*(self.dimension**i)
		return child_index

	def setComponent( self, indices, expr ):
		 # now since we have some content, this tensor is not collapsed during latex rendering
		self.collapsed = False
		return self.setChildren(self.getChildIndex(indices), expr)

	def getComponent(self, indices = []):
		return self.getChildren(self.getChildIndex(indices))

	def numComponents(self):
		return self.numChildren()


	def deep_copy(self):
		cpy = Tensor(self.getSymbol(), self.rank, self.dimension)
		for i in range(self.numChildren()):
			cpy.setChildren(i, self.getChildren(i))
		cpy.collapsed = self.collapsed
		return cpy

	def toLatex(self):
		if self.collapsed:
			return self.getSymbol() + " "

		if self.rank == 0:
			return self.getComponent().toLatex()
		elif self.rank == 1:
			result = ""
			result += "\\left["
			result += "\\begin{array}"
			result += "\\ " # this is due to a bug with latex rendering in jupyer
			for i in range(self.dimension):
				result += self.getComponent([i]).toLatex() + "\\\\"
			result += "\\end{array}"
			result += "\\right]"		
			return result
		elif self.rank == 2:
			result = ""
			result += "\\left["
			result += "\\begin{array}"
			result += "\\ " # this is due to a bug with latex rendering in jupyer
			for i in range(self.dimension):
				for j in range(self.dimension):
					result += self.getComponent([i, j]).toLatex()
					if j < self.dimension-1:
						result += " & "
				#result + "\\\\"
				result += "\\\\"
			result += "\\end{array}"
			result += "\\right]"		
			return result
		return self.getSymbol()
	def toHierarchy( self, indent = 0 ):
		return self.indent_string(indent) + "Tensor {}\n".format(self.getSymbol())



class TensorComponent(Variable):
	def __init__( self, tensor_symbol, component_symbol ):
		super().__init__(tensor_symbol)
		self.component_symbol = component_symbol
	def deep_copy(self):
		return TensorComponent(self.getSymbol(), self.component_symbol)
	def toLatex(self):
		return "{}_{{{}}}".format(self.getSymbol(), self.component_symbol)
	def toHierarchy( self, indent = 0 ):
		return self.indent_string(indent) + "TensorComponent {} {}\n".format(self.getSymbol(), self.component_symbol)


		
class Function( Expression ):
	def __init__( self, symbol, arguments, arg_symbols = None, body_expr = None ):
		# arg_symbols     : is list of variable symbols, which relate the function arguments
		# 					to variables in the expression for the function body
		#					this is only reqiuired, when this function should have a body for eval
		# body_expr       : an expression which expresses what the function is doing
		#					this is only reqiuired, when this function should be able to evaluate
		super().__init__(arguments)
		self.symbol = symbol
		self.arg_symbols = arg_symbols # 
		self.body_expr = body_expr # the function body expressed as expression
		self.latexArgumentPositions = [0 for arg in arguments]

		for arg in arguments:
			if isinstance(arg, str):
				print("===================================")
				raise ValueError()

	def getSymbol(self):
		return self.symbol

	def getFunctionBody(self):
		return self.body_expr

	def setFunctionBody( self, body_expr, *arg_symbols ):
		if len(arg_symbols) != self.numChildren():
			raise ValueError("number of argument symbols must match number of argument expressions")
		self.arg_symbols = arg_symbols
		self.body_expr = body_expr # the function body expressed as expression

	def numArguments(self):
		return self.numChildren()

	def getArgument(self, index):
		return self.getChildren(index)

	def getArgumentSymbol(self, index):
		return self.arg_symbols[index]

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

	def setAllSubScripts(self):
		self.latexArgumentPositions = [-1 for arg in self.getArguments()]

	def deep_copy(self):
		body_expr_cp = None
		if self.body_expr:
			body_expr_cp = self.body_expr.deep_copy()
		cp = Function(self.getSymbol(), [arg.deep_copy() for arg in self.getArguments()], self.arg_symbols, body_expr_cp)
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
	def toHierarchy( self, indent = 0 ):
		result = self.indent_string(indent) + "{} {}\n".format(self.__class__.__name__, self.getSymbol())
		for child in self.children:
			result += child.toHierarchy(indent+1)
		return result

	def __eq__(self, other):
		#print(other.__class__.__name__)
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
	def toHierarchy( self, indent = 0 ):
		return self.indent_string(indent) + "SHCoefficient {}\n".format(self.getSymbol())


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
	def getExpr(self):
		return self.getOperand(0)
	def setExpr(self, expr):
		self.setOperand(0, expr)
	def deep_copy(self):
		return Negate(self.getOperand(0).deep_copy())
	def toLatex(self):
		parentheses = False

		if self.getOperand(0).__class__ == Addition:
			parentheses = True

		if not parentheses:
			return "-" + self.getOperand(0).toLatex()
		else:
			return "-\\left(" + self.getOperand(0).toLatex() + "\\right)"


class Sqrt(Operator):
	def __init__(self, operand):
		super().__init__([operand])

	def deep_copy(self):
		return Sqrt(self.getOperand(0).deep_copy())
	def toLatex(self):
		return "\\sqrt{{{}}}".format(self.getOperand(0).toLatex())

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

		if len(flattened)==0:
			raise ValueError()

		if len(flattened)==1:
			raise ValueError("Multiplication requires more than one item")

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

			if op.__class__ == Number and op.getValue() < 0:
				parentheses = True

			if op.__class__ == Summation:
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

class Derivation( Operator ):
	def __init__( self, expr, variable, is_partial = False):
		super().__init__([expr, variable])
		self.is_partial = is_partial
	def getExpr(self):
		return self.getChildren(0)
	def setExpr(self, expr):
		return self.setChildren(0, expr)
	def getVariable(self):
		return self.getChildren(1)
	def toLatex(self):
		result = ""
		if self.is_partial:
			result += "\\partial_{}".format(self.getVariable().toLatex())

		if self.getExpr() != None:
			expr = self.getExpr()
			if isinstance(expr, Variable) or isinstance(expr, Function) or expr.__class__ == Derivation:
				result += "{}".format(expr.toLatex())
			else:
				result += "\\left({}\\right)".format(expr.toLatex())
		else:
			result += "\\left(\\bullet\\right)"
		return result
	def deep_copy(self):
		return Derivation(self.getExpr().deep_copy(), self.getVariable().deep_copy(), self.is_partial)


class Summation( Operator ):
	def __init__(self, expr, index, start, end):
		self.index = index
		super().__init__([expr, start, end])
	def getVariable(self):
		return self.index
	def getExpr(self):
		return self.getOperand(0)
	def setExpr(self, expr):
		self.setChildren(0, expr)
	def getStart(self):
		return self.getOperand(1)
	def getEnd(self):
		return self.getOperand(2)
	def deep_copy(self):
		return Summation( self.getExpr().deep_copy(), self.getVariable().deep_copy(), self.getStart().deep_copy(), self.getEnd().deep_copy() )
	def toLatex(self, parent=None):

		parentheses = False

		if parent != None and parent.__class__ == Multiplication:
			parentheses = True

		if parentheses == True:
			result = "\\sum_{{{}={}}}^{{{}}}{{\\left({}\\right)}}".format(self.index.toLatex(), self.getStart().toLatex(), self.getEnd().toLatex(), self.getExpr().toLatex())
		else:
			result = "\\sum_{{{}={}}}^{{{}}}{{{}}}".format(self.index.toLatex(), self.getStart().toLatex(), self.getEnd().toLatex(), self.getExpr().toLatex())
		 
		return result

class DotProduct(Operator):
	def __init__(self, left, right):
		super().__init__([left, right])
	def getLeft(self):
		return self.getOperand(0)
	def getRight(self):
		return self.getOperand(1)
	def deep_copy(self):
		return DotProduct(self.getLeft().deep_copy(), self.getRight().deep_copy())
	def toLatex(self):
		return "{}\\cdot{}".format(self.getLeft().toLatex(), self.getRight().toLatex())

def num( value ):
	if value < 0:
		return neg(num(-value))
	return Number(value)

# imaginary number
def imag( value ):
	if value == 1:
		return ImaginaryUnit()
	return mul(num(value), ImaginaryUnit())


def var( symbol ):
	return Variable(symbol)

def tensor( symbol, rank, dimension ):
	return Tensor(symbol, rank, dimension)

def tensor_component( tensor_symbol, component_symbol ):
	return TensorComponent(tensor_symbol, component_symbol)

def fun( symbol, *args, **kwargs ):
	f = Function(symbol, list(args))

	if 'arglevel' in kwargs:
		if kwargs['arglevel'] == -1:
			f.setAllSubScripts()
	return f

def add( *args ):
	return Addition(list(args))

def sub( a, b ):
	return add( a, neg(b) )

def mul( *args ):
	return Multiplication(list(args))

def div( numerator, denominator ):
	return Quotient(numerator, denominator)

def sqrt( expr ):
	return Sqrt(expr)

def pow( value ):
	pass

def neg( expr ):
	return Negate(expr)

def integrate( integrand, variable ):
	return Integration( integrand, variable )

def deriv( expr, variable, is_partial = False ):
	return Derivation(expr, variable, is_partial)

def sum( expr, index, start, end ):
	return Summation(expr, index, start, end)

def infty():
	return var("\\infty")

def kronecker( index_i, index_j ):
	return Kronecker(index_i, index_j)

def nabla():
	n = Tensor( "\\nabla", rank=1, dimension=3 )
	n.setComponent(0, var("\\partial_x"))
	n.setComponent(1, var("\\partial_y"))
	n.setComponent(2, var("\\partial_z"))
	n.collapsed = True
	return n
'''
def nabla_new( expr = None ):
	n = Tensor( "\\nabla", rank=1, dimension=3 )
	n.setComponent(0, deriv(expr, var("x"), is_partial = True))
	n.setComponent(1, deriv(expr, var("y"), is_partial = True))
	n.setComponent(2, deriv(expr, var("z"), is_partial = True))
	#n.collapsed = True
	return n
'''
def dot(left, right):
	return DotProduct(left, right)

def sh_expansion( fun, positional_variable, directional_variable ):
	return sum( sum( mul( SHCoefficient( fun.getSymbol(), var("l"), var("m"), positional_variable ), SHBasis(var("l"), var("m"), directional_variable, conjugate_complex=False) ), var('m'), neg(var('l')), var('l') ), var('l'), num(0), infty() )
 
def latex( expr ):
	return expr.toLatex()

def hierarchy( expr ):
	return expr.toHierarchy()

def eval(expr, variables, functions):
	result = apply_recursive(expr.deep_copy(), Evaluate(variables, functions))
	if result != None:
		return result
	return None



class FoldConstants(object):
	def __init__(self, debug = False):
		self.debug = debug

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

		if other_childs and folded_number != None:
			if folded_number == 0:
				return num(0)
			elif folded_number == 1:
				if len(other_childs)==1:
					folded_result = other_childs[0]
				else:
					folded_result = Multiplication( other_childs )
			else:
				# there are some expression childs and a folded number
				folded_result = Multiplication( [num(folded_number)] + other_childs )
		elif other_childs:
			# there are only other childs, no folded number
			if len(other_childs) == 1:
				folded_result = other_childs
			else:
				folded_result = Multiplication( other_childs )
		elif folded_number:
			# all childs are numbers, we can return a pure number
			folded_result = num(folded_number)


		if numNegates % 2 == 0:
			return folded_result
		else:
			return neg(folded_result)


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

class ImaginaryUnitProperty(object):
	def visit_Multiplication(self, expr):
		# we count the number of imaginary units
		numImaginaryUnits = 0
		numOperands = expr.numOperands()
		others = []
		for i in range(numOperands):
			child = expr.getOperand(i)
			if child.__class__ == ImaginaryUnit:
				numImaginaryUnits += 1
				continue
			others.append(child)

		if numImaginaryUnits < 2:
			return expr

		if numImaginaryUnits % 2 == 0:
			if len(others) == 0:
				return neg(num(1))
			elif len(others) == 1:
				return neg(others[0])
			else:
				return neg(Multiplication(others))
		else:
			if len(others) == 0:
				return imag(1)
			else:
				return Multiplication([imag(1)] + others)





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
				if isinstance(child.getArgument(2), Variable):
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

				# a
				expr0 = add(var('l'), neg(var('m')), num(1))
				expr1 = add(var('l'), var('m'), num(1))
				expr2 = add(mul(num(2), var('l')), num(3))
				expr3 = add(mul(num(2), var('l')), num(1))
				a_body = sqrt(div(mul(expr0, expr1), mul(expr2, expr3)))

				# b
				expr0 = add(var('l'), neg(var('m')))
				expr1 = add(var('l'), var('m'))
				expr2 = add(mul(num(2), var('l')), num(1))
				expr3 = add(mul(num(2), var('l')), neg(num(1)))
				b_body = sqrt(div(mul(expr0, expr1), mul(expr2, expr3)))

				# c
				expr0 = add(var('l'), var('m'), num(1))
				expr1 = add(var('l'), var('m'), num(2))
				expr2 = add(mul(num(2), var('l')), num(3))
				expr3 = add(mul(num(2), var('l')), num(1))
				c_body = sqrt(div(mul(expr0, expr1), mul(expr2, expr3)))

				# d
				expr0 = add(var('l'), neg(var('m')))
				expr1 = add(var('l'), neg(var('m')), neg(num(1)))
				expr2 = add(mul(num(2), var('l')), num(1))
				expr3 = add(mul(num(2), var('l')), neg(num(1)))
				d_body = sqrt(div(mul(expr0, expr1), mul(expr2, expr3)))

				# e
				expr0 = add(var('l'), neg(var('m')), num(1))
				expr1 = add(var('l'), neg(var('m')), num(2))
				expr2 = add(mul(num(2), var('l')), num(3))
				expr3 = add(mul(num(2), var('l')), num(1))
				e_body = sqrt(div(mul(expr0, expr1), mul(expr2, expr3)))

				# f
				expr0 = add(var('l'), var('m'))
				expr1 = add(var('l'), var('m'), neg(num(1)))
				expr2 = add(mul(num(2), var('l')), num(1))
				expr3 = add(mul(num(2), var('l')), neg(num(1)))
				f_body = sqrt(div(mul(expr0, expr1), mul(expr2, expr3)))

				a = fun( "a", sub(l, num(1)), m )
				a.setAllSuperScripts()
				a.setFunctionBody( a_body, 'l', 'm' )
				a_basis = SHBasis(sub(l, num(1)), m, sharg, conjugate_complex = True)
				b = fun( "b", add(l, num(1)), m )
				b.setAllSuperScripts()
				b.setFunctionBody( b_body, 'l', 'm' )
				b_basis = SHBasis(add(l, num(1)), m, sharg, conjugate_complex = True)
				c = fun( "c", sub(l, num(1)), sub(m, num(1)) )
				c.setAllSuperScripts()
				c.setFunctionBody( c_body, 'l', 'm' )
				c_basis = SHBasis(sub(l, num(1)), sub(m, num(1)), sharg, conjugate_complex = True)
				d = fun("d", add(l, num(1)), sub(m, num(1)))
				d.setAllSuperScripts()
				d.setFunctionBody( d_body, 'l', 'm' )
				d_basis = SHBasis(add(l, num(1)), sub(m, num(1)), sharg, conjugate_complex = True)
				e = fun("e", sub(l, num(1)), add(m, num(1)))
				e.setAllSuperScripts()
				e.setFunctionBody( e_body, 'l', 'm' )
				e_basis = SHBasis(sub(l, num(1)), add(m, num(1)), sharg, conjugate_complex = True)
				f = fun("f", add(l, num(1)), add(m, num(1)))
				f.setAllSuperScripts()
				f.setFunctionBody( f_body, 'l', 'm' )
				f_basis = SHBasis(add(l, num(1)), add(m, num(1)), sharg, conjugate_complex = True)

				# these functions are the coefficients for the recursive relation of the sh basis function
				# (see p. 4 in the starmap paper)

				'''
				def a_lm( l, m ):
				    return np.sqrt((l-m+1)*(l+m+1)/((2*l+3)*(2*l+1)))
				def b_lm( l, m ):
				    return np.sqrt((l-m)*(l+m)/((2*l+1)*(2*l-1)))
				def c_lm( l, m ):
				    return np.sqrt((l+m+1)*(l+m+2)/((2*l+3)*(2*l+1)))
				def d_lm( l, m ):
				    return np.sqrt((l-m)*(l-m-1)/((2*l+1)*(2*l-1)))
				def e_lm( l, m ):
				    return np.sqrt((l-m+1)*(l-m+2)/((2*l+3)*(2*l+1)))
				def f_lm( l, m ):
				    return np.sqrt((l+m)*(l+m-1)/((2*l+1)*(2*l-1)))
				'''
				a.body2 = lambda l,m: np.sqrt((l-m+1)*(l+m+1)/((2*l+3)*(2*l+1)))
				b.body2 = lambda l,m: np.sqrt((l-m)*(l+m)/((2*l+1)*(2*l-1)))
				c.body2 = lambda l,m: np.sqrt((l+m+1)*(l+m+2)/((2*l+3)*(2*l+1)))
				d.body2 = lambda l,m: np.sqrt((l-m)*(l-m-1)/((2*l+1)*(2*l-1)))
				e.body2 = lambda l,m: np.sqrt((l-m+1)*(l-m+2)/((2*l+3)*(2*l+1)))
				f.body2 = lambda l,m: np.sqrt((l+m)*(l+m-1)/((2*l+1)*(2*l-1)))



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
						children[pair[0]] = mul( div(imag(1), num(2)), add( term0, term1, term2, term3 ) )
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
		add_or_sum = None
		others = []
		factors_before = 0
		factors_after = 0
		# find addition and the rest
		for i in range(expr.numOperands()):
			op = expr.getOperand(i)
			if add_or_sum == None and (op.__class__ == Addition or op.__class__ == Summation):
				add_or_sum = op
			else:
				if not add_or_sum:
					factors_before += 1
				else:
					factors_after += 1
				others.append(op)

		#	print(others)
		if add_or_sum == None:
			# got nothing...return
			return expr

		if add_or_sum.__class__ == Addition:
			# we got an addition term
			add = add_or_sum
			terms = []
			for i in range(add.numOperands()):
				others_copies = [ op.deep_copy() for op in others]
				if factors_before >= factors_after:
					terms.append( Multiplication( others_copies + [add.getOperand(i)]) )
				else:
					terms.append( Multiplication([add.getOperand(i)] + others_copies) )
					
			return Addition(terms)
		elif add_or_sum.__class__ == Summation:
			sum = add_or_sum
			# we got an sum term
			# simply put all the other factors inside
			current = sum.getExpr()
			if factors_before >= factors_after:
				sum.setExpr(Multiplication( others+[current] ))
			else:
				sum.setExpr(Multiplication( [current]+others ))
			return sum
		else:
			raise ValueError("unexpected classtype")


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

	def visit_Quotient(self, expr):
		num = expr.getNumerator()
		denom = expr.getDenominator()
		if num.__class__ == Negate and denom.__class__ == Negate:
			return div(num.getOperand(0), denom.getOperand(0))
		elif num.__class__ == Negate:
			return neg(div(num.getOperand(0), denom))
		elif denom.__class__ == Negate:
			return neg(div(num, denom.getOperand(0)))
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
	def visit_Addition(self, expr):
		terms = []
			# find negations 
		for i in range(expr.numOperands()):
			op = expr.getOperand(i)
			if op.__class__ == Negate:
				if op.getOperand(0).__class__ == Addition:
					nested_addition = op.getOperand(0)
					# we have a negate term, which holds an addition
					## flatten the addition by negating each term of
					# the nested addition
					for j in range(nested_addition.numOperands()):
						terms.append( neg(nested_addition.getOperand(j)) )
					continue
			terms.append(op)
		return Addition(terms)
	def visit_Negate(self, expr):
		arg = expr.getOperand(0)
		if isinstance(arg, Negate):
			# double negation
			return arg.getOperand(0)
		elif arg.__class__ == Addition:
			# split negated addition by negating all its terms
			terms = []
			for i in range(arg.numOperands()):
				op = arg.getOperand(i)
				if op.__class__ == Negate:
					# double negation
					terms.append(op.getOperand(0))
				else:
					terms.append(neg(op))
			return Addition(terms)


		return expr
	def visit_DotProduct(self, expr):
		if expr.getLeft().__class__ == Negate and expr.getLeft().__class__ == Negate:
			# double negation
			return DotProduct(expr.getLeft().getOperand(0), expr.getRight().getOperand(0))
		elif expr.getLeft().__class__ == Negate:
			return neg(DotProduct(expr.getLeft().getOperand(0), expr.getRight()))
		elif expr.getRight().__class__ == Negate:
			return neg(DotProduct(expr.getLeft(), expr.getRight().getOperand(0)))
		else:
			return expr
	def visit_Summation(self, expr):
		if expr.getExpr().__class__ == Negate:
			# summation holds a negate
			neg_expr = expr.getExpr()
			# set the child of the negate to the child of the summation
			expr.setExpr( neg_expr.getExpr() )
			# and wrap the negate around the summation
			return neg(expr)
		return expr
	def visit_Derivation(self, expr):
		if expr.getExpr().__class__ == Negate:
			neg_expr = expr.getExpr()
			expr.setExpr(neg_expr.getExpr())
			return neg(expr)
		return expr



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
		
class SplitDerivatives(object):
	def visit_Derivation(self, expr):
		nested_expr = expr.getExpr()
		variable = expr.getVariable()
		terms = []

		if nested_expr.__class__ == Addition:
			for i in range(nested_expr.numOperands()):
				terms.append( Derivation(nested_expr.getOperand(i), variable, expr.is_partial) )
			return Addition(terms)
		return expr

class SplitSums(object):
	def visit_Summation(self, expr):
		nested_expr = expr.getOperand(0)
		variable = expr.getVariable()
		terms = []
		if nested_expr.__class__ == Addition:
			for i in range(nested_expr.numOperands()):
				terms.append( Summation(nested_expr.getOperand(i), variable, expr.getStart().deep_copy(), expr.getEnd().deep_copy()) )
			return Addition(terms)
		return expr




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
			nested_domain = None
			for i in range(m.numOperands()):
				op = m.getOperand(i)
				if not nested_domain and (op.__class__ == Summation or op.__class__ == Derivation):
					nested_domain = op
				else:
					outside_factors.append(op)
			if nested_domain:
				nested_domain_child = nested_domain.getOperand(0)
				# now we set the child of the parent domain as child of the nested domain
				# multiplied by the outside factors
				expr.setChildren(0, Multiplication(outside_factors + [nested_domain_child]))
				# and we set the integration to the child of the summation
				nested_domain.setChildren(0, expr)
				# the nested domain is now the head
				return nested_domain
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

			if not non_dependent_factors:
				# nothing can be extracted...keep everything as is
				return expr

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
	def visit_Derivation(self, expr):
		# if the derivation holds a multiplication operator,
		# then we will try to pull out all factors which
		# do not depend on the domain variable
		if expr.getOperand(0).__class__ == Multiplication:
			#print("!!!!!!!!!!!!!!!!!!!")
			multiplication = expr.getOperand(0)
			# iterate all factors 
			non_dependent_factors = []
			dependent_factors = []
			for i in range(multiplication.numOperands()):
				op = multiplication.getOperand(i)
				if not op.depends_on( expr.getVariable(), True ):
					'''
					print("{} ({}) does NOT depend on {}".format(op.toLatex(), op.__class__.__name__, expr.getVariable().toLatex()))
					if op.__class__ == SHCoefficient:
						print("numChilds={}".format(op.numChildren()))
						for j in range(op.numChildren()):
							c = op.getChildren(j)
							print("\t{}  {}".format(c.toLatex(), c.__class__.__name__))
							if c.__class__ == Tensor:
								for k in range(c.numChildren()):
									c2 = c.getChildren(k)
									print("\t\t{}  {}".format(c2.toLatex(), c2.__class__.__name__))
					'''					
					non_dependent_factors.append(op)
				else:
					dependent_factors.append(op)

			if not non_dependent_factors:
				# nothing can be extracted...keep everything as is
				return expr

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
	def visit_Summation(self, expr):
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

			if not non_dependent_factors:
				# nothing can be extracted...keep everything as is
				return expr

			# now we extract all non dependent factors and set only the remaining factors
			if len(dependent_factors) ==0:
				expr.setChildren(0, num(1))
			elif len(dependent_factors) ==1:
				expr.setChildren(0, dependent_factors[0])
			else:
				expr.setChildren(0, Multiplication(dependent_factors))

			#print("-------")
			#print(non_dependent_factors)
			#print(expr)
			return Multiplication( non_dependent_factors + [expr] )

		# do nothing
		return expr

class MergeQuotients(object):
	def visit_Multiplication(self, expr):
		numQuotients = 0
		quotients_numerators = []
		quotients_denominators = []
		others = []
		for i in range(expr.numOperands()):
			op = expr.getOperand(i)
			if op.__class__ == Quotient:
				quotients_numerators.append(op.getNumerator())
				quotients_denominators.append(op.getDenominator())
				numQuotients += 1
				continue
			others.append(op)

		if numQuotients == 0:
			return expr

		if len(quotients_numerators) == 1:
			num = quotients_numerators[0]
		else:
			num = Multiplication(quotients_numerators)

		if len(quotients_denominators) == 1:
			denom = quotients_denominators[0]
		else:
			denom = Multiplication(quotients_denominators)

		return Multiplication( [div(num,denom)] + others )


class SummationOverKronecker(object):
	def visit_Summation(self, expr):
		var = expr.getVariable()
		#print(expr.getVariable().__class__.__name__)
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

class Evaluate(object):
	def __init__(self, variables, functions):
		self.vars = variables
		self.funcs = functions

	def visit_Number(self, expr):
		return expr.getValue()
	def visit_Variable(self, expr):
		#print(expr.toLatex())
		if expr.getSymbol() in self.vars:
			return self.vars[expr.getSymbol()]
		else:
			raise ValueError("not a number")
	def visit_Negate(self, expr):
		#if expr.getOperand(0).__class__ == Number:
		return -expr.getOperand(0)
		#raise ValueError("not a number")
	def visit_Addition(self, expr):
		sum = 0
		for i in range(expr.numOperands()):
			op = expr.getOperand(i)
			#if op.__class__ == Number:
			sum += op
			#else:
			#	raise ValueError("not a number")
		return sum
	def visit_Multiplication(self, expr):
		prod = 1
		for i in range(expr.numOperands()):
			op = expr.getOperand(i)
			#if op.__class__ == Number:
			prod *= op
			#else:
			#	raise ValueError("not a number")
		return prod
	def visit_Quotient(self, expr):
		num = expr.getNumerator()
		denom = expr.getDenominator()
		#if num.__class__ == Number and denom.__class__ == Number:
		return num/denom
		#else:
		#	raise ValueError("not a number")
	def visit_ImaginaryUnit(self, expr):
		return complex(0,1)
	def visit_Sqrt(self, expr):
		#if expr.getOperand(0).__class__ == Number:
		return math.sqrt(expr.getOperand(0))
		#raise ValueError("not a number")
	def visit_Function(self, expr):

		if expr.getSymbol() in self.funcs:
			return self.funcs[expr.getSymbol()]( *expr.getArguments() )
		raise ValueError("function {} not defined for evaluation".format(expr.getSymbol()))
		'''
		# check if function has a function body which enables evaluation
		if expr.getFunctionBody() == None:
			raise ValueError("can not evaluate function without body")
		body_expr = expr.getFunctionBody().deep_copy()

		# check that all argument expressions are numbers
		numArgs = expr.numArguments()
		#for i in range(numArgs):
		#	if not expr.getArgument(i).__class__ == Number:
		#		raise ValueError()

		# all function args are numbers
		# we now take the function body and replace the argument symbols
		# with their evaluated numbers
		for i in range(numArgs):
			body_expr = apply_recursive(body_expr, Substitute(var(expr.getArgumentSymbol(i)), expr.getArgument(i)))

		# now all the variables within the function body have been replaced by the numbers from
		# the argument expressions. Lets evaluate the body
		#result = eval( body_expr )

		if result == None:
			raise ValueError("unable to evaluate function body")

		#print("function {} returned value:{}".format(expr.getSymbol(), result))

		return result
		'''
		return None


class ExpandDotProduct(object):
	# todo: generalize this to all kinds of tensors
	def visit_DotProduct(self, expr):
		l = expr.getLeft()
		r = expr.getRight()
		nc = l.numComponents()

		if l.numComponents() != r.numComponents():
			raise ValueError("number components has to match")
		terms = []
		for i in range(nc):
			terms.append(mul(l.getComponent(i).deep_copy(), r.getComponent(i).deep_copy()))

		return Addition( terms )




def apply( expr, cls, visitor ):
	visitor_function = None
	try:
		#print("asdasd {}".format(visitor.__class__.__name__))
		visitor_function = getattr(visitor.__class__, "visit_{}".format(cls.__name__))
	except AttributeError:
		pass

	if visitor_function:
		# there is a visitor function:call it
		return visitor_function(visitor, expr)

	# Iterates the class hierarchy of the Expression class (and its subclasses).
	# Tries to call visitor function for each subclass name (starting from the leave).
	bases = cls.__bases__
	for basecls in bases:
	#for i in range(len(bases)):
		#basecls = bases[i]
		# recursve up within the class hierarchy
		result = apply( expr, basecls, visitor )
		# we are done as soon as this function returns something valid
		if not result is None:
			return result

	# our visitor didnt have a visit function for the given expression type
	# or any of its base functions
	return expr


class Substitute(object):
	def __init__(self, expr, replacement):
		self.expr = expr
		self.replacement = replacement

	def visit_Expression(self, expr):
		if self.expr == expr:
			return self.replacement.deep_copy()
		return expr

def apply_recursive( expr, visitor ):
	for i in range(expr.numChildren()):
		expr.setChildren(i, apply_recursive(expr.getChildren(i), visitor))
	return apply(expr, expr.__class__, visitor)



if __name__ == "__main__":




	# setup equation

	#omega = Tensor("\\omega", rank=1, dimension=3)
	omega = var("\\omega")
	omega_x = tensor_component("\\omega", "x")
	omega_y = tensor_component("\\omega", "y")
	omega_z = tensor_component("\\omega", "z")
	dx = var("\\partial_x")
	dy = var("\\partial_y")
	dz = var("\\partial_z")
	L = fun( "L", var("\\vec{x}"), omega)
	# expression for the expanded radiance field
	L_expanded = sum( sum( mul( SHCoefficient( "L", var("l"), var("m"), var("\\vec{x}") ), SHBasis(var("l"), var("m"), omega, conjugate_complex=False) ), var('m'), neg(var('l')), var('l') ), var('l'), num(0), infty() )

	Ylm = SHBasis(var("l'"), var("m'"), omega, conjugate_complex=True)

	omega_gradient_L = add( mul(omega_x, dx, L), mul(omega_y, dy, L), mul(omega_z, dz, L))
	#expr = mul(dx, integrate( mul(mul( omega_x, Ylm ), omega_gradient_L), omega ))

	expr_x = mul(neg(dx), integrate( mul(mul( omega_x, Ylm ), omega_gradient_L), omega ))
	expr_y = mul(neg(dy), integrate( mul(mul( omega_y, Ylm ), omega_gradient_L), omega ))
	expr_z = mul(neg(dz), integrate( mul(mul( omega_z, Ylm ), omega_gradient_L), omega ))
	expr = add(expr_x, expr_y, expr_z)



	#expr = tensor("\\omega", rank=2, dimension=3)
	#expr.setComponent("xz", L)


	#print("\n----------------------------\n")
	#print("$$\n" + latex(expr) + "\n$$")

	expr = apply_recursive(expr, SHRecursiveRelation())
	expr = apply_recursive(expr, DistributiveLaw())
	expr = apply_recursive(expr, DistributiveLaw())
	expr = apply_recursive(expr, SplitIntegrals())
	expr = apply_recursive(expr, CleanupSigns())
	expr = apply_recursive(expr, Factorize())
	expr = apply_recursive(expr, SHRecursiveRelation())
	expr = apply_recursive(expr, FoldConstants())
	expr = apply_recursive(expr, Factorize())
	expr = apply_recursive(expr, DistributiveLaw())
	expr = apply_recursive(expr, SplitIntegrals())
	expr = apply_recursive(expr, CleanupSigns())
	expr = apply_recursive(expr, Factorize())
	expr = apply_recursive(expr, DistributiveLaw())
	expr = apply_recursive(expr, CleanupSigns())
	expr = apply_recursive(expr, Substitute(L, L_expanded))
	expr = apply_recursive(expr, SwitchDomains())
	expr = apply_recursive(expr, SwitchDomains())
	expr = apply_recursive(expr, Factorize())
	expr = apply_recursive(expr, SHOrthogonalityProperty())
	expr = apply_recursive(expr, SummationOverKronecker())
	expr = apply_recursive(expr, CleanupSigns())
	expr = apply_recursive(expr, MergeQuotients())
	expr = apply_recursive(expr, FoldConstants())
	expr = apply_recursive(expr, ImaginaryUnitProperty())
	expr = apply_recursive(expr, CleanupSigns())
	

	print("\n----------------------------\n")
	print("$$\n" + latex(expr) + "\n$$")
	exit(0)

	# factorize final expression according to unknowns
	numTerms = expr.numOperands()
	print( "numTerms={}".format(numTerms) )
	L_coefficients = {}
	for i in range(numTerms):
		term = expr.getOperand(i)
		is_negative = False

		if term.__class__ == Negate:
			is_negative = True
			term = term.getOperand(0)

		if term.__class__ != Multiplication:
			raise ValueError
		#print("achtung!> negative werte muessen beruecksichtigt werden!!!")
		#print("----------- term={} {}".format(i, hierarchy(term)))

		# find sh coefficient L^lm and its coefficients
		Llm = None
		coefficients = []
		partials = []

		for j in range(term.numOperands()):
			f = term.getOperand(j)
			if f.__class__ == SHCoefficient:
				if Llm != None:
					raise ValueError
				Llm = f
				continue
			elif f.__class__ == Variable and f.getSymbol().startswith("\\partial"):
				if f.getSymbol() == "\\partial_x":
					partials.append("dx")
				elif f.getSymbol() == "\\partial_y":
					partials.append("dy")
				elif f.getSymbol() == "\\partial_z":
					partials.append("dz")
				else:
					raise ValueError
				continue
			coefficients.append(f)

		# sanity check
		if Llm == None:
			raise ValueError
		if len(coefficients) == 0:
			raise ValueError
		if len(partials) != 2:
			raise ValueError

		# find the l,m arguments
		# this is tricky as we need to inspect the expressions
		# alternatively we could use the latex code as hash, but
		# I somehow didn't like it really
		l = Llm.getArgument(0)
		l_offset = 0
		m = Llm.getArgument(1)
		m_offset = 0

		if l.__class__ == Variable:
			pass
		elif l.__class__ == Addition and l.getOperand(1).__class__ == Number:
			l_offset = l.getOperand(1).getValue()
		elif l.__class__ == Addition and l.getOperand(1).__class__ == Negate and l.getOperand(1).getOperand(0).__class__ == Number:
			l_offset = -l.getOperand(1).getOperand(0).getValue()
		else:
			raise ValueError


		if m.__class__ == Variable:
			pass
		elif m.__class__ == Addition and m.getOperand(1).__class__ == Number:
			m_offset = m.getOperand(1).getValue()
		elif m.__class__ == Addition and m.getOperand(1).__class__ == Negate and m.getOperand(1).getOperand(0).__class__ == Number:
			m_offset = -m.getOperand(1).getOperand(0).getValue()
		else:
			raise ValueError

		key = (partials[0], partials[1], l_offset, m_offset)
		if not key in L_coefficients:
			L_coefficients[key] = []

		if is_negative:
			L_coefficients[key].append( neg(Multiplication(coefficients)) )
		else:
			L_coefficients[key].append( Multiplication(coefficients) )
		#L_coefficients[key] += 1

	#'''
	#print(L_coefficients)

	# ====================================================================
	# here we build the M matrices representing the second order transport operator ---------------
	partials = ["dx", "dy", "dz"]

	# this is the truncation level of the SH expansion
	order = 2
	# the number of SH coefficients per voxel
	numSHCoefficients = (order + 1) * (order + 1)

	# here we setup the M matrices for each second order derivative
	M_matrices = {}
	for d0 in partials:
		for d1 in partials:
			key = (d0, d1)
			M = [[None for j in range(numSHCoefficients)] for i in range(numSHCoefficients)]
			#M_matrices[key] = M
			M_matrices[key] = np.zeros((numSHCoefficients, numSHCoefficients), dtype = complex)
	
	
	# now we fill all M matrices with the coefficient expressions
	for d0 in partials:
		for d1 in partials:
			#print( "{} {} ---------------------------".format(d0,d1) )
			M = M_matrices[(d0, d1)]
			# now we iterate over all rows of M
			for i in range(numSHCoefficients):
				# find the l,m SH indices associated with the current row
				l,m = shtools.lmIndex(i)
				#print( "i={} l={} m={}".format(i,l,m) )

				# now we iterate all offsets
				for l_offset in [-2,0,2]:
					for m_offset in [-2,-1,0,1,2]:
						#the offset to l,m defines the column
						# offset 0,0 means: we want the main diagonal entry L_lm
						j = shtools.shIndex( l+l_offset, m+m_offset )
						# check bounds
						if j == None or j >=numSHCoefficients:
							continue

						#print( "j={} l_offset={} m_offset={}".format(j, l_offset, m_offset) )
						
						# now we evaluate the associated coefficient from the equation
						variables = {"l'":l, "m'":m}
						key = (d0, d1, l_offset, m_offset)
						if key in L_coefficients:
							coeff_list = L_coefficients[key]
							expr = None
							if len(coeff_list) > 1:
								expr = Addition(coeff_list)
							else:
								expr = coeff_list[0]
							#M[i][j] = eval(expr, **variables)
							M[i, j] = eval(expr, **variables)



	data = {}
	for d0 in partials:
		for d1 in partials:
			#print( "{} {} ---------------------------".format(d0,d1) )
			M = M_matrices[(d0, d1)]
			data['M_{}{}'.format(d0[1], d1[1])] = M
	scipy.io.savemat("C:/projects/epfl/epfl17/python/sopn/data.mat", data)




	


	'''
	for d0 in partials:
		for d1 in partials:
			#print("{}{} ----".format(d0, d1))
			for l_offset in [-2,0,2]:
				for m_offset in [-2,-1,0,1,2]:
					key = (d0, d1, l_offset, m_offset)
					if key in L_coefficients:
						#count = L_coefficients[key]
						#print( "L {} {} count={}".format(l_offset, m_offset, count) )
						coeffs = L_coefficients[key]
						print( "\n{}{}L {} {} ----------".format(d0, d1, l_offset, m_offset) )
						for coeff in coeffs:
							print(latex(coeff))
						#
						#print( "{} {} L {} {} {}".format(d0, d1, l_offset, m_offset, count) )
					else:
						pass
						#print( "{} {} L {} {} -----".format(d0, d1, l_offset, m_offset) )
			#print("")
						
	'''


