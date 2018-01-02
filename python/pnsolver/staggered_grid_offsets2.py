
import numpy as np
import util
import rte_terms
import meh
import itertools
import pnsolver
import os
import stencil2


def getUnknownIndex(unknown_node):
	args = unknown_node.getExpr().getArguments()
	# args[0] -> l
	l_expr = args[0]
	# args[1] -> m
	m_expr = args[1]

	if not (l_expr.canEvaluate() and m_expr.canEvaluate()):
		raise ValueError( "getUnknownIndex: unable to evaluate expressions for l or m" )

	l = l_expr.evaluate()
	m = m_expr.evaluate()
	return l*(l+1)+m

def getTermStructureId( node ):
	if node.__class__ == stencil2.Coefficient:
		return "C"
	elif node.__class__ == stencil2.Unknown:
		index = getUnknownIndex(node)
		return node.getSymbol() + "[{}]".format(index)
	elif node.__class__ == stencil2.Derivation:
		return "d{}({})".format(node.getVariable().getSymbol(), getTermStructureId(node.child))
	elif node.__class__ == stencil2.Product:
		return "({})*({})".format(getTermStructureId(node.a), getTermStructureId(node.b))
	else:
		raise ValueError("getTermStructureId: unknown term structure node type")


def factorize(equ):
	if equ.__class__ == meh.Addition:
		terms = equ.getOperands()
	else:
		terms = [equ]

	term_structures = {}
	for term in terms:
		term_structure = stencil2.TermStructure(term, "L")
		key = getTermStructureId(term_structure.getRoot())
		if not key in term_structures:
			term_structures[key] = term_structure
		else:
			# merge new term structure into existing one
			term_structures[key].merge(term_structure)
	return term_structures
	#factorized_terms = []
	#for key, value in term_structures.items():
	#	factorized_terms.append(value.getExpr())
	#return meh.Addition(factorized_terms)

if __name__ == "__main__":
	order = 6
	#staggered = True
	#filename = "c:/projects/epfl/epfl17/cpp/pnsolver/src/stencil.cpp"

	pn_equations = rte_terms.fopn_real.transport_term(order)
	numCoeffs = len(pn_equations)

	unknown_grid_locations = {}

	# we place the first unknown manually at the cell center
	unknown_grid_locations[0] = np.array([1,1,1])

	for coeff_index in range(numCoeffs):
	#for coeff_index in range(3):
		e = pn_equations[coeff_index]

		print("\n coeff_index={} ---------------------".format(coeff_index))
		#meh.print_expr(e)

		# factorize according to unknowns
		term_structures = factorize(e)

		# now iterate all terms and analyse coefficient and unknown
		# to see on which unknowns the coefficient according to the
		# current pn-equation depends
		for key, term_structure in term_structures.items():
			root = term_structure.getRoot()
			# we now analyse the term structure. We expect a product
			# between a coefficient and a derivative
			if not root.__class__ == stencil2.Product:
				raise ValueError("expected term to be a product between coefficient and derivative")

			if not root.a.__class__ == stencil2.Coefficient:
				raise ValueError("expected term to be a product between coefficient and derivative")

			if not root.b.__class__ == stencil2.Derivation:
				raise ValueError("expected term to be a product between coefficient and derivative")

			if not root.b.getChild().__class__ == stencil2.Unknown:
				raise ValueError("expected term to be a product between coefficient and derivative")

			# extract coefficient expression and collapse constants
			# to see which one reduce to 0
			coeff_expr = root.a.getExpr()
			coeff_expr = meh.apply_recursive(coeff_expr, meh.FoldConstants2())
			coeff_expr = meh.apply_recursive(coeff_expr, meh.CleanupSigns())

			# check if the expression evaluates to 0
			if coeff_expr.canEvaluate() and np.abs(coeff_expr.evaluate()) < 1.0e-8:
				# if the coefficient is 0, then the current equ does
				# not depend on this coefficient
				continue

			# get unknown sh index 
			unknown_index = getUnknownIndex(root.b.getChild())

			# get direction of derivation
			dimension = {"x":0, "y":1, "z":2}[root.b.getVariable().getSymbol()]

			# derive the location of this unknown from looking at
			# the location of the current equation index
			step = np.array([0,0,0])
			step[dimension] = 1

			new_location = np.mod(unknown_grid_locations[coeff_index]+step, np.array([2, 2, 2]))
			#print(new_location)

			if not unknown_index in unknown_grid_locations:
				#print("placing {} at ({}, {}, {})".format(unknown_index, new_location[0], new_location[1], new_location[2]))
				unknown_grid_locations[unknown_index] = new_location
			else:
				# make sure its consistent
				#print(np.sum(new_location-unknown_grid_locations[unknown_index]))
				if np.sum(new_location-unknown_grid_locations[unknown_index]) != 0:
					raise ValueError("safasfasf")


	# generate code from dictionary
	for o in range(order):
		#numCoeffs = 
		#print("order={}".format(o))
		print("if self.order >= {}:".format(o))
		min_index = o*o
		max_index = (o + 1) * (o + 1)
		for i in range(min_index, max_index):
			loc = unknown_grid_locations[i]
			print("\tself.place_unknown( {}, ({}, {}, {}) )".format(i, loc[0], loc[1], loc[2]))
		#print("{} {}".format(min_index, max_index))

	#print(unknown_grid_locations)

	#print(pn_equations)







	#terms = []
	#terms.append(rte_terms.fopn.transport_term())
	#terms = rte_terms.splitAddition( terms )

	#analyse(order, terms)

