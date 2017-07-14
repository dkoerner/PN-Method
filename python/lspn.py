import numpy as np
import util
import time
import shtools
import scipy.io
import cas
import pnbuilder



# define problem ---
def source( pWS ):
	x = pWS[0]
	y = pWS[1]
	if x > 3.0 and x < 4.0 and y > 3.0 and y < 4.0:
		return 1.0
	return 0.0

def sigma_a( pWS ):
	x = pWS[0]
	y = pWS[1]
	cx = np.ceil(x)
	cy = np.ceil(y)
	g = 0
	if np.ceil((x+y)/2.0)*2.0 == (cx+cy) and cx > 1.0 and cx < 7.0 and cy > 1.0 and cy-2.0*np.abs(cx-4.0) < 4:
		g = 1
	return (1.0-g)*0 + g*10

def sigma_s( pWS ):
	x = pWS[0]
	y = pWS[1]
	cx = np.ceil(x)
	cy = np.ceil(y)
	g = 0
	if np.ceil((x+y)/2.0)*2.0 == (cx+cy) and cx > 1.0 and cx < 7.0 and cy > 1.0 and cy-2.0*np.abs(cx-4.0) < 4:
		g = 1
	return (1.0-g)*1 + g*0


def phase_shcoeffs( l, m, pWS ):
	if l == 0:
		return 1.0
	return 0.0

def source_shcoeffs( l, m, pWS ):
	if l==0:
		x = pWS[0]
		y = pWS[1]
		if x > 3.0 and x < 4.0 and y > 3.0 and y < 4.0:
			return 1.0
		return 0.0
	return 0.0


def lspn_sotransport_term():
	# setup equation
	omega = cas.tensor("\\omega", rank=1, dimension=3)
	omega_x = omega.getComponent(0)
	omega_y = omega.getComponent(1)
	omega_z = omega.getComponent(2)

	x = cas.tensor("\\vec{x}", rank=1, dimension=3)
	x.setComponent(0, cas.var("x"))
	x.setComponent(1, cas.var("y"))
	x.setComponent(2, cas.var("z"))
	x.collapsed = True

	#L = cas.fun( "L", cas.var("\\vec{x}"), omega)
	L = cas.fun( "L", x, omega)
	# expression for the expanded radiance field
	L_expanded = cas.sh_expansion(L, x, omega)

	# extinction coefficient field
	sigma_t = cas.fun( "\\sigma_t", x)
	sigma_s = cas.fun( "\\sigma_s", x)

	# direction-dependent emission field
	Q = cas.fun( "Q", x, omega)
	Q_expanded = cas.sh_expansion(Q, x, omega)

	Ylm = cas.SHBasis(cas.var("l'"), cas.var("m'"), omega, conjugate_complex=True)
	dx_L = cas.deriv(L, x.getComponent(0), is_partial = True)
	dy_L = cas.deriv(L, x.getComponent(1), is_partial = True)
	dz_L = cas.deriv(L, x.getComponent(2), is_partial = True)
	omega__dot_nablaL = cas.add( cas.mul(omega_x, dx_L), cas.mul(omega_y, dy_L), cas.mul(omega_z, dz_L))

	expr_x = cas.neg(cas.deriv(cas.integrate( cas.mul( omega_x, Ylm, omega__dot_nablaL), omega ), x.getComponent(0), is_partial=True))
	expr_y = cas.neg(cas.deriv(cas.integrate( cas.mul( omega_y, Ylm, omega__dot_nablaL), omega ), x.getComponent(1), is_partial=True))
	expr_z = cas.neg(cas.deriv(cas.integrate( cas.mul( omega_z, Ylm, omega__dot_nablaL), omega ), x.getComponent(2), is_partial=True))
	expr = cas.add(expr_x, expr_y, expr_z)

	# now expr holds the equation. we now perform a series of operations
	# to bring it into a simple factor form from which we can read of the matrix coefficients
	expr = cas.apply_recursive(expr, cas.SHRecursiveRelation())
	expr = cas.apply_recursive(expr, cas.DistributiveLaw())
	expr = cas.apply_recursive(expr, cas.DistributiveLaw())
	expr = cas.apply_recursive(expr, cas.CleanupSigns())
	expr = cas.apply_recursive(expr, cas.SHRecursiveRelation())
	expr = cas.apply_recursive(expr, cas.MergeQuotients())
	expr = cas.apply_recursive(expr, cas.ImaginaryUnitProperty())
	expr = cas.apply_recursive(expr, cas.FoldConstants())
	expr = cas.apply_recursive(expr, cas.DistributiveLaw())
	expr = cas.apply_recursive(expr, cas.CleanupSigns())
	expr = cas.apply_recursive(expr, cas.SplitIntegrals())
	expr = cas.apply_recursive(expr, cas.SplitDerivatives())
	expr = cas.apply_recursive(expr, cas.Factorize())
	expr = cas.apply_recursive(expr, cas.CleanupSigns())
	expr = cas.apply_recursive(expr, cas.Substitute(L, L_expanded))
	expr = cas.apply_recursive(expr, cas.SwitchDomains())
	expr = cas.apply_recursive(expr, cas.SwitchDomains())
	expr = cas.apply_recursive(expr, cas.SwitchDomains())
	expr = cas.apply_recursive(expr, cas.Factorize())
	expr = cas.apply_recursive(expr, cas.SHOrthogonalityProperty())
	expr = cas.apply_recursive(expr, cas.SummationOverKronecker())
	expr = cas.apply_recursive(expr, cas.CleanupSigns())

	return expr

def fo_transport_term():

	omega = cas.tensor("\\omega", rank=1, dimension=3)
	omega_x = omega.getComponent(0)
	omega_y = omega.getComponent(1)
	omega_z = omega.getComponent(2)

	x = cas.tensor("\\vec{x}", rank=1, dimension=3)
	x.setComponent(0, cas.var("x"))
	x.setComponent(1, cas.var("y"))
	x.setComponent(2, cas.var("z"))
	x.collapsed = True

	L = cas.fun( "L", x, omega)
	L_expanded = cas.sh_expansion(L, x, omega)
	#L_coeffs = cas.SHCoefficient( "L", cas.var("l"), cas.var("m"), x )

	dx_L = cas.deriv(L, x.getComponent(0), is_partial = True)
	dy_L = cas.deriv(L, x.getComponent(1), is_partial = True)
	dz_L = cas.deriv(L, x.getComponent(2), is_partial = True)
	omega_dot_nablaL = cas.add( cas.mul(omega_x, dx_L), cas.mul(omega_y, dy_L), cas.mul(omega_z, dz_L))


	expr = omega_dot_nablaL


	expr = cas.integrate(cas.mul( cas.SHBasis(cas.var("l'"), cas.var("m'"), omega, conjugate_complex=True), expr), omega) 
	expr = cas.apply_recursive(expr, cas.DistributiveLaw())
	expr = cas.apply_recursive(expr, cas.SHRecursiveRelation())
	expr = cas.apply_recursive(expr, cas.SplitIntegrals())
	expr = cas.apply_recursive(expr, cas.DistributiveLaw())
	expr = cas.apply_recursive(expr, cas.CleanupSigns())
	expr = cas.apply_recursive(expr, cas.SplitIntegrals())
	expr = cas.apply_recursive(expr, cas.CleanupSigns())
	expr = cas.apply_recursive(expr, cas.Factorize())
	expr = cas.apply_recursive(expr, cas.Substitute(L, L_expanded))
	expr = cas.apply_recursive(expr, cas.SwitchDomains())
	expr = cas.apply_recursive(expr, cas.SwitchDomains())
	expr = cas.apply_recursive(expr, cas.SwitchDomains())
	expr = cas.apply_recursive(expr, cas.Factorize())
	expr = cas.apply_recursive(expr, cas.SHOrthogonalityProperty())
	expr = cas.apply_recursive(expr, cas.SummationOverKronecker())

	#print("\n----------------------------\n")
	#print("$$\n" + cas.latex(expr) + "\n$$")
	#print(cas.hierarchy(expr))
	return expr


def fo_collision_term():
	omega = cas.tensor("\\omega", rank=1, dimension=3)
	omega_x = omega.getComponent(0)
	omega_y = omega.getComponent(1)
	omega_z = omega.getComponent(2)

	x = cas.tensor("\\vec{x}", rank=1, dimension=3)
	x.setComponent(0, cas.var("x"))
	x.setComponent(1, cas.var("y"))
	x.setComponent(2, cas.var("z"))
	x.collapsed = True

	L = cas.fun( "L", x, omega)
	L_expanded = cas.sh_expansion(L, x, omega)
	L_coeffs = cas.SHCoefficient( "L", cas.var("l'"), cas.var("m'"), x )


	sigma_t = cas.fun( "\\sigma_t", x)
	expr = cas.mul(sigma_t, L_coeffs)

	#print("\n----------------------------\n")
	#print("$$\n" + cas.latex(expr) + "\n$$")
	#print(cas.hierarchy(expr))

	return expr

def fo_scattering_term():
	omega = cas.tensor("\\omega", rank=1, dimension=3)
	omega_x = omega.getComponent(0)
	omega_y = omega.getComponent(1)
	omega_z = omega.getComponent(2)

	x = cas.tensor("\\vec{x}", rank=1, dimension=3)
	x.setComponent(0, cas.var("x"))
	x.setComponent(1, cas.var("y"))
	x.setComponent(2, cas.var("z"))
	x.collapsed = True

	L = cas.fun( "L", x, omega)
	L_expanded = cas.sh_expansion(L, x, omega)
	L_coeffs = cas.SHCoefficient( "L", cas.var("l'"), cas.var("m'"), x )

	sigma_s = cas.fun( "\\sigma_s", x)

	f_p = cas.fun( "f_p", cas.var("l'"), cas.num(0), x)
	# TODO: we skip lambda_l as this was not used in the starmap paper either
	# 		however, I think it needs to be here...

	# the negation comes from the fact that we bring this to the other side
	expr = cas.neg(cas.mul(sigma_s, f_p, L_coeffs))
	#print("\n----------------------------\n")
	#print("$$\n" + cas.latex(expr) + "\n$$")
	#print(cas.hierarchy(expr))
	return expr

def fo_source_term():
	x = cas.tensor("\\vec{x}", rank=1, dimension=3)
	x.setComponent(0, cas.var("x"))
	x.setComponent(1, cas.var("y"))
	x.setComponent(2, cas.var("z"))
	x.collapsed = True

	q = cas.fun( "q", cas.var("l'"), cas.var("m'"), x)
	expr = q
	return expr




if __name__ == "__main__":

	order = 1
	domain = util.Domain2D(7.0, 70)
	#solver = SimplePN2D(order, )
	#solver = LSPN2D(order, util.Domain2D(7.0, 50))
	#u = solver.run(sigma_a, sigma_s, source)


	pnb = pnbuilder.PNBuilder(order, domain)
	#pnb.add_lhs(lspn_sotransport_term()) 
	pnb.add_terms(fo_transport_term()) 
	pnb.add_terms(fo_collision_term()) 
	pnb.add_terms(fo_scattering_term()) 
	pnb.add_terms(fo_source_term())
	#exit(1)
	(A,b) = pnb.build_global( sigma_a, sigma_s, phase_shcoeffs, source_shcoeffs )


	data = {}
	data['A_new'] = A
	data['b_new'] = b.reshape((domain.numVoxels*pnb.numCoeffs, 1))
	scipy.io.savemat("C:/projects/epfl/epfl17/python/sopn/data_new.mat", data)


	'''
	
	pnb.add_rhs(expr)
	#	stores coefficient chain with coefficient matrices and information about partial derivatives internally

	(A, b) = pnb.build_global()


	


	'''

	#check in matlab
	#result = (abs((A == A_new)-1));
	#max(result(:)) <- should be zero

