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


def print_expr(expr):
	#pass
	print("\n----------------------------\n")
	print("$$\n" + cas.latex(expr) + "\n$$")



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
	Q = cas.fun( "q", x, omega)
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
	print_expr(expr)
	expr = cas.apply_recursive(expr, cas.SHRecursiveRelation())
	print_expr(expr)
	expr = cas.apply_recursive(expr, cas.DistributiveLaw())
	print_expr(expr)
	expr = cas.apply_recursive(expr, cas.DistributiveLaw())
	print_expr(expr)
	expr = cas.apply_recursive(expr, cas.CleanupSigns())
	print_expr(expr)
	expr = cas.apply_recursive(expr, cas.SHRecursiveRelation())
	print_expr(expr)
	expr = cas.apply_recursive(expr, cas.MergeQuotients())
	print_expr(expr)
	expr = cas.apply_recursive(expr, cas.ImaginaryUnitProperty())
	print_expr(expr)
	expr = cas.apply_recursive(expr, cas.FoldConstants())
	print_expr(expr)
	expr = cas.apply_recursive(expr, cas.DistributiveLaw())
	print_expr(expr)
	expr = cas.apply_recursive(expr, cas.CleanupSigns())
	print_expr(expr)
	expr = cas.apply_recursive(expr, cas.SplitIntegrals())
	print_expr(expr)
	expr = cas.apply_recursive(expr, cas.SplitDerivatives())
	print_expr(expr)
	expr = cas.apply_recursive(expr, cas.Factorize())
	print_expr(expr)
	expr = cas.apply_recursive(expr, cas.CleanupSigns())
	print_expr(expr)
	expr = cas.apply_recursive(expr, cas.Substitute(L, L_expanded))
	print_expr(expr)
	expr = cas.apply_recursive(expr, cas.SwitchDomains())
	print_expr(expr)
	expr = cas.apply_recursive(expr, cas.SwitchDomains())
	print_expr(expr)
	expr = cas.apply_recursive(expr, cas.SwitchDomains())
	print_expr(expr)
	expr = cas.apply_recursive(expr, cas.Factorize())
	print_expr(expr)
	expr = cas.apply_recursive(expr, cas.SHOrthogonalityProperty())
	print_expr(expr)
	expr = cas.apply_recursive(expr, cas.SummationOverKronecker())
	print_expr(expr)
	expr = cas.apply_recursive(expr, cas.CleanupSigns())
	print_expr(expr)

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

	print_expr(expr)
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
	print_expr(expr)
	expr = cas.apply_recursive(expr, cas.SHRecursiveRelation())
	print_expr(expr)
	expr = cas.apply_recursive(expr, cas.DistributiveLaw())
	print_expr(expr)
	expr = cas.apply_recursive(expr, cas.DistributiveLaw())
	print_expr(expr)
	expr = cas.apply_recursive(expr, cas.CleanupSigns())
	print_expr(expr)
	expr = cas.apply_recursive(expr, cas.SHRecursiveRelation())
	print_expr(expr)
	expr = cas.apply_recursive(expr, cas.MergeQuotients())
	print_expr(expr)
	expr = cas.apply_recursive(expr, cas.ImaginaryUnitProperty())
	print_expr(expr)
	expr = cas.apply_recursive(expr, cas.FoldConstants())
	print_expr(expr)
	expr = cas.apply_recursive(expr, cas.DistributiveLaw())
	print_expr(expr)
	expr = cas.apply_recursive(expr, cas.CleanupSigns())
	print_expr(expr)
	expr = cas.apply_recursive(expr, cas.SplitIntegrals())
	print_expr(expr)
	expr = cas.apply_recursive(expr, cas.SplitDerivatives())
	print_expr(expr)
	expr = cas.apply_recursive(expr, cas.Factorize())
	print_expr(expr)
	expr = cas.apply_recursive(expr, cas.CleanupSigns())
	print_expr(expr)
	expr = cas.apply_recursive(expr, cas.Substitute(L, L_expanded))
	print_expr(expr)
	expr = cas.apply_recursive(expr, cas.SwitchDomains())
	print_expr(expr)
	expr = cas.apply_recursive(expr, cas.SwitchDomains())
	print_expr(expr)
	expr = cas.apply_recursive(expr, cas.SwitchDomains())
	print_expr(expr)
	expr = cas.apply_recursive(expr, cas.Factorize())
	print_expr(expr)
	expr = cas.apply_recursive(expr, cas.SHOrthogonalityProperty())
	print_expr(expr)
	expr = cas.apply_recursive(expr, cas.SummationOverKronecker())
	print_expr(expr)
	expr = cas.apply_recursive(expr, cas.CleanupSigns())
	print_expr(expr)

	return expr

def lspn_extinction_directional_derivative_term():
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

	# extinction coefficient field
	sigma_t = cas.fun( "\\sigma_t", x)
	sigma_s = cas.fun( "\\sigma_s", x)

	omegaL = cas.tensor("", rank=1, dimension=3)
	omegaL.setComponent(0, cas.mul(omega_x, L))
	omegaL.setComponent(1, cas.mul(omega_y, L))
	omegaL.setComponent(2, cas.mul(omega_z, L))

	nabla_sigma_t = cas.tensor("", rank=1, dimension=3)
	nabla_sigma_t.setComponent(0, cas.deriv(sigma_t, x.getComponent(0), is_partial = True))
	nabla_sigma_t.setComponent(1, cas.deriv(sigma_t, x.getComponent(1), is_partial = True))
	nabla_sigma_t.setComponent(2, cas.deriv(sigma_t, x.getComponent(2), is_partial = True))

	expr = cas.neg(cas.dot(omegaL, nabla_sigma_t))

	print_expr(expr)
	expr = cas.apply_recursive(expr, cas.ExpandDotProduct())
	print_expr(expr)
	expr = cas.integrate(cas.mul( cas.SHBasis(cas.var("l'"), cas.var("m'"), omega, conjugate_complex=True), expr), omega) 
	print_expr(expr)
	expr = cas.apply_recursive(expr, cas.CleanupSigns())
	print_expr(expr)
	expr = cas.apply_recursive(expr, cas.DistributiveLaw())
	print_expr(expr)
	expr = cas.apply_recursive(expr, cas.CleanupSigns())
	print_expr(expr)
	expr = cas.apply_recursive(expr, cas.SplitIntegrals())
	print_expr(expr)
	expr = cas.apply_recursive(expr, cas.CleanupSigns())
	print_expr(expr)
	expr = cas.apply_recursive(expr, cas.SHRecursiveRelation())
	print_expr(expr)
	expr = cas.apply_recursive(expr, cas.DistributiveLaw())
	print_expr(expr)
	expr = cas.apply_recursive(expr, cas.CleanupSigns())
	print_expr(expr)
	expr = cas.apply_recursive(expr, cas.SplitIntegrals())
	print_expr(expr)
	expr = cas.apply_recursive(expr, cas.CleanupSigns())
	print_expr(expr)
	expr = cas.apply_recursive(expr, cas.Factorize())
	print_expr(expr)
	expr = cas.apply_recursive(expr, cas.Substitute(L, L_expanded))
	print_expr(expr)
	expr = cas.apply_recursive(expr, cas.SwitchDomains())
	print_expr(expr)
	expr = cas.apply_recursive(expr, cas.SwitchDomains())
	print_expr(expr)
	expr = cas.apply_recursive(expr, cas.Factorize())
	print_expr(expr)
	expr = cas.apply_recursive(expr, cas.SHOrthogonalityProperty())
	print_expr(expr)
	expr = cas.apply_recursive(expr, cas.SummationOverKronecker())
	print_expr(expr)

	return expr

def lspn_squared_extinction_term():
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

	# extinction coefficient field
	sigma_t = cas.fun( "\\sigma_t", x)
	sigma_s = cas.fun( "\\sigma_s", x)

	expr = cas.mul(cas.pow(sigma_t, cas.num(2)), L)
	print_expr(expr)
	expr = cas.integrate(cas.mul( cas.SHBasis(cas.var("l'"), cas.var("m'"), omega, conjugate_complex=True), expr), omega) 
	print_expr(expr)
	expr = cas.apply_recursive(expr, cas.Substitute(L, L_expanded))
	print_expr(expr)
	expr = cas.apply_recursive(expr, cas.SwitchDomains())
	print_expr(expr)
	expr = cas.apply_recursive(expr, cas.SwitchDomains())
	print_expr(expr)
	expr = cas.apply_recursive(expr, cas.Factorize())
	print_expr(expr)
	expr = cas.apply_recursive(expr, cas.SHOrthogonalityProperty())
	print_expr(expr)
	expr = cas.apply_recursive(expr, cas.SummationOverKronecker())
	print_expr(expr)

	return expr

def lspn_directional_derivative_scattering_term():
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

	sigma_s = cas.fun( "\\sigma_s", x)

	lambda_l = cas.fun( "\\lambda", cas.var("l"), arglevel=-1)
	# TODO: use correct value
	lambda_l.body2 = lambda l:1.0

	SL_isotropic_expanded = cas.mul( sigma_s ,cas.sum( cas.sum( cas.mul( lambda_l, cas.SHCoefficient( "f_p", cas.var("l"), cas.num(0), x ), cas.SHCoefficient( "L", cas.var("l"), cas.var("m"), x ), cas.SHBasis(cas.var("l"), cas.var("m"), omega, conjugate_complex=False) ), cas.var('m'), cas.neg(cas.var('l')), cas.var('l') ), cas.var('l'), cas.num(0), cas.infty() ) )

	nabla_SL = cas.tensor("", rank=1, dimension=3)
	nabla_SL.setComponent(0, cas.deriv(SL_isotropic_expanded, x.getComponent(0), is_partial = True))
	nabla_SL.setComponent(1, cas.deriv(SL_isotropic_expanded, x.getComponent(1), is_partial = True))
	nabla_SL.setComponent(2, cas.deriv(SL_isotropic_expanded, x.getComponent(2), is_partial = True))


	expr = cas.dot(omega, nabla_SL)
	print_expr(expr)
	expr = cas.integrate(cas.mul( cas.SHBasis(cas.var("l'"), cas.var("m'"), omega, conjugate_complex=True), expr), omega) 
	##print_expr(expr)
	expr = cas.apply_recursive(expr, cas.CleanupSigns())
	#print_expr(expr)
	expr = cas.apply_recursive(expr, cas.ExpandDotProduct())
	#print_expr(expr)
	expr = cas.apply_recursive(expr, cas.DistributiveLaw())
	#print_expr(expr)
	expr = cas.apply_recursive(expr, cas.SplitIntegrals())
	#print_expr(expr)
	expr = cas.apply_recursive(expr, cas.CleanupSigns())
	##print_expr(expr)
	expr = cas.apply_recursive(expr, cas.SwitchDomains())
	##print_expr(expr)
	expr = cas.apply_recursive(expr, cas.SwitchDomains())
	##print_expr(expr)
	expr = cas.apply_recursive(expr, cas.SwitchDomains())
	#print_expr(expr)
	expr = cas.apply_recursive(expr, cas.Factorize())
	#print_expr(expr)
	expr = cas.apply_recursive(expr, cas.SHRecursiveRelation())
	#print_expr(expr)
	expr = cas.apply_recursive(expr, cas.DistributiveLaw())
	#print_expr(expr)
	expr = cas.apply_recursive(expr, cas.SplitIntegrals())
	#print_expr(expr)
	expr = cas.apply_recursive(expr, cas.CleanupSigns())
	#print_expr(expr)
	expr = cas.apply_recursive(expr, cas.Factorize())
	#print_expr(expr)
	expr = cas.apply_recursive(expr, cas.DistributiveLaw())
	#print_expr(expr)
	expr = cas.apply_recursive(expr, cas.SHOrthogonalityProperty())
	#print_expr(expr)
	expr = cas.apply_recursive(expr, cas.DistributiveLaw())
	expr = cas.apply_recursive(expr, cas.DistributiveLaw())
	#print_expr(expr)
	expr = cas.apply_recursive(expr, cas.CleanupSigns())
	#print_expr(expr)
	expr = cas.apply_recursive(expr, cas.SplitSums())
	#print_expr(expr)
	expr = cas.apply_recursive(expr, cas.CleanupSigns())
	#print_expr(expr)
	expr = cas.apply_recursive(expr, cas.SummationOverKronecker())
	#print_expr(expr)
	expr = cas.apply_recursive(expr, cas.SplitDerivatives())
	#print_expr(expr)
	expr = cas.apply_recursive(expr, cas.CleanupSigns())
	#print_expr(expr)
	expr = cas.apply_recursive(expr, cas.Factorize())
	#print_expr(expr)
	expr = cas.apply_recursive(expr, cas.ProductRule())
	print_expr(expr)
	expr = cas.apply_recursive(expr, cas.DistributiveLaw())
	print_expr(expr)
	expr = cas.apply_recursive(expr, cas.CleanupSigns())
	print_expr(expr)


	return expr


def lspn_directional_derivative_scattering_term2():
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

	sigma_s = cas.fun( "\\sigma_s", x)

	lambda_l = cas.fun( "\\lambda", cas.var("l"), arglevel=-1)
	# TODO: use correct value
	lambda_l.body2 = lambda l:1.0

	SL_isotropic_expanded = cas.mul( sigma_s ,cas.sum( cas.sum( cas.mul( lambda_l, cas.SHCoefficient( "f_p", cas.var("l"), cas.num(0), x ), cas.SHCoefficient( "L", cas.var("l"), cas.var("m"), x ), cas.SHBasis(cas.var("l"), cas.var("m"), omega, conjugate_complex=False) ), cas.var('m'), cas.neg(cas.var('l')), cas.var('l') ), cas.var('l'), cas.num(0), cas.infty() ) )

	nabla_SL = cas.tensor("", rank=1, dimension=3)
	nabla_SL.setComponent(0, cas.deriv(SL_isotropic_expanded, x.getComponent(0), is_partial = True))
	nabla_SL.setComponent(1, cas.deriv(SL_isotropic_expanded, x.getComponent(1), is_partial = True))
	nabla_SL.setComponent(2, cas.deriv(SL_isotropic_expanded, x.getComponent(2), is_partial = True))

	# NB: we should have negative omega, as we use the adjoint of the transport operator
	# however, since we move the term to the lhs the negative signs cancel out
	expr = cas.dot( omega, cas.grad(SL_isotropic_expanded) )
	#print_expr(expr)
	expr = cas.integrate(cas.mul( cas.SHBasis(cas.var("l'"), cas.var("m'"), omega, conjugate_complex=True), expr), omega) 
	#print_expr(expr)
	expr = cas.apply_recursive(expr, cas.CleanupSigns())
	#print_expr(expr)
	expr = cas.apply_recursive(expr, cas.ExpandDotProduct())
	#print_expr(expr)
	expr = cas.apply_recursive(expr, cas.DistributiveLaw())
	#print_expr(expr)
	expr = cas.apply_recursive(expr, cas.SplitIntegrals())
	#print_expr(expr)
	expr = cas.apply_recursive(expr, cas.CleanupSigns())
	#print_expr(expr)
	expr = cas.apply_recursive(expr, cas.SHRecursiveRelation())
	#print_expr(expr)
	expr = cas.apply_recursive(expr, cas.DistributiveLaw())
	#print_expr(expr)
	expr = cas.apply_recursive(expr, cas.SplitIntegrals())
	#print_expr(expr)
	expr = cas.apply_recursive(expr, cas.CleanupSigns())
	#print_expr(expr)
	expr = cas.apply_recursive(expr, cas.SwitchDomains())
	expr = cas.apply_recursive(expr, cas.SwitchDomains())
	expr = cas.apply_recursive(expr, cas.SwitchDomains())
	#print_expr(expr)
	expr = cas.apply_recursive(expr, cas.Factorize())
	#print_expr(expr)
	expr = cas.apply_recursive(expr, cas.SHOrthogonalityProperty())
	#print_expr(expr)
	expr = cas.apply_recursive(expr, cas.SummationOverKronecker())
	#print_expr(expr)
	expr = cas.apply_recursive(expr, cas.Factorize())
	#print_expr(expr)
	expr = cas.apply_recursive(expr, cas.ProductRule())
	#print_expr(expr)
	expr = cas.apply_recursive(expr, cas.DistributiveLaw())
	#print_expr(expr)
	expr = cas.apply_recursive(expr, cas.CleanupSigns())
	print_expr(expr)

	return expr


def lspn_extinction_scattering_term():
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

	lambda_l = cas.fun( "\\lambda", cas.var("l"), arglevel=-1)
	# TODO: use correct value
	lambda_l.body2 = lambda l:1.0

	SL_isotropic_expanded = cas.sum( cas.sum( cas.mul( lambda_l, cas.SHCoefficient( "f_p", cas.var("l"), cas.num(0), x ), cas.SHCoefficient( "L", cas.var("l"), cas.var("m"), x ), cas.SHBasis(cas.var("l"), cas.var("m"), omega, conjugate_complex=False) ), cas.var('m'), cas.neg(cas.var('l')), cas.var('l') ), cas.var('l'), cas.num(0), cas.infty() )

	nabla_SL = cas.tensor("", rank=1, dimension=3)
	nabla_SL.setComponent(0, cas.deriv(SL_isotropic_expanded, x.getComponent(0), is_partial = True))
	nabla_SL.setComponent(1, cas.deriv(SL_isotropic_expanded, x.getComponent(1), is_partial = True))
	nabla_SL.setComponent(2, cas.deriv(SL_isotropic_expanded, x.getComponent(2), is_partial = True))

	# extinction coefficient field
	sigma_t = cas.fun( "\\sigma_t", x)
	sigma_s = cas.fun( "\\sigma_s", x)

	# we negate to move it onto the lefthandside
	expr = cas.neg(cas.mul( sigma_t, sigma_s, SL_isotropic_expanded ))
	print_expr(expr)
	expr = cas.integrate(cas.mul( cas.SHBasis(cas.var("l'"), cas.var("m'"), omega, conjugate_complex=True), expr), omega) 
	print_expr(expr)
	expr = cas.apply_recursive(expr, cas.CleanupSigns())
	print_expr(expr)
	expr = cas.apply_recursive(expr, cas.SwitchDomains())
	expr = cas.apply_recursive(expr, cas.SwitchDomains())
	print_expr(expr)
	expr = cas.apply_recursive(expr, cas.Factorize())
	print_expr(expr)
	expr = cas.apply_recursive(expr, cas.SHOrthogonalityProperty())
	print_expr(expr)
	expr = cas.apply_recursive(expr, cas.SummationOverKronecker())
	print_expr(expr)

	return expr



def lspn_directional_derivative_source_term():
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

	# direction-dependent emission field
	Q = cas.fun( "q", x, omega)
	Q_expanded = cas.sh_expansion(Q, x, omega)

	nabla_Q = cas.tensor("", rank=1, dimension=3)
	nabla_Q.setComponent(0, cas.deriv(Q, x.getComponent(0), is_partial = True))
	nabla_Q.setComponent(1, cas.deriv(Q, x.getComponent(1), is_partial = True))
	nabla_Q.setComponent(2, cas.deriv(Q, x.getComponent(2), is_partial = True))

	expr = cas.neg(cas.dot(omega, nabla_Q))
	#print_expr(expr)
	expr = cas.apply_recursive(expr, cas.ExpandDotProduct())
	#print_expr(expr)
	expr = cas.apply_recursive(expr, cas.CleanupSigns())
	#print_expr(expr)
	expr = cas.integrate(cas.mul( cas.SHBasis(cas.var("l'"), cas.var("m'"), omega, conjugate_complex=True), expr), omega) 
	#print_expr(expr)
	expr = cas.apply_recursive(expr, cas.DistributiveLaw())
	#print_expr(expr)
	expr = cas.apply_recursive(expr, cas.CleanupSigns())
	#print_expr(expr)
	expr = cas.apply_recursive(expr, cas.SplitIntegrals())
	#print_expr(expr)
	expr = cas.apply_recursive(expr, cas.CleanupSigns())
	#print_expr(expr)
	expr = cas.apply_recursive(expr, cas.SHRecursiveRelation())
	#print_expr(expr)
	expr = cas.apply_recursive(expr, cas.DistributiveLaw())
	#print_expr(expr)
	expr = cas.apply_recursive(expr, cas.SplitIntegrals())
	#print_expr(expr)
	expr = cas.apply_recursive(expr, cas.CleanupSigns())
	#print_expr(expr)
	expr = cas.apply_recursive(expr, cas.Factorize())
	#print_expr(expr)
	expr = cas.apply_recursive(expr, cas.Substitute(Q, Q_expanded))
	#print_expr(expr)
	expr = cas.apply_recursive(expr, cas.SwitchDomains())
	expr = cas.apply_recursive(expr, cas.SwitchDomains())
	expr = cas.apply_recursive(expr, cas.SwitchDomains())
	#print_expr(expr)
	expr = cas.apply_recursive(expr, cas.Factorize())
	#print_expr(expr)
	expr = cas.apply_recursive(expr, cas.SHOrthogonalityProperty())
	#print_expr(expr)
	expr = cas.apply_recursive(expr, cas.SummationOverKronecker())

	#expr = expr.getOperand(0)
	print_expr(expr)

	return expr

def lspn_extinction_source_term():
	omega = cas.tensor("\\omega", rank=1, dimension=3)
	omega_x = omega.getComponent(0)
	omega_y = omega.getComponent(1)
	omega_z = omega.getComponent(2)

	x = cas.tensor("\\vec{x}", rank=1, dimension=3)
	x.setComponent(0, cas.var("x"))
	x.setComponent(1, cas.var("y"))
	x.setComponent(2, cas.var("z"))
	x.collapsed = True
	
	sigma_t = cas.fun( "\\sigma_t", x)

	# direction-dependent emission field
	Q = cas.fun( "q", x, omega)
	Q_expanded = cas.sh_expansion(Q, x, omega)

	expr = cas.mul(sigma_t, Q)
	print_expr(expr)
	expr = cas.integrate(cas.mul( cas.SHBasis(cas.var("l'"), cas.var("m'"), omega, conjugate_complex=True), expr), omega) 
	print_expr(expr)
	expr = cas.apply_recursive(expr, cas.Substitute(Q, Q_expanded))
	print_expr(expr)
	expr = cas.apply_recursive(expr, cas.SwitchDomains())
	expr = cas.apply_recursive(expr, cas.SwitchDomains())
	expr = cas.apply_recursive(expr, cas.SwitchDomains())
	print_expr(expr)
	expr = cas.apply_recursive(expr, cas.Factorize())
	print_expr(expr)
	expr = cas.apply_recursive(expr, cas.SHOrthogonalityProperty())
	print_expr(expr)
	expr = cas.apply_recursive(expr, cas.SummationOverKronecker())
	print_expr(expr)

	return expr

def diffusion_terms():
	x = cas.tensor("\\vec{x}", rank=1, dimension=3)
	x.setComponent(0, cas.var("x"))
	x.setComponent(1, cas.var("y"))
	x.setComponent(2, cas.var("z"))
	x.collapsed = True

	sigma_t = cas.fun("\\sigma_t", x)

	phi = cas.SHCoefficient("L", cas.num(0), cas.num(0), x )

	D = cas.frac(cas.num(1), cas.mul(cas.num(3), sigma_t))

	Q = cas.fun( "q", x)

	expr = cas.add( cas.div(cas.mul( D, cas.grad(phi) )), cas.neg(Q) )
	cas.print_expr(expr)
	expr = cas.apply_recursive(expr, cas.ExpandGradient())
	cas.print_expr(expr)
	expr = cas.apply_recursive(expr, cas.ExpandDivergence())
	#print(cas.hierarchy(expr))
	#expr = expr.getOperand(3)

	cas.print_expr(expr)
	return  expr

def point_source(pWS):
	if np.linalg.norm(pWS - np.array([3.5, 3.5])) < 0.2:
		return 1.0
	return 0.0



def build_diffusion():
	order = 0
	domain = util.Domain2D(7.0, 70)

	pnb = pnbuilder.PNBuilder(order, domain)

	pnb.add_terms(diffusion_terms())

	functions = {}
	functions["\\sigma_t"] = lambda pWS: 10.0
	functions["q"] = point_source

	(A,b) = pnb.build_global( functions )

	data = {}
	data['A'] = A
	data['b'] = b.reshape((domain.numVoxels*pnb.numCoeffs, 1))
	filename = "C:/projects/epfl/epfl17/python/sopn/data_diffusion.mat"
	scipy.io.savemat(filename, data)



if __name__ == "__main__":
	## simple diffusion
	#build_diffusion()
	#exit(0)


	order = 1
	domain = util.Domain2D(7.0, 70)

	functions = {}
	functions["\\sigma_t"] = lambda pWS: sigma_a(pWS) + sigma_s(pWS)
	functions["\\sigma_a"] = sigma_a
	functions["\\sigma_s"] = sigma_s
	functions["f_p"] = phase_shcoeffs
	functions["q"] = source_shcoeffs

	pnb = pnbuilder.PNBuilder(order, domain)

	# staggered grid (and 3 point stencil)
	pnb.place_unknown( 0, (1,1) )
	pnb.place_unknown( 1, (1,0) )
	pnb.place_unknown( 2, (0,1) )
	pnb.set_stencil_half_steps(1)

	# non-staggered grid (all unknowns placed at cell center and 5 point stencils)
	#pnb.place_unknown( 0, (1,1) )
	#pnb.place_unknown( 1, (1,1) )
	#pnb.place_unknown( 2, (1,1) )
	#pnb.set_stencil_half_steps(2)

	# first order form ---
	#pnb.add_terms(fo_transport_term()) 
	#pnb.add_terms(fo_collision_term()) 
	#pnb.add_terms(fo_scattering_term()) 
	#pnb.add_terms(fo_source_term())

	# second order form ---
	pnb.add_terms(lspn_sotransport_term())
	#pnb.add_terms(lspn_extinction_directional_derivative_term())
	pnb.add_terms(lspn_squared_extinction_term())
	#pnb.add_terms(lspn_directional_derivative_scattering_term())
	pnb.add_terms(lspn_directional_derivative_scattering_term2())
	pnb.add_terms(lspn_extinction_scattering_term())
	pnb.add_terms(lspn_directional_derivative_source_term())
	pnb.add_terms(lspn_extinction_source_term())


	(A,b) = pnb.build_global( functions )

	#print( "check: {}".format( A[6447, 6448] ) )

	#'''
	data = {}
	data['A'] = A
	data['b'] = b.reshape((domain.numVoxels*pnb.numCoeffs, 1))
	data['numCoeffs'] = pnb.numCoeffs
	#filename = "C:/projects/epfl/epfl17/python/sopn/data_firstorder_corrected_staggered.mat"
	#filename = "C:/projects/epfl/epfl17/python/sopn/data_firstorder_staggered.mat"
	#filename = "C:/projects/epfl/epfl17/python/sopn/data_lspn_corrected.mat"
	#filename = "C:/projects/epfl/epfl17/python/sopn/data_lspn_staggered.mat"
	#filename = "C:/projects/epfl/epfl17/python/sopn/data_lspn_noddt.mat"
	filename = "C:/projects/epfl/epfl17/python/sopn/data_lspn_test2.mat"
	scipy.io.savemat(filename, data)

	#test: all terms with new scattering derivative term

	#'''
	
	'''
	data = {}
	data['sigma_s'] = domain.rasterize2(functions["\\sigma_s"])
	data['sigma_a'] = domain.rasterize2(functions["\\sigma_a"])
	data['sigma_t'] = domain.rasterize2(functions["\\sigma_t"])
	data['q'] = domain.rasterize2(lambda pWS: source_shcoeffs(0,0,pWS))
	scipy.io.savemat("C:/projects/epfl/epfl17/python/sopn/problem.mat", data)
	'''


