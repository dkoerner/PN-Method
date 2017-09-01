# This module contains the terms for the least squares form of the RTE
# (see http://www.tandfonline.com/doi/abs/10.1080/00411450.2014.927364).
# These terms are discretized automatically by the PNBuilder into Ax=b form.

import numpy as np
import meh
import util




def term0_projected_expr(debug = False):
	# setup equation
	omega = meh.tensor("\\omega", rank=1, dimension=3)
	omega_x = omega.getComponent(0)
	omega_y = omega.getComponent(1)
	omega_z = omega.getComponent(2)

	x = meh.tensor("\\vec{x}", rank=1, dimension=3)
	x.setComponent(0, meh.var("x"))
	x.setComponent(1, meh.var("y"))
	x.setComponent(2, meh.var("z"))
	x.collapsed = True

	#L = meh.fun( "L", meh.var("\\vec{x}"), omega)
	L = meh.fun( "L", x, omega)
	# expression for the expanded radiance field
	L_expanded = meh.sh_expansion(L, x, omega)

	# extinction coefficient field
	sigma_t = meh.fun( "\\sigma_t", x)
	sigma_s = meh.fun( "\\sigma_s", x)

	# direction-dependent emission field
	Q = meh.fun( "Q", x, omega)
	Q_expanded = meh.sh_expansion(Q, x, omega)

	Ylm = meh.SHBasis(meh.var("l'"), meh.var("m'"), omega, conjugate_complex=True)
	dx_L = meh.deriv(L, x.getComponent(0), is_partial = True)
	dy_L = meh.deriv(L, x.getComponent(1), is_partial = True)
	dz_L = meh.deriv(L, x.getComponent(2), is_partial = True)
	omega__dot_nablaL = meh.add( meh.mul(omega_x, dx_L), meh.mul(omega_y, dy_L), meh.mul(omega_z, dz_L))

	expr_x = meh.neg(meh.deriv(meh.integrate( meh.mul( omega_x, Ylm, omega__dot_nablaL), omega ), x.getComponent(0), is_partial=True))
	expr_y = meh.neg(meh.deriv(meh.integrate( meh.mul( omega_y, Ylm, omega__dot_nablaL), omega ), x.getComponent(1), is_partial=True))
	expr_z = meh.neg(meh.deriv(meh.integrate( meh.mul( omega_z, Ylm, omega__dot_nablaL), omega ), x.getComponent(2), is_partial=True))
	#'''
	expr = meh.add(expr_x, expr_y, expr_z)

	# now expr holds the equation. we now perform a series of operations
	# to bring it into a simple factor form from which we can read of the matrix coefficients
	if debug == True:
		meh.print_expr(expr)
	expr = meh.apply_recursive(expr, meh.SHRecursiveRelation())
	if debug == True:
		meh.print_expr(expr)
	expr = meh.apply_recursive(expr, meh.DistributiveLaw())
	if debug == True:
		meh.print_expr(expr)
	expr = meh.apply_recursive(expr, meh.DistributiveLaw())
	if debug == True:
		meh.print_expr(expr)
	expr = meh.apply_recursive(expr, meh.CleanupSigns())
	if debug == True:
		meh.print_expr(expr)
	expr = meh.apply_recursive(expr, meh.SHRecursiveRelation())
	if debug == True:
		meh.print_expr(expr)
	expr = meh.apply_recursive(expr, meh.MergeQuotients())
	if debug == True:
		meh.print_expr(expr)
	expr = meh.apply_recursive(expr, meh.ImaginaryUnitProperty())
	if debug == True:
		meh.print_expr(expr)
	expr = meh.apply_recursive(expr, meh.FoldConstants())
	if debug == True:
		meh.print_expr(expr)
	expr = meh.apply_recursive(expr, meh.DistributiveLaw())
	if debug == True:
		meh.print_expr(expr)
	expr = meh.apply_recursive(expr, meh.CleanupSigns())
	if debug == True:
		meh.print_expr(expr)
	expr = meh.apply_recursive(expr, meh.SplitIntegrals())
	if debug == True:
		meh.print_expr(expr)
	expr = meh.apply_recursive(expr, meh.SplitDerivatives())
	if debug == True:
		meh.print_expr(expr)
	expr = meh.apply_recursive(expr, meh.Factorize())
	if debug == True:
		meh.print_expr(expr)
	expr = meh.apply_recursive(expr, meh.CleanupSigns())
	if debug == True:
		meh.print_expr(expr)
	expr = meh.apply_recursive(expr, meh.Substitute(L, L_expanded))
	if debug == True:
		meh.print_expr(expr)
	expr = meh.apply_recursive(expr, meh.SwitchDomains())
	if debug == True:
		meh.print_expr(expr)
	expr = meh.apply_recursive(expr, meh.SwitchDomains())
	if debug == True:
		meh.print_expr(expr)
	expr = meh.apply_recursive(expr, meh.SwitchDomains())
	if debug == True:
		meh.print_expr(expr)
	expr = meh.apply_recursive(expr, meh.Factorize())
	if debug == True:
		meh.print_expr(expr)
	expr = meh.apply_recursive(expr, meh.SHOrthogonalityProperty())
	if debug == True:
		meh.print_expr(expr)
	expr = meh.apply_recursive(expr, meh.SummationOverKronecker())
	if debug == True:
		meh.print_expr(expr)
	expr = meh.apply_recursive(expr, meh.CleanupSigns())
	if debug == True:
		meh.print_expr(expr)
	#'''

	'''
	expr =  meh.mul( omega_x, dx_L )
	print_expr(expr, debug)
	expr = meh.integrate(meh.mul( meh.SHBasis(meh.var("l'"), meh.var("m'"), omega, conjugate_complex=True), expr), omega) 
	print_expr(expr, debug)
	expr = meh.apply_recursive(expr, meh.SHRecursiveRelation())
	print_expr(expr, debug)
	expr = meh.apply_recursive(expr, meh.DistributiveLaw())
	print_expr(expr, debug)
	expr = meh.apply_recursive(expr, meh.CleanupSigns())
	print_expr(expr, debug)
	expr = meh.apply_recursive(expr, meh.SplitIntegrals())
	print_expr(expr, debug)
	expr = meh.apply_recursive(expr, meh.Substitute(L, L_expanded))
	print_expr(expr, debug)
	expr = meh.apply_recursive(expr, meh.CleanupSigns())
	print_expr(expr, debug)
	expr = meh.apply_recursive(expr, meh.Factorize())
	print_expr(expr, debug)
	expr = meh.apply_recursive(expr, meh.SwitchDomains())
	print_expr(expr, debug)
	expr = meh.apply_recursive(expr, meh.SwitchDomains())
	print_expr(expr, debug)
	expr = meh.apply_recursive(expr, meh.SwitchDomains())
	print_expr(expr, debug)
	expr = meh.apply_recursive(expr, meh.Factorize())
	print_expr(expr, debug)
	expr = meh.apply_recursive(expr, meh.SHOrthogonalityProperty())
	print_expr(expr, debug)
	expr = meh.apply_recursive(expr, meh.SummationOverKronecker())
	print_expr(expr, debug)
	'''

	'''
	expr = meh.SHCoefficient( "L", meh.add( meh.var("l'"), meh.num(1)), meh.add( meh.var("m'"), meh.num(1)), x )
	expr = meh.deriv(expr, x.getComponent(0), is_partial=True)
	expr = meh.mul( meh.num(np.sqrt(4*np.pi)), expr )
	print_expr(expr, debug)
	'''

	return expr


def term0(pWS, theta, phi, problem):
	omega = util.sphericalDirection(theta, phi)
	L = problem["L"]

	result = 0.0
	
	#'''
	for i in range(2):
		for j in range(2):
			c = 0.0
			if i ==0 and j == 0:
				c = L.dxdx(pWS, omega)
			elif i ==0 and j == 1:
				c = L.dxdy(pWS, omega)
			elif i ==1 and j == 0:
				c = L.dydx(pWS, omega)
			elif i ==1 and j == 1:
				c = L.dydy(pWS, omega)
			else:
				raise ValueError("sdadasdasdsd")
			result += -omega[i]*omega[j]*c
	#'''

	#t = L.coeff_functions[0](location)
	#print("t={}".format(t))
	#result += 
	#L.coeff_functions[0].test()
	
	#print(result)
	#result += omega[0]*L.dx(x, omega)
	#result += L.coeff_functions[shtools.shIndex(1,1)].dx(location)

	return result


def term1_projected_expr( debug = False):
	omega = meh.tensor("\\omega", rank=1, dimension=3)
	omega_x = omega.getComponent(0)
	omega_y = omega.getComponent(1)
	omega_z = omega.getComponent(2)

	x = meh.tensor("\\vec{x}", rank=1, dimension=3)
	x.setComponent(0, meh.var("x"))
	x.setComponent(1, meh.var("y"))
	x.setComponent(2, meh.var("z"))
	x.collapsed = True

	L = meh.fun( "L", x, omega)
	L_expanded = meh.sh_expansion(L, x, omega)

	# extinction coefficient field
	sigma_t = meh.fun( "\\sigma_t", x)
	sigma_s = meh.fun( "\\sigma_s", x)

	omegaL = meh.tensor("", rank=1, dimension=3)
	omegaL.setComponent(0, meh.mul(omega_x, L))
	omegaL.setComponent(1, meh.mul(omega_y, L))
	omegaL.setComponent(2, meh.mul(omega_z, L))

	nabla_sigma_t = meh.tensor("", rank=1, dimension=3)
	nabla_sigma_t.setComponent(0, meh.deriv(sigma_t, x.getComponent(0), is_partial = True))
	nabla_sigma_t.setComponent(1, meh.deriv(sigma_t, x.getComponent(1), is_partial = True))
	nabla_sigma_t.setComponent(2, meh.deriv(sigma_t, x.getComponent(2), is_partial = True))

	expr = meh.neg(meh.dot(omegaL, nabla_sigma_t))

	if debug == True:
		meh.print_expr(expr)
	expr = meh.apply_recursive(expr, meh.ExpandDotProduct())
	if debug == True:
		meh.print_expr(expr)
	expr = meh.integrate(meh.mul( meh.SHBasis(meh.var("l'"), meh.var("m'"), omega, conjugate_complex=True), expr), omega) 
	if debug == True:
		meh.print_expr(expr)
	expr = meh.apply_recursive(expr, meh.CleanupSigns())
	if debug == True:
		meh.print_expr(expr)
	expr = meh.apply_recursive(expr, meh.DistributiveLaw())
	if debug == True:
		meh.print_expr(expr)
	expr = meh.apply_recursive(expr, meh.CleanupSigns())
	if debug == True:
		meh.print_expr(expr)
	expr = meh.apply_recursive(expr, meh.SplitIntegrals())
	if debug == True:
		meh.print_expr(expr)
	expr = meh.apply_recursive(expr, meh.CleanupSigns())
	if debug == True:
		meh.print_expr(expr)
	expr = meh.apply_recursive(expr, meh.SHRecursiveRelation())
	if debug == True:
		meh.print_expr(expr)
	expr = meh.apply_recursive(expr, meh.DistributiveLaw())
	if debug == True:
		meh.print_expr(expr)
	expr = meh.apply_recursive(expr, meh.CleanupSigns())
	if debug == True:
		meh.print_expr(expr)
	expr = meh.apply_recursive(expr, meh.SplitIntegrals())
	if debug == True:
		meh.print_expr(expr)
	expr = meh.apply_recursive(expr, meh.CleanupSigns())
	if debug == True:
		meh.print_expr(expr)
	expr = meh.apply_recursive(expr, meh.Factorize())
	if debug == True:
		meh.print_expr(expr)
	expr = meh.apply_recursive(expr, meh.Substitute(L, L_expanded))
	if debug == True:
		meh.print_expr(expr)
	expr = meh.apply_recursive(expr, meh.SwitchDomains())
	if debug == True:
		meh.print_expr(expr)
	expr = meh.apply_recursive(expr, meh.SwitchDomains())
	if debug == True:
		meh.print_expr(expr)
	expr = meh.apply_recursive(expr, meh.Factorize())
	if debug == True:
		meh.print_expr(expr)
	expr = meh.apply_recursive(expr, meh.SHOrthogonalityProperty())
	if debug == True:
		meh.print_expr(expr)
	expr = meh.apply_recursive(expr, meh.SummationOverKronecker())
	if debug == True:
		meh.print_expr(expr)


	return expr

def term1( pWS, theta, phi, problem ):
	sigma_t = problem["\\sigma_t"]
	L = problem["L"]
	omega = util.sphericalDirection(theta, phi)
	return -(omega[0]*sigma_t.dx(pWS) + omega[1]*sigma_t.dy(pWS) + omega[2]*sigma_t.dz(pWS))*L(pWS, omega)


def term2_projected_expr( debug = False ):
	omega = meh.tensor("\\omega", rank=1, dimension=3)
	omega_x = omega.getComponent(0)
	omega_y = omega.getComponent(1)
	omega_z = omega.getComponent(2)

	x = meh.tensor("\\vec{x}", rank=1, dimension=3)
	x.setComponent(0, meh.var("x"))
	x.setComponent(1, meh.var("y"))
	x.setComponent(2, meh.var("z"))
	x.collapsed = True

	L = meh.fun( "L", x, omega)
	L_expanded = meh.sh_expansion(L, x, omega)

	# extinction coefficient field
	sigma_t = meh.fun("\\sigma_t", x)
	sigma_s = meh.fun("\\sigma_s", x)

	expr = meh.mul(meh.pow(sigma_t, meh.num(2)), L)
	#print_expr(expr, debug)
	expr = meh.integrate(meh.mul( meh.SHBasis(meh.var("l'"), meh.var("m'"), omega, conjugate_complex=True), expr), omega) 
	#print_expr(expr, debug)
	expr = meh.apply_recursive(expr, meh.Substitute(L, L_expanded))
	#print_expr(expr, debug)
	expr = meh.apply_recursive(expr, meh.SwitchDomains())
	#print_expr(expr, debug)
	expr = meh.apply_recursive(expr, meh.SwitchDomains())
	#print_expr(expr, debug)
	expr = meh.apply_recursive(expr, meh.Factorize())
	#print_expr(expr, debug)
	expr = meh.apply_recursive(expr, meh.SHOrthogonalityProperty())
	#print_expr(expr, debug)
	expr = meh.apply_recursive(expr, meh.SummationOverKronecker())
	#print_expr(expr, debug)

	return expr

def term2( pWS, theta, phi, problem ):
	sigma_t = problem["\\sigma_t"](pWS)
	L = problem["L"]
	omega = util.sphericalDirection(theta, phi)
	return sigma_t*sigma_t*L(pWS, omega)


def term3_projected_expr(debug = False):
	omega = meh.tensor("\\omega", rank=1, dimension=3)
	omega_x = omega.getComponent(0)
	omega_y = omega.getComponent(1)
	omega_z = omega.getComponent(2)

	x = meh.tensor("\\vec{x}", rank=1, dimension=3)
	x.setComponent(0, meh.var("x"))
	x.setComponent(1, meh.var("y"))
	x.setComponent(2, meh.var("z"))
	x.collapsed = True

	L = meh.fun( "L", x, omega)

	sigma_s = meh.fun( "\\sigma_s", x)

	lambda_l = meh.fun( "\\lambda", meh.var("l"), arglevel=-1)
	# TODO: use correct value
	lambda_l.body2 = lambda l:1.0

	SL_isotropic_expanded = meh.mul( sigma_s ,meh.sum( meh.sum( meh.mul( lambda_l, meh.SHCoefficient( "f_p", meh.var("l"), meh.num(0), x ), meh.SHCoefficient( "L", meh.var("l"), meh.var("m"), x ), meh.SHBasis(meh.var("l"), meh.var("m"), omega, conjugate_complex=False) ), meh.var('m'), meh.neg(meh.var('l')), meh.var('l') ), meh.var('l'), meh.num(0), meh.infty() ) )

	nabla_SL = meh.tensor("", rank=1, dimension=3)
	nabla_SL.setComponent(0, meh.deriv(SL_isotropic_expanded, x.getComponent(0), is_partial = True))
	nabla_SL.setComponent(1, meh.deriv(SL_isotropic_expanded, x.getComponent(1), is_partial = True))
	nabla_SL.setComponent(2, meh.deriv(SL_isotropic_expanded, x.getComponent(2), is_partial = True))


	expr = meh.dot(omega, nabla_SL)
	#print_expr(expr)
	expr = meh.integrate(meh.mul( meh.SHBasis(meh.var("l'"), meh.var("m'"), omega, conjugate_complex=True), expr), omega) 
	##print_expr(expr)
	expr = meh.apply_recursive(expr, meh.CleanupSigns())
	#print_expr(expr)
	expr = meh.apply_recursive(expr, meh.ExpandDotProduct())
	#print_expr(expr)
	expr = meh.apply_recursive(expr, meh.DistributiveLaw())
	#print_expr(expr)
	expr = meh.apply_recursive(expr, meh.SplitIntegrals())
	#print_expr(expr)
	expr = meh.apply_recursive(expr, meh.CleanupSigns())
	##print_expr(expr)
	expr = meh.apply_recursive(expr, meh.SwitchDomains())
	##print_expr(expr)
	expr = meh.apply_recursive(expr, meh.SwitchDomains())
	##print_expr(expr)
	expr = meh.apply_recursive(expr, meh.SwitchDomains())
	#print_expr(expr)
	expr = meh.apply_recursive(expr, meh.Factorize())
	#print_expr(expr)
	expr = meh.apply_recursive(expr, meh.SHRecursiveRelation())
	#print_expr(expr)
	expr = meh.apply_recursive(expr, meh.DistributiveLaw())
	#print_expr(expr)
	expr = meh.apply_recursive(expr, meh.SplitIntegrals())
	#print_expr(expr)
	expr = meh.apply_recursive(expr, meh.CleanupSigns())
	#print_expr(expr)
	expr = meh.apply_recursive(expr, meh.Factorize())
	#print_expr(expr)
	expr = meh.apply_recursive(expr, meh.DistributiveLaw())
	#print_expr(expr)
	expr = meh.apply_recursive(expr, meh.SHOrthogonalityProperty())
	#print_expr(expr)
	expr = meh.apply_recursive(expr, meh.DistributiveLaw())
	expr = meh.apply_recursive(expr, meh.DistributiveLaw())
	#print_expr(expr)
	expr = meh.apply_recursive(expr, meh.CleanupSigns())
	#print_expr(expr)
	expr = meh.apply_recursive(expr, meh.SplitSums())
	#print_expr(expr)
	expr = meh.apply_recursive(expr, meh.CleanupSigns())
	#print_expr(expr)
	expr = meh.apply_recursive(expr, meh.SummationOverKronecker())
	#print_expr(expr)
	expr = meh.apply_recursive(expr, meh.SplitDerivatives())
	#print_expr(expr)
	expr = meh.apply_recursive(expr, meh.CleanupSigns())
	#print_expr(expr)
	expr = meh.apply_recursive(expr, meh.Factorize())
	#print_expr(expr)
	expr = meh.apply_recursive(expr, meh.ProductRule())
	#print_expr(expr)
	expr = meh.apply_recursive(expr, meh.DistributiveLaw())
	#print_expr(expr)
	expr = meh.apply_recursive(expr, meh.CleanupSigns())
	#print_expr(expr)
	expr = meh.apply_recursive(expr, meh.CleanupSigns())
	#print_expr(expr, debug)


	return expr


def term3( pWS, theta, phi, problem ):
	# NB: assuming constant phase function....
	omega = util.sphericalDirection(theta, phi)
	sigma_t = problem["\\sigma_t"]
	sigma_s = problem["\\sigma_s"]

	L00 = problem["L"].get_coeff(0)

	result = 0.0
	result += omega[0]*(sigma_s.dx(pWS)*L00(pWS) + sigma_s(pWS)*L00.dx(pWS))
	result += omega[1]*(sigma_s.dy(pWS)*L00(pWS) + sigma_s(pWS)*L00.dy(pWS))
	result += omega[2]*(sigma_s.dz(pWS)*L00(pWS) + sigma_s(pWS)*L00.dz(pWS))
	result *= np.sqrt(4.0*np.pi)/(4.0*np.pi)

	return result


def term4_projected_expr( debug = False):
	omega = meh.tensor("\\omega", rank=1, dimension=3)
	omega_x = omega.getComponent(0)
	omega_y = omega.getComponent(1)
	omega_z = omega.getComponent(2)

	x = meh.tensor("\\vec{x}", rank=1, dimension=3)
	x.setComponent(0, meh.var("x"))
	x.setComponent(1, meh.var("y"))
	x.setComponent(2, meh.var("z"))
	x.collapsed = True

	L = meh.fun( "L", x, omega)

	lambda_l = meh.fun( "\\lambda", meh.var("l"), arglevel=-1)
	# TODO: use correct value
	lambda_l.body2 = lambda l:1.0

	SL_isotropic_expanded = meh.sum( meh.sum( meh.mul( lambda_l, meh.SHCoefficient( "f_p", meh.var("l"), meh.num(0), x ), meh.SHCoefficient( "L", meh.var("l"), meh.var("m"), x ), meh.SHBasis(meh.var("l"), meh.var("m"), omega, conjugate_complex=False) ), meh.var('m'), meh.neg(meh.var('l')), meh.var('l') ), meh.var('l'), meh.num(0), meh.infty() )

	nabla_SL = meh.tensor("", rank=1, dimension=3)
	nabla_SL.setComponent(0, meh.deriv(SL_isotropic_expanded, x.getComponent(0), is_partial = True))
	nabla_SL.setComponent(1, meh.deriv(SL_isotropic_expanded, x.getComponent(1), is_partial = True))
	nabla_SL.setComponent(2, meh.deriv(SL_isotropic_expanded, x.getComponent(2), is_partial = True))

	# extinction coefficient field
	sigma_t = meh.fun( "\\sigma_t", x)
	sigma_s = meh.fun( "\\sigma_s", x)

	# we negate to move it onto the lefthandside
	expr = meh.neg(meh.mul( sigma_t, sigma_s, SL_isotropic_expanded ))
	#print_expr(expr,debug)
	expr = meh.integrate(meh.mul( meh.SHBasis(meh.var("l'"), meh.var("m'"), omega, conjugate_complex=True), expr), omega) 
	#print_expr(expr,debug)
	expr = meh.apply_recursive(expr, meh.CleanupSigns())
	#print_expr(expr,debug)
	expr = meh.apply_recursive(expr, meh.SwitchDomains())
	expr = meh.apply_recursive(expr, meh.SwitchDomains())
	#print_expr(expr,debug)
	expr = meh.apply_recursive(expr, meh.Factorize())
	#print_expr(expr,debug)
	expr = meh.apply_recursive(expr, meh.SHOrthogonalityProperty())
	#print_expr(expr,debug)
	expr = meh.apply_recursive(expr, meh.SummationOverKronecker())
	#print_expr(expr,debug)

	return expr
	'''
	omega = meh.tensor("\\omega", rank=1, dimension=3)
	omega_x = omega.getComponent(0)
	omega_y = omega.getComponent(1)
	omega_z = omega.getComponent(2)

	x = meh.tensor("\\vec{x}", rank=1, dimension=3)
	x.setComponent(0, meh.var("x"))
	x.setComponent(1, meh.var("y"))
	x.setComponent(2, meh.var("z"))
	x.collapsed = True

	L = meh.fun( "L", x, omega)

	sigma_s = meh.fun( "\\sigma_s", x)

	lambda_l = meh.fun( "\\lambda", meh.var("l"), arglevel=-1)
	# TODO: use correct value
	lambda_l.body2 = lambda l:1.0

	SL_isotropic_expanded = meh.mul( sigma_s ,meh.sum( meh.sum( meh.mul( lambda_l, meh.SHCoefficient( "f_p", meh.var("l"), meh.num(0), x ), meh.SHCoefficient( "L", meh.var("l"), meh.var("m"), x ), meh.SHBasis(meh.var("l"), meh.var("m"), omega, conjugate_complex=False) ), meh.var('m'), meh.neg(meh.var('l')), meh.var('l') ), meh.var('l'), meh.num(0), meh.infty() ) )

	nabla_SL = meh.tensor("", rank=1, dimension=3)
	nabla_SL.setComponent(0, meh.deriv(SL_isotropic_expanded, x.getComponent(0), is_partial = True))
	nabla_SL.setComponent(1, meh.deriv(SL_isotropic_expanded, x.getComponent(1), is_partial = True))
	nabla_SL.setComponent(2, meh.deriv(SL_isotropic_expanded, x.getComponent(2), is_partial = True))

	# NB: we should have negative omega, as we use the adjoint of the transport operator
	# however, since we move the term to the lhs the negative signs cancel out
	expr = meh.dot( omega, meh.grad(SL_isotropic_expanded) )
	#print_expr(expr,debug)
	expr = meh.integrate(meh.mul( meh.SHBasis(meh.var("l'"), meh.var("m'"), omega, conjugate_complex=True), expr), omega) 
	#print_expr(expr,debug)
	expr = meh.apply_recursive(expr, meh.CleanupSigns())
	#print_expr(expr,debug)
	expr = meh.apply_recursive(expr, meh.ExpandDotProduct())
	#print_expr(expr,debug)
	expr = meh.apply_recursive(expr, meh.DistributiveLaw())
	#print_expr(expr,debug)
	expr = meh.apply_recursive(expr, meh.SplitIntegrals())
	#print_expr(expr,debug)
	expr = meh.apply_recursive(expr, meh.CleanupSigns())
	#print_expr(expr,debug)
	expr = meh.apply_recursive(expr, meh.SHRecursiveRelation())
	#print_expr(expr,debug)
	expr = meh.apply_recursive(expr, meh.DistributiveLaw())
	#print_expr(expr,debug)
	expr = meh.apply_recursive(expr, meh.SplitIntegrals())
	#print_expr(expr,debug)
	expr = meh.apply_recursive(expr, meh.CleanupSigns())
	#print_expr(expr,debug)
	expr = meh.apply_recursive(expr, meh.SwitchDomains())
	expr = meh.apply_recursive(expr, meh.SwitchDomains())
	expr = meh.apply_recursive(expr, meh.SwitchDomains())
	#print_expr(expr,debug)
	expr = meh.apply_recursive(expr, meh.Factorize())
	#print_expr(expr,debug)
	expr = meh.apply_recursive(expr, meh.SHOrthogonalityProperty())
	#print_expr(expr,debug)
	expr = meh.apply_recursive(expr, meh.SummationOverKronecker())
	#print_expr(expr,debug)
	expr = meh.apply_recursive(expr, meh.Factorize())
	#print_expr(expr,debug)
	expr = meh.apply_recursive(expr, meh.ProductRule())
	#print_expr(expr,debug)
	expr = meh.apply_recursive(expr, meh.DistributiveLaw())
	#print_expr(expr,debug)
	expr = meh.apply_recursive(expr, meh.CleanupSigns())
	#print_expr(expr,debug)

	return expr
	'''

def term4( pWS, theta, phi, problem ):
	# NB: assuming constant phase function....
	omega = util.sphericalDirection(theta, phi)
	sigma_t = problem["\\sigma_t"](pWS)
	sigma_s = problem["\\sigma_s"](pWS)
	L = problem["L"]
	int_L = L.integral_over_solid_angle(pWS)
	return -sigma_t*sigma_s*(1.0/(4.0*np.pi))*int_L

'''
def extinction_scattering_term(debug=False):
	omega = meh.tensor("\\omega", rank=1, dimension=3)
	omega_x = omega.getComponent(0)
	omega_y = omega.getComponent(1)
	omega_z = omega.getComponent(2)

	x = meh.tensor("\\vec{x}", rank=1, dimension=3)
	x.setComponent(0, meh.var("x"))
	x.setComponent(1, meh.var("y"))
	x.setComponent(2, meh.var("z"))
	x.collapsed = True

	L = meh.fun( "L", x, omega)

	lambda_l = meh.fun( "\\lambda", meh.var("l"), arglevel=-1)
	# TODO: use correct value
	lambda_l.body2 = lambda l:1.0

	SL_isotropic_expanded = meh.sum( meh.sum( meh.mul( lambda_l, meh.SHCoefficient( "f_p", meh.var("l"), meh.num(0), x ), meh.SHCoefficient( "L", meh.var("l"), meh.var("m"), x ), meh.SHBasis(meh.var("l"), meh.var("m"), omega, conjugate_complex=False) ), meh.var('m'), meh.neg(meh.var('l')), meh.var('l') ), meh.var('l'), meh.num(0), meh.infty() )

	nabla_SL = meh.tensor("", rank=1, dimension=3)
	nabla_SL.setComponent(0, meh.deriv(SL_isotropic_expanded, x.getComponent(0), is_partial = True))
	nabla_SL.setComponent(1, meh.deriv(SL_isotropic_expanded, x.getComponent(1), is_partial = True))
	nabla_SL.setComponent(2, meh.deriv(SL_isotropic_expanded, x.getComponent(2), is_partial = True))

	# extinction coefficient field
	sigma_t = meh.fun( "\\sigma_t", x)
	sigma_s = meh.fun( "\\sigma_s", x)

	# we negate to move it onto the lefthandside
	expr = meh.neg(meh.mul( sigma_t, sigma_s, SL_isotropic_expanded ))
	print_expr(expr,debug)
	expr = meh.integrate(meh.mul( meh.SHBasis(meh.var("l'"), meh.var("m'"), omega, conjugate_complex=True), expr), omega) 
	print_expr(expr,debug)
	expr = meh.apply_recursive(expr, meh.CleanupSigns())
	print_expr(expr,debug)
	expr = meh.apply_recursive(expr, meh.SwitchDomains())
	expr = meh.apply_recursive(expr, meh.SwitchDomains())
	print_expr(expr,debug)
	expr = meh.apply_recursive(expr, meh.Factorize())
	print_expr(expr,debug)
	expr = meh.apply_recursive(expr, meh.SHOrthogonalityProperty())
	print_expr(expr,debug)
	expr = meh.apply_recursive(expr, meh.SummationOverKronecker())
	print_expr(expr,debug)

	return expr
'''


def term5_projected_expr(debug = False):
	omega = meh.tensor("\\omega", rank=1, dimension=3)
	omega_x = omega.getComponent(0)
	omega_y = omega.getComponent(1)
	omega_z = omega.getComponent(2)

	x = meh.tensor("\\vec{x}", rank=1, dimension=3)
	x.setComponent(0, meh.var("x"))
	x.setComponent(1, meh.var("y"))
	x.setComponent(2, meh.var("z"))
	x.collapsed = True

	L = meh.fun( "L", x, omega)

	# direction-dependent emission field
	Q = meh.fun( "q", x, omega)
	Q_expanded = meh.sh_expansion(Q, x, omega)

	nabla_Q = meh.tensor("", rank=1, dimension=3)
	nabla_Q.setComponent(0, meh.deriv(Q, x.getComponent(0), is_partial = True))
	nabla_Q.setComponent(1, meh.deriv(Q, x.getComponent(1), is_partial = True))
	nabla_Q.setComponent(2, meh.deriv(Q, x.getComponent(2), is_partial = True))

	expr = meh.neg(meh.dot(omega, nabla_Q))
	#print_expr(expr)
	expr = meh.apply_recursive(expr, meh.ExpandDotProduct())
	#print_expr(expr)
	expr = meh.apply_recursive(expr, meh.CleanupSigns())
	#print_expr(expr)
	expr = meh.integrate(meh.mul( meh.SHBasis(meh.var("l'"), meh.var("m'"), omega, conjugate_complex=True), expr), omega) 
	#print_expr(expr)
	expr = meh.apply_recursive(expr, meh.DistributiveLaw())
	#print_expr(expr)
	expr = meh.apply_recursive(expr, meh.CleanupSigns())
	#print_expr(expr)
	expr = meh.apply_recursive(expr, meh.SplitIntegrals())
	#print_expr(expr)
	expr = meh.apply_recursive(expr, meh.CleanupSigns())
	#print_expr(expr)
	expr = meh.apply_recursive(expr, meh.SHRecursiveRelation())
	#print_expr(expr)
	expr = meh.apply_recursive(expr, meh.DistributiveLaw())
	#print_expr(expr)
	expr = meh.apply_recursive(expr, meh.SplitIntegrals())
	#print_expr(expr)
	expr = meh.apply_recursive(expr, meh.CleanupSigns())
	#print_expr(expr)
	expr = meh.apply_recursive(expr, meh.Factorize())
	#print_expr(expr)
	expr = meh.apply_recursive(expr, meh.Substitute(Q, Q_expanded))
	#print_expr(expr)
	expr = meh.apply_recursive(expr, meh.SwitchDomains())
	expr = meh.apply_recursive(expr, meh.SwitchDomains())
	expr = meh.apply_recursive(expr, meh.SwitchDomains())
	#print_expr(expr)
	expr = meh.apply_recursive(expr, meh.Factorize())
	#print_expr(expr)
	expr = meh.apply_recursive(expr, meh.SHOrthogonalityProperty())
	#print_expr(expr)
	expr = meh.apply_recursive(expr, meh.SummationOverKronecker())

	#expr = expr.getOperand(0)
	#print_expr(expr, debug)

	return expr

def term6_projected_expr( debug = False ):
	omega = meh.tensor("\\omega", rank=1, dimension=3)
	omega_x = omega.getComponent(0)
	omega_y = omega.getComponent(1)
	omega_z = omega.getComponent(2)

	x = meh.tensor("\\vec{x}", rank=1, dimension=3)
	x.setComponent(0, meh.var("x"))
	x.setComponent(1, meh.var("y"))
	x.setComponent(2, meh.var("z"))
	x.collapsed = True
	
	sigma_t = meh.fun( "\\sigma_t", x)

	# direction-dependent emission field
	Q = meh.fun( "q", x, omega)
	Q_expanded = meh.sh_expansion(Q, x, omega)

	expr = meh.mul(sigma_t, Q)
	#print_expr(expr, debug)
	expr = meh.integrate(meh.mul( meh.SHBasis(meh.var("l'"), meh.var("m'"), omega, conjugate_complex=True), expr), omega) 
	#print_expr(expr, debug)
	expr = meh.apply_recursive(expr, meh.Substitute(Q, Q_expanded))
	#print_expr(expr, debug)
	expr = meh.apply_recursive(expr, meh.SwitchDomains())
	expr = meh.apply_recursive(expr, meh.SwitchDomains())
	expr = meh.apply_recursive(expr, meh.SwitchDomains())
	#print_expr(expr, debug)
	expr = meh.apply_recursive(expr, meh.Factorize())
	#print_expr(expr, debug)
	expr = meh.apply_recursive(expr, meh.SHOrthogonalityProperty())
	#print_expr(expr, debug)
	expr = meh.apply_recursive(expr, meh.SummationOverKronecker())
	#print_expr(expr, debug)

	return expr











if __name__ == "__main__":
	pass