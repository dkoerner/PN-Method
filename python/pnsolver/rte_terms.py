# This module contains the terms for the least squares form of the RTE
# (see http://www.tandfonline.com/doi/abs/10.1080/00411450.2014.927364).

import numpy as np
import meh
import util

class cda(object):
	@staticmethod
	def term0_projected_expr(debug = False):
		pass


class lspn(object):
	@staticmethod
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

		return expr

	@staticmethod
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

	@staticmethod
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

	@staticmethod
	def term1( pWS, theta, phi, problem ):
		sigma_t = problem["\\sigma_t"]
		L = problem["L"]
		omega = util.sphericalDirection(theta, phi)
		return -(omega[0]*sigma_t.dx(pWS) + omega[1]*sigma_t.dy(pWS) + omega[2]*sigma_t.dz(pWS))*L(pWS, omega)

	@staticmethod
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

	@staticmethod
	def term2( pWS, theta, phi, problem ):
		sigma_t = problem["\\sigma_t"](pWS)
		L = problem["L"]
		omega = util.sphericalDirection(theta, phi)
		return sigma_t*sigma_t*L(pWS, omega)

	@staticmethod
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

	@staticmethod
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

	@staticmethod
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

	@staticmethod
	def term4( pWS, theta, phi, problem ):
		# NB: assuming constant phase function....
		omega = util.sphericalDirection(theta, phi)
		sigma_t = problem["\\sigma_t"](pWS)
		sigma_s = problem["\\sigma_s"](pWS)
		L = problem["L"]
		int_L = L.integral_over_solid_angle(pWS)
		return -sigma_t*sigma_s*(1.0/(4.0*np.pi))*int_L


	@staticmethod
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

	@staticmethod
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





class fopn(object):
	@staticmethod
	def transport_term( debug = False ):

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
		#L_coeffs = meh.SHCoefficient( "L", meh.var("l"), meh.var("m"), x )

		dx_L = meh.deriv(L, x.getComponent(0), is_partial = True)
		dy_L = meh.deriv(L, x.getComponent(1), is_partial = True)
		dz_L = meh.deriv(L, x.getComponent(2), is_partial = True)
		omega_dot_nablaL = meh.add( meh.mul(omega_x, dx_L), meh.mul(omega_y, dy_L), meh.mul(omega_z, dz_L))


		expr = omega_dot_nablaL

		#if debug == True:
		#	meh.print_expr(expr)



		expr = meh.integrate(meh.mul( meh.SHBasis(meh.var("l'"), meh.var("m'"), omega, conjugate_complex=True), expr), omega) 
		#if debug == True:
		#	print("------------------------------")
		#	meh.print_expr(expr)

		expr = meh.apply_recursive(expr, meh.DistributiveLaw())
		#if debug == True:
		#	print("------------------------------")
		#	meh.print_expr(expr)

		expr = meh.apply_recursive(expr, meh.SHRecursiveRelation())
		expr = meh.apply_recursive(expr, meh.SplitIntegrals())
		expr = meh.apply_recursive(expr, meh.DistributiveLaw())
		expr = meh.apply_recursive(expr, meh.CleanupSigns())
		expr = meh.apply_recursive(expr, meh.SplitIntegrals())
		expr = meh.apply_recursive(expr, meh.CleanupSigns())
		expr = meh.apply_recursive(expr, meh.Factorize())
		expr = meh.apply_recursive(expr, meh.Substitute(L, L_expanded))
		expr = meh.apply_recursive(expr, meh.SwitchDomains())
		expr = meh.apply_recursive(expr, meh.SwitchDomains())
		expr = meh.apply_recursive(expr, meh.SwitchDomains())
		expr = meh.apply_recursive(expr, meh.Factorize())
		expr = meh.apply_recursive(expr, meh.SHOrthogonalityProperty())
		expr = meh.apply_recursive(expr, meh.SummationOverKronecker())

		#meh.print_expr(expr)
		return expr

	@staticmethod
	def collision_term():
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
		L_coeffs = meh.SHCoefficient( "L", meh.var("l'"), meh.var("m'"), x )


		sigma_t = meh.fun( "\\sigma_t", x)
		expr = meh.mul(sigma_t, L_coeffs)

		#print("\n----------------------------\n")
		#print("$$\n" + meh.latex(expr) + "\n$$")
		#print(meh.hierarchy(expr))

		return expr

	@staticmethod
	def scattering_term():
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
		L_coeffs = meh.SHCoefficient( "L", meh.var("l'"), meh.var("m'"), x )

		sigma_s = meh.fun( "\\sigma_s", x)

		f_p = meh.fun( "f_p", meh.var("l'"), meh.num(0), x)
		# TODO: we skip lambda_l as this was not used in the starmap paper either
		# 		however, I think it needs to be here...

		# the negation comes from the fact that we bring this to the other side
		expr = meh.neg(meh.mul(sigma_s, f_p, L_coeffs))
		#print("\n----------------------------\n")
		#print("$$\n" + meh.latex(expr) + "\n$$")
		#print(meh.hierarchy(expr))
		return expr

	@staticmethod
	def source_term():
		x = meh.tensor("\\vec{x}", rank=1, dimension=3)
		x.setComponent(0, meh.var("x"))
		x.setComponent(1, meh.var("y"))
		x.setComponent(2, meh.var("z"))
		x.collapsed = True

		q = meh.fun( "q", meh.var("l'"), meh.var("m'"), x)
		q.setLatexArgumentPosition(0, -1)
		q.setLatexArgumentPosition(1, -1)
		expr = q
		return expr



def splitAddition( terms ):
	final_terms = []
	for term in terms:
		if term.__class__ == meh.Addition:
			#for i in range(0,5):
			#print(term.numOperands())
			for i in range(term.numOperands()):
			#for i in range(0,1):
				final_terms.append(term.getOperand(i))
		else:
			final_terms.append(term)
	return final_terms


class CountImaginaryUnits(object):
	def __init__(self):
		self.count = 0
	def visit_ImaginaryUnit(self, expr):
		self.count += 1
		return expr

class fopn_real(object):



	@staticmethod
	def transport_term_expand( m, order, debug = False, gg = None ):
		sh_real_basis = fopn_real.conj_real_shbasis(m)

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
		#L_expanded = meh.sh_expansion(L, x, omega)

		#L_expanded = meh.sh_expansion_real2(L, x, omega, meh.var('N'))

		#'''
		L_expanded = meh.sh_expansion_real(L, x, omega, meh.var('N'))
		#if debug == True:
		#	print("------------------------------")
		#	meh.print_expr(L_expanded)

		L_expanded = meh.apply_recursive(L_expanded, meh.DistributiveLaw())
		#if debug == True:
		#	print("------------------------------")
		#	meh.print_expr(L_expanded)

		L_expanded = meh.apply_recursive(L_expanded, meh.CleanupSigns())
		#if debug == True:
		#	print("------------------------------")
		#	meh.print_expr(L_expanded)

		L_expanded = meh.apply_recursive(L_expanded, meh.SplitSums())
		#if debug == True:
		#	print("------------------------------")
		#	meh.print_expr(L_expanded)

		L_expanded = meh.apply_recursive(L_expanded, meh.SplitSums())
		#if debug == True:
		#	print("------------------------------")
		#	meh.print_expr(L_expanded)

		L_expanded = meh.apply_recursive(L_expanded, meh.CleanupSigns())
		#if debug == True:
		#	print("------------------------------")
		#	meh.print_expr(L_expanded)

		L_expanded = meh.apply_recursive(L_expanded, meh.Factorize())
		#if debug == True:
		#	print("------------------------------")
		#	meh.print_expr(L_expanded)
		#'''

		#L_coeffs = meh.SHCoefficient( "L", meh.var("l"), meh.var("m"), x )

		dx_L = meh.deriv(L, x.getComponent(0), is_partial = True)
		dy_L = meh.deriv(L, x.getComponent(1), is_partial = True)
		dz_L = meh.deriv(L, x.getComponent(2), is_partial = True)
		omega_dot_nablaL = meh.add( meh.mul(omega_x, dx_L), meh.mul(omega_y, dy_L), meh.mul(omega_z, dz_L))


		expr = omega_dot_nablaL

		#if debug == True:
		#	meh.print_expr(expr)

		expr = meh.integrate(meh.mul( sh_real_basis, expr), omega) 
		
		#if debug == True:
		#	print("------------------------------")
		#	meh.print_expr(expr)

		expr = meh.apply_recursive(expr, meh.DistributiveLaw())
		#if debug == True:
		#	print("\n------------------------------")
		#	meh.print_expr(expr)

		expr = meh.apply_recursive(expr, meh.DistributiveLaw())
		#if debug == True:
		#	print("\n------------------------------")
		#	meh.print_expr(expr)

		expr = meh.apply_recursive(expr, meh.CleanupSigns())
		#if debug == True:
		#	print("\n------------------------------")
		#	meh.print_expr(expr)

		expr = meh.apply_recursive(expr, meh.SHRecursiveRelation())
		#if debug == True:
		#	print("\n------------------------------")
		#	meh.print_expr(expr)

		expr = meh.apply_recursive(expr, meh.DistributiveLaw())
		#if debug == True:
		#	print("\n------------------------------")
		#	meh.print_expr(expr)

		expr = meh.apply_recursive(expr, meh.CleanupSigns())
		#if debug == True:
		#	print("\n------------------------------")
		#	meh.print_expr(expr)

		expr = meh.apply_recursive(expr, meh.SplitIntegrals())
		#if debug == True:
		#	print("\n------------------------------")
		#	meh.print_expr(expr)

		expr = meh.apply_recursive(expr, meh.CleanupSigns())
		#if debug == True:
		#	print("\n------------------------------")
		#	meh.print_expr(expr)


		expr = meh.apply_recursive(expr, meh.Factorize())
		#if m == 0:
		#	if debug == True:
		#		print("\n------------------------------")
		#		meh.print_expr(expr)
		#gg = []
		#for i in range(expr.numOperands()):
		#	if i != 4 and i != 16:
		#		expr.setOperand(i, meh.num(0))
		#'''
		#'''

		if m < 0:
		#if False:
			def beta_m_body( m ):
				if m == -1:
					return 2.0/np.sqrt(2.0)
				else:
					return 1.0
			beta_m = meh.fun( "\\beta", meh.var("m'") )
			beta_m.setAllSuperScripts()
			beta_m.body2 = beta_m_body

			def gamma_m_body( m ):
				if m == -1:
					return 0.0
				else:
					return 1.0
			gamma_m = meh.fun( "\\gamma", meh.var("m'") )
			gamma_m.setAllSuperScripts()
			gamma_m.body2 = gamma_m_body

			# replace term
			coeff_c = meh.SHRecursiveRelation.c( meh.sub(meh.var("l'"), meh.num(1)), meh.sub( meh.var("m'"), meh.num(1) ) )
			coeff_d = meh.SHRecursiveRelation.d( meh.add(meh.var("l'"), meh.num(1)), meh.sub( meh.var("m'"), meh.num(1) ) )
			coeff_e = meh.SHRecursiveRelation.e( meh.sub(meh.var("l'"), meh.num(1)), meh.add( meh.var("m'"), meh.num(1) ) )
			coeff_f = meh.SHRecursiveRelation.f( meh.add(meh.var("l'"), meh.num(1)), meh.add( meh.var("m'"), meh.num(1) ) )
			coeff_a = meh.SHRecursiveRelation.a( meh.sub(meh.var("l'"), meh.num(1)), meh.var("m'"))
			coeff_b = meh.SHRecursiveRelation.b( meh.add(meh.var("l'"), meh.num(1)), meh.var("m'"))

			L0 = meh.SHCoefficient( 'L', meh.sub(meh.var("l'"), meh.num(1)), meh.add(meh.neg(meh.var("m'")), meh.num(1)), x )
			L1 = meh.SHCoefficient( 'L', meh.add(meh.var("l'"), meh.num(1)), meh.add(meh.neg(meh.var("m'")), meh.num(1)), x )
			L2 = meh.SHCoefficient( 'L', meh.sub(meh.var("l'"), meh.num(1)), meh.sub(meh.neg(meh.var("m'")), meh.num(1)), x )
			L3 = meh.SHCoefficient( 'L', meh.add(meh.var("l'"), meh.num(1)), meh.sub(meh.neg(meh.var("m'")), meh.num(1)), x )

			L4 = meh.SHCoefficient( 'L', meh.sub(meh.var("l'"), meh.num(1)), meh.sub(meh.var("m'"), meh.num(1)), x )
			L5 = meh.SHCoefficient( 'L', meh.sub(meh.var("l'"), meh.num(1)), meh.add(meh.var("m'"), meh.num(1)), x )
			L6 = meh.SHCoefficient( 'L', meh.add(meh.var("l'"), meh.num(1)), meh.add(meh.var("m'"), meh.num(1)), x )
			L7 = meh.SHCoefficient( 'L', meh.add(meh.var("l'"), meh.num(1)), meh.sub(meh.var("m'"), meh.num(1)), x )

			L8 = meh.SHCoefficient( 'L', meh.sub(meh.var("l'"), meh.num(1)), meh.var("m'"), x )
			L9 = meh.SHCoefficient( 'L', meh.add(meh.var("l'"), meh.num(1)), meh.var("m'"), x )

			terms = []
			# y ----------
			terms.append(meh.neg(meh.mul(meh.frac(meh.num(1), meh.num(2)), coeff_c, meh.deriv(L0, x.getComponent(1), is_partial=True))))

			terms.append(meh.mul(meh.frac(meh.num(1), meh.num(2)), coeff_d, meh.deriv(L1, x.getComponent(1), is_partial=True)))

			term2 = meh.neg(meh.mul(meh.frac(meh.num(1), meh.num(2)), coeff_e, meh.deriv(L2, x.getComponent(1), is_partial=True)))
			terms.append( meh.mul( beta_m, term2))

			term3 = meh.mul(meh.frac(meh.num(1), meh.num(2)), coeff_f, meh.deriv(L3, x.getComponent(1), is_partial=True))
			terms.append(meh.mul( beta_m, term3))

			# x ----------
			terms.append(meh.mul(meh.frac(meh.num(1), meh.num(2)), coeff_c, meh.deriv(L4, x.getComponent(0), is_partial=True)))

			term5 = meh.neg(meh.mul(meh.frac(meh.num(1), meh.num(2)), coeff_e, meh.deriv(L5, x.getComponent(0), is_partial=True)))
			terms.append(meh.mul(gamma_m, term5))

			term6 = meh.mul(meh.frac(meh.num(1), meh.num(2)), coeff_f, meh.deriv(L6, x.getComponent(0), is_partial=True))
			terms.append(meh.mul(gamma_m, term6))

			terms.append(meh.neg(meh.mul(meh.frac(meh.num(1), meh.num(2)), coeff_d, meh.deriv(L7, x.getComponent(0), is_partial=True))))

			# z ----------
			terms.append(meh.mul(coeff_a, meh.deriv(L8, x.getComponent(2), is_partial=True)))

			terms.append(meh.mul(coeff_b, meh.deriv(L9, x.getComponent(2), is_partial=True)))

			return meh.Addition(terms)
			#return meh.num(0)
		#'''
		if m == 0:
			coeff_c = meh.SHRecursiveRelation.c( meh.sub(meh.var("l'"), meh.num(1)), meh.num(-1) )
			coeff_d = meh.neg(meh.SHRecursiveRelation.d( meh.add(meh.var("l'"), meh.num(1)), meh.num(-1) ))
			coeff_a = meh.SHRecursiveRelation.a( meh.sub(meh.var("l'"), meh.num(1)), meh.num(0))
			coeff_b = meh.SHRecursiveRelation.b( meh.add(meh.var("l'"), meh.num(1)), meh.num(0))
			L0 = meh.SHCoefficient( 'L', meh.sub(meh.var("l'"), meh.num(1)), meh.num(1), x )
			L1 = meh.SHCoefficient( 'L', meh.add(meh.var("l'"), meh.num(1)), meh.num(1), x )
			L2 = meh.SHCoefficient( 'L', meh.sub(meh.var("l'"), meh.num(1)), meh.num(-1), x )
			L3 = meh.SHCoefficient( 'L', meh.add(meh.var("l'"), meh.num(1)), meh.num(-1), x )
			L4 = meh.SHCoefficient( 'L', meh.sub(meh.var("l'"), meh.num(1)), meh.num(0), x )
			L5 = meh.SHCoefficient( 'L', meh.add(meh.var("l'"), meh.num(1)), meh.num(0), x )
			c = meh.frac( meh.num(1), meh.sqrt(meh.num(2.0)) )

			terms = []
			# x ----------
			terms.append( meh.mul( c, coeff_c, meh.deriv( L0, x.getComponent(0), is_partial=True ) ) )
			terms.append( meh.mul( c, coeff_d, meh.deriv( L1, x.getComponent(0), is_partial=True ) ) )
			# y ----------
			terms.append( meh.mul( c, coeff_c, meh.deriv( L2, x.getComponent(1), is_partial=True ) ) )
			terms.append( meh.mul( c, coeff_d, meh.deriv( L3, x.getComponent(1), is_partial=True ) ) )
			# z ----------
			terms.append( meh.mul( coeff_a, meh.deriv( L4, x.getComponent(2), is_partial=True ) ) )
			terms.append( meh.mul( coeff_b, meh.deriv( L5, x.getComponent(2), is_partial=True ) ) )

			return meh.Addition(terms)
		#'''

		#'''
		if m > 0:
			#print("asfasfasf afasf")
			#pass
		#if False:
			def beta_m_body( m ):
				if m == 1:
					return 2.0/np.sqrt(2.0)
				else:
					return 1.0
			beta_m = meh.fun( "\\beta", meh.var("m'") )
			beta_m.setAllSuperScripts()
			beta_m.body2 = beta_m_body

			def gamma_m_body( m ):
				if m == 1:
					return 0.0
				else:
					return 1.0
			gamma_m = meh.fun( "\\gamma", meh.var("m'") )
			gamma_m.setAllSuperScripts()
			gamma_m.body2 = gamma_m_body

			# replace term
			coeff_c = meh.SHRecursiveRelation.c( meh.sub(meh.var("l'"), meh.num(1)), meh.sub( meh.neg(meh.var("m'")), meh.num(1) ) )
			coeff_d = meh.SHRecursiveRelation.d( meh.add(meh.var("l'"), meh.num(1)), meh.sub( meh.neg(meh.var("m'")), meh.num(1) ) )
			coeff_e = meh.SHRecursiveRelation.e( meh.sub(meh.var("l'"), meh.num(1)), meh.add( meh.neg(meh.var("m'")), meh.num(1) ) )
			coeff_f = meh.SHRecursiveRelation.f( meh.add(meh.var("l'"), meh.num(1)), meh.add( meh.neg(meh.var("m'")), meh.num(1) ) )
			coeff_a = meh.SHRecursiveRelation.a( meh.sub(meh.var("l'"), meh.num(1)), meh.neg(meh.var("m'")))
			coeff_b = meh.SHRecursiveRelation.b( meh.add(meh.var("l'"), meh.num(1)), meh.neg(meh.var("m'")))

			L0 = meh.SHCoefficient( 'L', meh.sub(meh.var("l'"), meh.num(1)), meh.add(meh.var("m'"), meh.num(1)), x )
			L1 = meh.SHCoefficient( 'L', meh.add(meh.var("l'"), meh.num(1)), meh.add(meh.var("m'"), meh.num(1)), x )
			L2 = meh.SHCoefficient( 'L', meh.sub(meh.var("l'"), meh.num(1)), meh.sub(meh.var("m'"), meh.num(1)), x )
			L3 = meh.SHCoefficient( 'L', meh.add(meh.var("l'"), meh.num(1)), meh.sub(meh.var("m'"), meh.num(1)), x )

			L4 = meh.SHCoefficient( 'L', meh.sub(meh.var("l'"), meh.num(1)), meh.sub(meh.neg(meh.var("m'")), meh.num(1)), x )
			L5 = meh.SHCoefficient( 'L', meh.add(meh.var("l'"), meh.num(1)), meh.sub(meh.neg(meh.var("m'")), meh.num(1)), x )
			L6 = meh.SHCoefficient( 'L', meh.sub(meh.var("l'"), meh.num(1)), meh.add(meh.neg(meh.var("m'")), meh.num(1)), x )
			L7 = meh.SHCoefficient( 'L', meh.add(meh.var("l'"), meh.num(1)), meh.add(meh.neg(meh.var("m'")), meh.num(1)), x )

			L8 = meh.SHCoefficient( 'L', meh.sub(meh.var("l'"), meh.num(1)), meh.var("m'"), x )
			L9 = meh.SHCoefficient( 'L', meh.add(meh.var("l'"), meh.num(1)), meh.var("m'"), x )

			terms = []

			half = meh.frac( meh.num(1), meh.num(2) )

			# x ----------
			terms.append(meh.mul(half, coeff_c, meh.deriv(L0, x.getComponent(0), is_partial=True)))
			#meh.print_expr(terms[0])

			terms.append(meh.neg(meh.mul(half, coeff_d, meh.deriv(L1, x.getComponent(0), is_partial=True))))

			term2 = meh.neg(meh.mul(half, coeff_e, meh.deriv(L2, x.getComponent(0), is_partial=True)))
			terms.append(meh.mul(beta_m, term2))

			term3 = meh.mul(half, coeff_f, meh.deriv(L3, x.getComponent(0), is_partial=True))
			terms.append(meh.mul(beta_m, term3))

			# y ----------
			terms.append(meh.mul(half, coeff_c, meh.deriv(L4, x.getComponent(1), is_partial=True)))

			terms.append(meh.neg(meh.mul(half, coeff_d, meh.deriv(L5, x.getComponent(1), is_partial=True))))

			term6 = meh.mul(half, coeff_e, meh.deriv(L6, x.getComponent(1), is_partial=True))
			terms.append( meh.mul( gamma_m, term6))

			term7 = meh.neg(meh.mul(half, coeff_f, meh.deriv(L7, x.getComponent(1), is_partial=True)))
			terms.append(meh.mul( gamma_m, term7))

			# z ----------
			terms.append(meh.mul(coeff_a, meh.deriv(L8, x.getComponent(2), is_partial=True)))

			terms.append(meh.mul(coeff_b, meh.deriv(L9, x.getComponent(2), is_partial=True)))

			tt = meh.Addition(terms)
			#return terms[7]
			#meh.print_expr(tt)
			return tt
			#return terms[0]
		#'''
		#if not gg is None:
		#	expr = meh.Addition( cancelation_terms[gg] )

		# this is for checking term 15 and 16 from the paper
		#if m > 0:
		if False:
			#pass
			#ops = [0, 12, 2, 10, 3, 11, 1, 13]
			#ops = [0, 12]
			ops = [7, 15]
			#ops = [0]
			#ops = [12]
			if len(ops) == 1:
				expr = expr.getOperand(ops[0])
			else:
				expr = meh.Addition( [expr.getOperand(op) for op in ops] )
			meh.print_expr(expr)
			#expr = meh.Addition( [expr.getOperand(4), expr.getOperand(16)] )
			#expr = meh.Addition( [expr.getOperand(5), expr.getOperand(17)] )
			#expr = meh.Addition( [expr.getOperand(6), expr.getOperand(14)] )
			#expr = meh.Addition( [expr.getOperand(7), expr.getOperand(15)] )
			#expr = meh.Addition( [expr.getOperand(0), expr.getOperand(12)] )
			#expr = meh.Addition( [expr.getOperand(2), expr.getOperand(10)] )
			#expr = meh.Addition( [expr.getOperand(3), expr.getOperand(11)] )
			#expr = meh.Addition( [expr.getOperand(1), expr.getOperand(13)] )
			#expr = meh.Addition( [expr.getOperand(8), expr.getOperand(18)] )
			#expr = meh.Addition( [expr.getOperand(9), expr.getOperand(19)] )
			#expr = expr.getOperand(0)

		#if debug == True:
		#	print("\n------------------------------")
		#	meh.print_expr(expr)

		expr = meh.apply_recursive(expr, meh.Substitute(L, L_expanded))
		#if debug == True:
		#	print("\n------------------------------")
		#	split = meh.split(expr.deep_copy())
		#	for s in split:
		#		meh.print_expr(s)

		#if debug == True:
		#	print("\n------------------------------")
		#	meh.print_expr(expr)


		#if debug == True:
		#	print("\n------------------------------")
		#	split = meh.split(expr)
		#	for s in split:
		#		meh.print_expr(s)
		expr = meh.apply_recursive(expr, meh.SplitDerivatives())

		#if debug == True:
		#	print("\n------------------------------")
		#	split = meh.split(expr)
		#	for s in split:
		#		meh.print_expr(s)

		#if debug == True:
		#	print("\n------------------------------")
		#	meh.print_expr(expr)

		expr = meh.apply_recursive(expr, meh.DistributiveLaw())

		#if debug == True:
		#	print("\n------------------------------")
		#	split = meh.split(expr)
		#	for s in split:
		#		meh.print_expr(s)


		#if debug == True:
		#	print("\n------------------------------")
		#	meh.print_expr(expr)


		expr = meh.apply_recursive(expr, meh.SplitIntegrals())
		#if debug == True:
		#	print("\n------------------------------")
		#	meh.print_expr(expr)

		expr = meh.apply_recursive(expr, meh.DistributiveLaw())

		#if debug == True:
		#	print("\n------------------------------")
		#	meh.print_expr(expr)

		expr = meh.apply_recursive(expr, meh.CleanupSigns())

		#if debug == True:
		#	print("\n------------------------------")
		#	split = meh.split(expr)
		#	for s in split:
		#		meh.print_expr(s)


		#if debug == True:
		#	print("\n------------------------------")
		#	meh.print_expr(expr)


		expr = meh.apply_recursive(expr, meh.SwitchDomains())
		expr = meh.apply_recursive(expr, meh.SwitchDomains())
		expr = meh.apply_recursive(expr, meh.SwitchDomains())


		#if debug == True:
		#	print("\n------------------------------")
		#	split = meh.split(expr)
		#	for s in split:
		#		meh.print_expr(s)
		#if debug == True:
		#	print("\n------------------------------")
		#	meh.print_expr(expr)

		expr = meh.apply_recursive(expr, meh.Factorize())

		#if debug == True:
		#	print("\n------------------------------")
		#	meh.print_expr(expr)
		#if debug == True:
		#	print("\n------------------------------")
		#	split = meh.split(expr)
		#	for s in split:
		#		meh.print_expr(s)
		#print(expr.numOperands())
		#meh.print_expr(expr.getOperand(12))
		#meh.print_expr(expr.getOperand(52))
		#meh.print_expr(expr.getOperand(8))
		#meh.print_expr(expr.getOperand(69))
		#meh.print_expr(expr.getOperand(48))
		#meh.print_expr(expr.getOperand(99))
		#meh.print_expr(expr.getOperand(18))
		#meh.print_expr(expr.getOperand(59))
		#meh.print_expr(expr.getOperand(13))
		#meh.print_expr(expr.getOperand(54))

		#if not gg is None:
		#	num = expr.numOperands()
		#	for i in range(num):
		#		if i != gg:
		#			expr.setChildren( i, meh.num(0) )

		expr = meh.apply_recursive(expr, meh.SHOrthogonalityProperty())

		#expr = expr.getOperand(0)
		#print(expr.numOperands())
		#if debug == True:
		#	print("\n------------------------------")
		#	split = meh.split(expr.deep_copy())
		#	for s in split:
		#		meh.print_expr(s)
		#if debug == True:
		#	print("\n------------------------------")
		#	split = meh.split(expr)
		#	for s in split:
		#		meh.print_expr(s)

		'''
		numOps = expr.numOperands()
		for i in range(numOps):
			op = expr.getOperand(i)	
			counter = CountImaginaryUnits()
			op = meh.apply_recursive(op, counter)
			#if counter.count % 2  == 0 :
			print("i={} #imag={}".format(i, counter.count))
		'''


		# replace N by the SH truncation order
		expr = meh.apply_recursive(expr, meh.Substitute(meh.var('N'), meh.num(order)))

		#expr = expr.getOperand(0)

		#if debug == True:
		#	print("\n------------------------------")
		#	meh.print_expr(expr)

		#ops = expr.getOperands()

		#expr.setChildren( 12, meh.num(0) )
		#expr.setChildren( 52, meh.num(0) )


		#for i in range(13,100):
		#	#if i != 12:
		#	expr.setChildren( i, meh.num(0) )
		#expr = meh.add(ops[3], ops[15], ops[21], ops[29])
		#expr = ops[2]

		#meh.print_expr(expr.getOperands()[3])





		# now expand the summations into individual terms
		expr = meh.apply_recursive(expr, meh.ExpandSums())
		expr = meh.apply_recursive(expr, meh.ExpandSums())
		

		#if debug == True:
		#	print("\n------------------------------")
		#	meh.print_expr(expr)

		# final cleanup
		expr = meh.apply_recursive(expr, meh.CleanupSigns())


		'''
		# temp ------------------------
		ll = 3
		mm = 0
		#print("l={} m={}".format(ll,mm))

		expr = meh.apply_recursive(expr, meh.Substitute(meh.var("l'"), meh.num(ll)))
		expr = meh.apply_recursive(expr, meh.Substitute(meh.var("m'"), meh.num(mm)))
		#expr = meh.apply_recursive(expr, meh.FoldConstants())
		#print(expr.numOperands())
		if debug == True:
			print("\n------------------------------")
			split = meh.split(expr)
			for s in split:
				ss = meh.apply_recursive(s, meh.FoldConstants())
				#ss = s
				meh.print_expr(ss)
		# /temp -------------------------
		'''
		
		#print(expr.numOperands())
		# split derivatives into seperate terms
		expr = meh.apply_recursive(expr, meh.SplitDerivatives())
		expr = meh.apply_recursive(expr, meh.DistributiveLaw())
		#meh.print_expr(expr)
		#print(expr.numOperands()) # 110
		#print( meh.track_term(expr, 14, meh.CleanupSigns()) )
		#for i in range(56, 57):
		#	expr.setChildren(i, meh.num(0))
		#expr.setChildren(13, meh.num(0))
		#expr.setChildren(56, meh.num(0))
		expr = meh.apply_recursive(expr, meh.CleanupSigns())
		#expr.setChildren(14, meh.num(0))
		#print(expr.numOperands()) # 120

		return expr

	def transport_term( order, debug = False ):
		equ_a = fopn_real.transport_term_expand(-1, order, debug)
		equ_b = fopn_real.transport_term_expand(0, order, debug)
		equ_c = fopn_real.transport_term_expand(1, order, debug)
		#return

		'''
		# temp ------------------------
		ll = 1
		mm = 1
		#print("l={} m={}".format(ll,mm))

		equ_c = meh.apply_recursive(equ_c, meh.Substitute(meh.var("l'"), meh.num(ll)))
		equ_c = meh.apply_recursive(equ_c, meh.Substitute(meh.var("m'"), meh.num(mm)))
		#meh.print_expr(equ_c)
		#expr = meh.apply_recursive(expr, meh.FoldConstants())
		#print(expr.numOperands())
		if debug == True:
			print("\n------------------------------")
			ss = meh.apply_recursive(equ_c, meh.FoldConstants())
			meh.print_expr(ss)
		# /temp -------------------------
		'''
		#return

		#equ_c = equ_c.getOperand(4)

		#'''
		pn_equations = []
		for l in range(0, order+1):
			for m in range(-l, l+1):

				#if not (l==2 and m==1):
				#	pn_equations.append(meh.num(0))
				#	continue

				if debug == True:
					print("l={} m={} -------------".format(l,m))

				# instantiate real PN equation for concrete l,m combination
				if m < 0:
					equ = equ_a.deep_copy()
				elif m == 0:
					equ = equ_b.deep_copy()
				elif m > 0:
					equ = equ_c.deep_copy()
				else:
					raise ValueError("unexpected")

				if debug == True:
					print("\n before ------------------------------")
					meh.print_expr(equ)

				equ = meh.apply_recursive(equ, meh.Substitute(meh.var("l'"), meh.num(l)))
				equ = meh.apply_recursive(equ, meh.Substitute(meh.var("m'"), meh.num(m)))
				equ = meh.apply_recursive(equ, meh.FoldConstants())

				# temp
				#gg = []
				#for op in equ.getOperands():
				#	gg.append(meh.apply_recursive(op, meh.FoldConstants()))
				#equ = meh.Addition(gg)
				#equ = gg[14]
				#equ = gg[38]
				#equ = gg[62]
				#equ = gg[86]

				if debug == True:
					print("\nafter ------------------------------")
					meh.print_expr(equ)

				pn_equations.append(equ)

		return pn_equations
		#'''

	def collision_term_expand(m, order, debug = False):
		#sh_real_basis_a, sh_real_basis_b, sh_real_basis_c = 
		sh_real_basis = fopn_real.conj_real_shbasis(m)

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
		#L_expanded = meh.sh_expansion(L, x, omega)
		L_expanded = meh.sh_expansion_real(L, x, omega, meh.var('N'))
		#if debug == True:
		#	print("------------------------------")
		#	meh.print_expr(L_expanded)

		L_expanded = meh.apply_recursive(L_expanded, meh.DistributiveLaw())
		#if debug == True:
		#	print("------------------------------")
		#	meh.print_expr(L_expanded)

		L_expanded = meh.apply_recursive(L_expanded, meh.CleanupSigns())
		#if debug == True:
		#	print("------------------------------")
		#	meh.print_expr(L_expanded)

		L_expanded = meh.apply_recursive(L_expanded, meh.SplitSums())
		#if debug == True:
		#	print("------------------------------")
		#	meh.print_expr(L_expanded)

		L_expanded = meh.apply_recursive(L_expanded, meh.SplitSums())
		#if debug == True:
		#	print("------------------------------")
		#	meh.print_expr(L_expanded)

		L_expanded = meh.apply_recursive(L_expanded, meh.CleanupSigns())
		#if debug == True:
		#	print("------------------------------")
		#	meh.print_expr(L_expanded)

		L_expanded = meh.apply_recursive(L_expanded, meh.Factorize())

		sigma_t = meh.fun( "\\sigma_t", x)
		# there is no negative sign, because we moved the term onto the LHS
		expr = meh.mul(sigma_t, L)

		return meh.mul( sigma_t, meh.SHCoefficient( "L", meh.var("l'"), meh.var("m'"), x ) )

		#if debug == True:
		#	meh.print_expr(expr)

		expr = meh.integrate(meh.mul( sh_real_basis, expr), omega) 

		#if debug == True:
		#	print("\n------------------------------")
		#	meh.print_expr(expr)

		expr = meh.apply_recursive(expr, meh.DistributiveLaw())

		#if debug == True:
		#	print("\n------------------------------")
		#	meh.print_expr(expr)

		expr = meh.apply_recursive(expr, meh.CleanupSigns())

		#if debug == True:
		#	print("\n------------------------------")
		#	meh.print_expr(expr)

		expr = meh.apply_recursive(expr, meh.SplitIntegrals())
		expr = meh.apply_recursive(expr, meh.CleanupSigns())

		#if debug == True:
		#	print("\n------------------------------")
		#	meh.print_expr(expr)

		expr = meh.apply_recursive(expr, meh.Factorize())
		#if debug == True:
		#	print("\n------------------------------")
		#	meh.print_expr(expr)

		expr = meh.apply_recursive(expr, meh.Substitute(L, L_expanded))
		#if debug == True:
		#	print("\n------------------------------")
		#	meh.print_expr(expr)

		expr = meh.apply_recursive(expr, meh.DistributiveLaw())
		#if debug == True:
		#	print("\n------------------------------")
		#	meh.print_expr(expr)

		expr = meh.apply_recursive(expr, meh.SplitIntegrals())
		#if debug == True:
		#	print("\n------------------------------")
		#	meh.print_expr(expr)

		expr = meh.apply_recursive(expr, meh.CleanupSigns())
		#if debug == True:
		#	print("\n------------------------------")
		#	meh.print_expr(expr)

		expr = meh.apply_recursive(expr, meh.SwitchDomains())
		expr = meh.apply_recursive(expr, meh.SwitchDomains())
		expr = meh.apply_recursive(expr, meh.Factorize())

		expr = meh.apply_recursive(expr, meh.DistributiveLaw())
		expr = meh.apply_recursive(expr, meh.CleanupSigns())
		expr = meh.apply_recursive(expr, meh.Factorize())


		expr = meh.apply_recursive(expr, meh.SHOrthogonalityProperty())

		#if debug == True:
		#	print("\n------------------------------")
		#	meh.print_expr(expr)

		# replace N by the SH truncation order
		expr = meh.apply_recursive(expr, meh.Substitute(meh.var('N'), meh.num(order)))
		#if debug == True:
		#	print("\n------------------------------")
		#	meh.print_expr(expr)

		# now expand the summations into individual terms
		expr = meh.apply_recursive(expr, meh.ExpandSums())
		expr = meh.apply_recursive(expr, meh.ExpandSums())
		#if debug == True:
		#	print("\n------------------------------")
		#	meh.print_expr(expr)

		# final cleanup
		#expr = meh.apply_recursive(expr, meh.DistributiveLaw())
		#expr = meh.apply_recursive(expr, meh.CleanupSigns())
		#if debug == True:
		#	print("\n------------------------------")
		#	meh.print_expr(expr)

		'''
		# temp ------------------------
		ll = 3
		mm = 0
		#print("l={} m={}".format(ll,mm))

		expr = meh.apply_recursive(expr, meh.Substitute(meh.var("l'"), meh.num(ll)))
		expr = meh.apply_recursive(expr, meh.Substitute(meh.var("m'"), meh.num(mm)))
		#expr = meh.apply_recursive(expr, meh.FoldConstants())
		print(expr.numOperands())
		if debug == True:
			print("\n------------------------------")
			split = meh.split(expr)
			for s in split:
				ss = meh.apply_recursive(s, meh.FoldConstants())
				#ss = s
				meh.print_expr(ss)
		# /temp -------------------------
		'''

		return expr


	def collision_term(order, debug=False):
		equ_a = fopn_real.collision_term_expand(-1, order, debug)
		equ_b = fopn_real.collision_term_expand(0, order, debug)
		equ_c = fopn_real.collision_term_expand(1, order, debug)
		#return

		pn_equations = []
		for l in range(0, order+1):
			for m in range(-l, l+1):

				#if not (l==0 and m==0):
				#	continue

				# instantiate real PN equation for concrete l,m combination
				if m < 0:
					equ = equ_a.deep_copy()
				elif m == 0:
					equ = equ_b.deep_copy()
				elif m > 0:
					equ = equ_c.deep_copy()
				else:
					#continue
					raise ValueError("unexpected")

				#if debug == True and l == 1 and m == -1:
				#	print("\n------------------------------")
				#	print("l={} m={}".format(l,m))
				#	meh.print_expr(equ)


				equ = meh.apply_recursive(equ, meh.Substitute(meh.var("l'"), meh.num(l)))
				equ = meh.apply_recursive(equ, meh.Substitute(meh.var("m'"), meh.num(m)))

				#if debug == True and l == 1 and m == -1:
				#	#equ = meh.add( equ.getOperand(0), equ.getOperand(1), equ.getOperand(6), equ.getOperand(7) )
				#	print("\n------------------------------")
				#	print("l={} m={}".format(l,m))
				#	meh.print_expr(equ)

				equ = meh.apply_recursive(equ, meh.FoldConstants())
				equ = meh.apply_recursive(equ, meh.CleanupSigns())

				#if debug == True and l == 1 and m == -1:
				#	print("\n------------------------------")
				#	print("l={} m={}".format(l,m))
				#	meh.print_expr(equ)


				#if debug == True:
				#	meh.print_expr(equ)

				pn_equations.append(equ)
		return pn_equations



	def scattering_term_expand(sh_real_basis, order, debug = False):
		omega = meh.tensor("\\omega", rank=1, dimension=3)
		omega_x = omega.getComponent(0)
		omega_y = omega.getComponent(1)
		omega_z = omega.getComponent(2)

		x = meh.tensor("\\vec{x}", rank=1, dimension=3)
		x.setComponent(0, meh.var("x"))
		x.setComponent(1, meh.var("y"))
		x.setComponent(2, meh.var("z"))
		x.collapsed = True

		#L = meh.fun( "L", x, omega)
		sigma_s = meh.fun( "\\sigma_s", x)

		lambda_l = meh.fun( "\\lambda", meh.var("l'") )
		lambda_l.setLatexArgumentPosition(0, -1)
		lambda_l.setBody( lambda l: np.sqrt(4.0*np.pi/(2*l+1)) )

		L = meh.SHCoefficient( "L", meh.var("l'"), meh.var("m'"), x )
		phase = meh.SHCoefficient( "f", meh.var("l'"), meh.num(0), x )

		# This is the real-valued scattering term...
		return meh.neg( meh.mul( sigma_s, lambda_l, phase, L ) )



		L = meh.SHCoefficient( "L", meh.var("l"), meh.var("m"), x )
		phase = meh.SHCoefficient( "f", meh.var("l"), meh.num(0), x )
		basis = meh.SHBasis(meh.var("l"), meh.var("m"), omega, conjugate_complex=False)
		term = meh.mul(L, phase, lambda_l, basis)
		sum_m = meh.sum( term, meh.var("m"), meh.neg(meh.var("l")), meh.num(-1) )
		sum_l = meh.sum( sum_m, meh.var("l"), meh.num(0), meh.var("N") )
		t0 = meh.mul(sigma_s, meh.frac(meh.imag(1), meh.sqrt(meh.num(2))), sum_l)


		L = meh.SHCoefficient( "L", meh.var("l"), meh.var("m"), x )
		phase = meh.SHCoefficient( "f", meh.var("l"), meh.num(0), x )
		basis = meh.SHBasis(meh.var("l"), meh.neg(meh.var("m")), omega, conjugate_complex=False)
		term = meh.mul(meh.pow(meh.num(-1), meh.var("m")), L, phase, lambda_l, basis)
		sum_m = meh.sum( term, meh.var("m"), meh.neg(meh.var("l")), meh.num(-1) )
		sum_l = meh.sum( sum_m, meh.var("l"), meh.num(0), meh.var("N") )
		t1 = meh.neg(meh.mul(sigma_s, meh.frac(meh.imag(1), meh.sqrt(meh.num(2))), sum_l))

		L = meh.SHCoefficient( "L", meh.var("l"), meh.num(0), x )
		phase = meh.SHCoefficient( "f", meh.var("l"), meh.num(0), x )
		basis = meh.SHBasis(meh.var("l"), meh.num(0), omega, conjugate_complex=False)
		term = meh.mul(L, phase, lambda_l, basis)
		sum_l = meh.sum( term, meh.var("l"), meh.num(0), meh.var("N") )
		t2 = meh.mul(sigma_s, sum_l)

		L = meh.SHCoefficient( "L", meh.var("l"), meh.var("m"), x )
		phase = meh.SHCoefficient( "f", meh.var("l"), meh.num(0), x )
		basis = meh.SHBasis(meh.var("l"), meh.neg(meh.var("m")), omega, conjugate_complex=False)
		term = meh.mul(L, phase, lambda_l, basis)
		sum_m = meh.sum( term, meh.var("m"), meh.num(1), meh.var("l") )
		sum_l = meh.sum( sum_m, meh.var("l"), meh.num(0), meh.var("N") )
		t3 = meh.mul(sigma_s, meh.frac(meh.num(1), meh.sqrt(meh.num(2))), sum_l)


		L = meh.SHCoefficient( "L", meh.var("l"), meh.var("m"), x )
		phase = meh.SHCoefficient( "f", meh.var("l"), meh.num(0), x )
		basis = meh.SHBasis(meh.var("l"), meh.var("m"), omega, conjugate_complex=False)
		term = meh.mul(meh.pow(meh.num(-1), meh.var("m")), L, phase, lambda_l, basis)
		sum_m = meh.sum( term, meh.var("m"), meh.num(1), meh.var("l") )
		sum_l = meh.sum( sum_m, meh.var("l"), meh.num(0), meh.var("N") )
		t4 = meh.mul(sigma_s, meh.frac(meh.num(1), meh.sqrt(meh.num(2))), sum_l)

		# we negate everythin, since the scattering term is being brought
		# from the RHS to the LHS of the RTE which results in a sign change
		expr = meh.neg(meh.add(t0, t1, t2, t3, t4))
		expr = meh.apply_recursive(expr, meh.CleanupSigns())

		#meh.print_expr(expr)

		#if debug == True:
		#	print("\n------------------------------")
		#	meh.print_expr(expr)

		# now multiply with basis and integrate over solid angle
		expr = meh.integrate(meh.mul( sh_real_basis, expr), omega) 

		#if debug == True:
		#	print("\n------------------------------")
		#	meh.print_expr(expr)

		expr = meh.apply_recursive(expr, meh.DistributiveLaw())
		expr = meh.apply_recursive(expr, meh.DistributiveLaw())
		expr = meh.apply_recursive(expr, meh.SplitIntegrals())
		expr = meh.apply_recursive(expr, meh.CleanupSigns())
		#if debug == True:
		#	print("\n------------------------------")
		#	meh.print_expr(expr)

		expr = meh.apply_recursive(expr, meh.SwitchDomains())
		expr = meh.apply_recursive(expr, meh.SwitchDomains())
		expr = meh.apply_recursive(expr, meh.Factorize())

		#if debug == True:
		#	print("\n------------------------------")
		#	meh.print_expr(expr)

		expr = meh.apply_recursive(expr, meh.SHOrthogonalityProperty())

		#expr = expr.getOperand(0)
		#if debug == True:
		#	print("\n------------------------------")
		#	meh.print_expr(expr)

		expr = meh.apply_recursive(expr, meh.Substitute(meh.var('N'), meh.num(order)))
		expr = meh.apply_recursive(expr, meh.ExpandSums())
		expr = meh.apply_recursive(expr, meh.ExpandSums())

		'''
		# temp ------------------------
		ll = 3
		mm = 0
		#print("l={} m={}".format(ll,mm))

		expr = meh.apply_recursive(expr, meh.Substitute(meh.var("l'"), meh.num(ll)))
		expr = meh.apply_recursive(expr, meh.Substitute(meh.var("m'"), meh.num(mm)))
		#expr = meh.apply_recursive(expr, meh.FoldConstants())
		print(expr.numOperands())
		if debug == True:
			print("\n------------------------------")
			split = meh.split(expr)
			for s in split:
				ss = meh.apply_recursive(s, meh.FoldConstants())
				#ss = s
				meh.print_expr(ss)
		# /temp -------------------------
		'''

		expr = meh.apply_recursive(expr, meh.DistributiveLaw())
		expr = meh.apply_recursive(expr, meh.CleanupSigns())




		#if debug == True:
		#	print("\n------------------------------")
		#	meh.print_expr(expr)

		#if debug == True:
		#	print("\n------------------------------")
		#	meh.print_expr(expr)

		'''
		expr = meh.apply_recursive(expr, meh.SwitchDomains())
		expr = meh.apply_recursive(expr, meh.Factorize())
		expr = meh.apply_recursive(expr, meh.SHOrthogonalityProperty())
		expr = meh.apply_recursive(expr, meh.Substitute(meh.var('N'), meh.num(order)))
		expr = meh.apply_recursive(expr, meh.ExpandSums())
		expr = meh.apply_recursive(expr, meh.DistributiveLaw())
		expr = meh.apply_recursive(expr, meh.CleanupSigns())
		'''


		return expr


	def conj_real_shbasis( m = None):
		omega = meh.tensor("\\omega", rank=1, dimension=3)
		omega_x = omega.getComponent(0)
		omega_y = omega.getComponent(1)
		omega_z = omega.getComponent(2)

		f1 = meh.mul(meh.frac(meh.imag(-1), meh.sqrt(meh.num(2))), meh.SHBasis(meh.var("l'"), meh.var("m'"), omega, conjugate_complex=True))
		f2 = meh.mul(meh.frac(meh.imag(-1), meh.sqrt(meh.num(2))), meh.pow(meh.num(-1), meh.var("m'")), meh.SHBasis(meh.var("l'"), meh.neg(meh.var("m'")), omega, conjugate_complex=True))
		sh_real_basis_a = meh.sub(f1, f2)
		sh_real_basis_b = meh.SHBasis(meh.var("l'"), meh.var("m'"), omega, conjugate_complex=True)
		f1 = meh.mul(meh.frac(meh.num(1), meh.sqrt(meh.num(2))), meh.SHBasis(meh.var("l'"), meh.neg(meh.var("m'")), omega, conjugate_complex=True))
		f2 = meh.mul(meh.frac(meh.num(1), meh.sqrt(meh.num(2))), meh.pow(meh.num(-1), meh.var("m'")), meh.SHBasis(meh.var("l'"), meh.var("m'"), omega, conjugate_complex=True))
		sh_real_basis_c = meh.add(f1, f2)

		if m is None:
			return sh_real_basis_a, sh_real_basis_b, sh_real_basis_c
		elif m<0:
			return sh_real_basis_a
		elif m==0:
			return sh_real_basis_b
		elif m>0:
			return sh_real_basis_c



	def scattering_term(order, debug=False):

		sh_real_basis_a, sh_real_basis_b, sh_real_basis_c = fopn_real.conj_real_shbasis()

		#if debug == True:
		#	print("\n------------------------------")
		#	meh.print_expr(sh_real_basis_c)


		equ_a = fopn_real.scattering_term_expand(sh_real_basis_a, order, debug)
		equ_b = fopn_real.scattering_term_expand(sh_real_basis_b, order, debug)
		equ_c = fopn_real.scattering_term_expand(sh_real_basis_c, order, debug)
		#return
		

		#'''
		pn_equations = []
		for l in range(0, order+1):
			for m in range(-l, l+1):

				#if not (l==0 and m==0):
				#	continue
				#if debug == True:
				#	print("\n------------------------------")
				#	print("l={} m={}".format(l,m))

				# instantiate real PN equation for concrete l,m combination
				if m < 0:
					equ = equ_a.deep_copy()
				elif m == 0:
					equ = equ_b.deep_copy()
				elif m > 0:
					equ = equ_c.deep_copy()
				else:
					#continue
					raise ValueError("unexpected")


				equ = meh.apply_recursive(equ, meh.Substitute(meh.var("l'"), meh.num(l)))
				equ = meh.apply_recursive(equ, meh.Substitute(meh.var("m'"), meh.num(m)))
				equ = meh.apply_recursive(equ, meh.FoldConstants())
				equ = meh.apply_recursive(equ, meh.CleanupSigns())

				#if debug == True:
				#	meh.print_expr(equ)

				pn_equations.append(equ)
		return pn_equations
		#'''

	def source_term_expand(sh_real_basis, order, debug = False):
		omega = meh.tensor("\\omega", rank=1, dimension=3)
		omega_x = omega.getComponent(0)
		omega_y = omega.getComponent(1)
		omega_z = omega.getComponent(2)

		x = meh.tensor("\\vec{x}", rank=1, dimension=3)
		x.setComponent(0, meh.var("x"))
		x.setComponent(1, meh.var("y"))
		x.setComponent(2, meh.var("z"))
		x.collapsed = True

		return meh.SHCoefficient( "Q", meh.var("l'"), meh.var("m'"), x )

		Q = meh.fun( "Q", x, omega)
		#L_expanded = meh.sh_expansion(L, x, omega)
		Q_expanded = meh.sh_expansion_real(Q, x, omega, meh.var('N'))
		#if debug == True:
		#	print("------------------------------")
		#	meh.print_expr(L_expanded)

		Q_expanded = meh.apply_recursive(Q_expanded, meh.DistributiveLaw())
		#if debug == True:
		#	print("------------------------------")
		#	meh.print_expr(L_expanded)

		Q_expanded = meh.apply_recursive(Q_expanded, meh.CleanupSigns())
		#if debug == True:
		#	print("------------------------------")
		#	meh.print_expr(L_expanded)

		Q_expanded = meh.apply_recursive(Q_expanded, meh.SplitSums())
		#if debug == True:
		#	print("------------------------------")
		#	meh.print_expr(L_expanded)

		Q_expanded = meh.apply_recursive(Q_expanded, meh.SplitSums())
		#if debug == True:
		#	print("------------------------------")
		#	meh.print_expr(L_expanded)

		Q_expanded = meh.apply_recursive(Q_expanded, meh.CleanupSigns())
		#if debug == True:
		#	print("------------------------------")
		#	meh.print_expr(L_expanded)

		Q_expanded = meh.apply_recursive(Q_expanded, meh.Factorize())
		#if debug == True:
		#	print("------------------------------")
		#	meh.print_expr(Q_expanded)

		expr = Q

		#if debug == True:
		#	meh.print_expr(expr)

		expr = meh.integrate(meh.mul( sh_real_basis, expr), omega) 

		#if debug == True:
		#	print("------------------------------")
		#	meh.print_expr(expr)

		expr = meh.apply_recursive(expr, meh.DistributiveLaw())

		#if debug == True:
		#	print("------------------------------")
		#	meh.print_expr(expr)

		expr = meh.apply_recursive(expr, meh.CleanupSigns())

		#if debug == True:
		#	print("------------------------------")
		#	meh.print_expr(expr)

		expr = meh.apply_recursive(expr, meh.SplitIntegrals())
		expr = meh.apply_recursive(expr, meh.CleanupSigns())

		#if debug == True:
		#	print("------------------------------")
		#	meh.print_expr(expr)

		expr = meh.apply_recursive(expr, meh.Factorize())
		#if debug == True:
		#	print("------------------------------")
		#	meh.print_expr(expr)

		expr = meh.apply_recursive(expr, meh.Substitute(Q, Q_expanded))
		#if debug == True:
		#	print("------------------------------")
		#	meh.print_expr(expr)

		expr = meh.apply_recursive(expr, meh.DistributiveLaw())
		#if debug == True:
		#	print("------------------------------")
		#	meh.print_expr(expr)

		expr = meh.apply_recursive(expr, meh.SplitIntegrals())
		#if debug == True:
		#	print("------------------------------")
		#	meh.print_expr(expr)

		expr = meh.apply_recursive(expr, meh.CleanupSigns())
		#if debug == True:
		#	print("------------------------------")
		#	meh.print_expr(expr)

		expr = meh.apply_recursive(expr, meh.SwitchDomains())
		expr = meh.apply_recursive(expr, meh.SwitchDomains())
		expr = meh.apply_recursive(expr, meh.Factorize())
		#if debug == True:
		#	print("------------------------------")
		#	meh.print_expr(expr)

		expr = meh.apply_recursive(expr, meh.SHOrthogonalityProperty())

		expr = meh.apply_recursive(expr, meh.DistributiveLaw())
		expr = meh.apply_recursive(expr, meh.CleanupSigns())
		expr = meh.apply_recursive(expr, meh.Factorize())
		#if debug == True:
		#	print("------------------------------")
		#	meh.print_expr(expr)		

		# replace N by the SH truncation order
		expr = meh.apply_recursive(expr, meh.Substitute(meh.var('N'), meh.num(order)))
		#if debug == True:
		#	print("------------------------------")
		#	meh.print_expr(expr)

		# now expand the summations into individual terms
		expr = meh.apply_recursive(expr, meh.ExpandSums())
		expr = meh.apply_recursive(expr, meh.ExpandSums())
		#if debug == True:
		#	print("------------------------------")
		#	meh.print_expr(expr)


		'''
		# temp ------------------------
		ll = 3
		mm = -2
		#print("l={} m={}".format(ll,mm))

		expr = meh.apply_recursive(expr, meh.Substitute(meh.var("l'"), meh.num(ll)))
		expr = meh.apply_recursive(expr, meh.Substitute(meh.var("m'"), meh.num(mm)))
		#expr = meh.apply_recursive(expr, meh.FoldConstants())
		print(expr.numOperands())
		if debug == True:
			print("\n------------------------------")
			split = meh.split(expr)
			for s in split:
				ss = meh.apply_recursive(s, meh.FoldConstants())
				#ss = s
				meh.print_expr(ss)
		# /temp -------------------------
		'''

		# final cleanup
		expr = meh.apply_recursive(expr, meh.DistributiveLaw())
		expr = meh.apply_recursive(expr, meh.CleanupSigns())
		#if debug == True:
		#	print("------------------------------")
		#	meh.print_expr(expr)

		return expr



	def source_term(order, debug=False):

		sh_real_basis_a, sh_real_basis_b, sh_real_basis_c = fopn_real.conj_real_shbasis()

		#if debug == True:
		#	print("\n------------------------------")
		#	meh.print_expr(sh_real_basis_c)


		equ_a = fopn_real.source_term_expand(sh_real_basis_a, order, debug)
		equ_b = fopn_real.source_term_expand(sh_real_basis_b, order, debug)
		equ_c = fopn_real.source_term_expand(sh_real_basis_c, order, debug)
		#return
		

		#meh.print_expr(equ_a)
		#return equ_a

		#'''
		pn_equations = []
		for l in range(0, order+1):
			for m in range(-l, l+1):

				#if not (l==0 and m==0):
				#	continue
				if debug == True:
					print("\n------------------------------")
					print("l={} m={}".format(l,m))

				# instantiate real PN equation for concrete l,m combination
				if m < 0:
					equ = equ_a.deep_copy()
				elif m == 0:
					equ = equ_b.deep_copy()
				elif m > 0:
					equ = equ_c.deep_copy()
				else:
					#continue
					raise ValueError("unexpected")


				equ = meh.apply_recursive(equ, meh.Substitute(meh.var("l'"), meh.num(l)))
				equ = meh.apply_recursive(equ, meh.Substitute(meh.var("m'"), meh.num(m)))
				equ = meh.apply_recursive(equ, meh.FoldConstants())
				equ = meh.apply_recursive(equ, meh.CleanupSigns())

				if debug == True:
					meh.print_expr(equ)

				pn_equations.append(equ)
		return pn_equations
		#'''

if __name__ == "__main__":
	#pass
	#test = lspn.term0_projected_expr(True)
	order = 2
	test = fopn_real.transport_term( order, True)
	#test = fopn_real.collision_term( order, True)
	#test = fopn_real.scattering_term( order, True)
	#test = fopn_real.source_term( order, True)

	#print("test ------------")
	#for p in test:
	#	meh.print_expr(p)

	#l = 1
	#m = 0

	'''
	# find term origins
	for i in range(100):

		equ = fopn_real.transport_term_expand(-1, order, True, gg = i)


		equ = meh.apply_recursive(equ, meh.Substitute(meh.var("l'"), meh.num(l)))
		equ = meh.apply_recursive(equ, meh.Substitute(meh.var("m'"), meh.num(m)))
		equ = meh.apply_recursive(equ, meh.FoldConstants())

		print("\n i={} ------------------------------".format(i))
		meh.print_expr(equ)
	exit(1)
	'''



	#for i in range(10):
	#i = 4
	#print("i={} -----------------------".format(i))
	#equ = fopn_real.transport_term_expand2(-1, order, True)
	
	#equ = fopn_real.transport_term_expand(-1, order, True)
	#equ = meh.apply_recursive(equ, meh.Substitute(meh.var("l'"), meh.num(l)))
	#equ = meh.apply_recursive(equ, meh.Substitute(meh.var("m'"), meh.num(m)))
	#equ = meh.apply_recursive(equ, meh.FoldConstants())
	#meh.print_expr(equ)


	'''
	equ = meh.apply_recursive(equ, meh.Substitute(meh.var("l'"), meh.num(l)))
	equ = meh.apply_recursive(equ, meh.Substitute(meh.var("m'"), meh.num(m)))
	equ = meh.apply_recursive(equ, meh.FoldConstants())


	# temp
	#gg = []
	#for op in equ.getOperands():
	#	gg.append(meh.apply_recursive(op, meh.FoldConstants()))
	#equ = meh.Addition(gg)
	#equ = gg[14]
	#equ = gg[38]
	#equ = gg[62]
	#equ = gg[86]
	#print("\n------------------------------")
	meh.print_expr(equ)
	'''