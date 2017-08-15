# This module contains the terms for the standard first order form of the RTE
# These terms are discretized automatically by the PNBuilder into Ax=b form.

import meh





def transport_term():

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


	expr = meh.integrate(meh.mul( meh.SHBasis(meh.var("l'"), meh.var("m'"), omega, conjugate_complex=True), expr), omega) 
	expr = meh.apply_recursive(expr, meh.DistributiveLaw())
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

	print_expr(expr)
	return expr


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

def source_term():
	x = meh.tensor("\\vec{x}", rank=1, dimension=3)
	x.setComponent(0, meh.var("x"))
	x.setComponent(1, meh.var("y"))
	x.setComponent(2, meh.var("z"))
	x.collapsed = True

	q = meh.fun( "q", meh.var("l'"), meh.var("m'"), x)
	expr = q
	return expr








if __name__ == "__main__":
	pass