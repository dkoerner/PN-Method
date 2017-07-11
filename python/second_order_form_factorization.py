import cas




if __name__ == "__main__":

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

	#
	# second order transport term -------------------------------------------
	#
	'''
	Ylm = cas.SHBasis(cas.var("l'"), cas.var("m'"), omega, conjugate_complex=True)
	dx_L = cas.deriv(L, x.getComponent(0), is_partial = True)
	dy_L = cas.deriv(L, x.getComponent(1), is_partial = True)
	dz_L = cas.deriv(L, x.getComponent(2), is_partial = True)
	omega__dot_nablaL = cas.add( cas.mul(omega_x, dx_L), cas.mul(omega_y, dy_L), cas.mul(omega_z, dz_L))
	expr_x = cas.neg(cas.deriv(cas.integrate( cas.mul( omega_x, Ylm, omega__dot_nablaL), omega ), x.getComponent(0), is_partial=True))
	expr_y = cas.neg(cas.deriv(cas.integrate( cas.mul( omega_y, Ylm, omega__dot_nablaL), omega ), x.getComponent(1), is_partial=True))
	expr_z = cas.neg(cas.deriv(cas.integrate( cas.mul( omega_z, Ylm, omega__dot_nablaL), omega ), x.getComponent(2), is_partial=True))
	expr = cas.add(expr_x, expr_y, expr_z)

	#print("\n----------------------------\n")
	#print("$$\n" + cas.latex(expr) + "\n$$")
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
	#print("\n----------------------------\n")
	#print("$$\n" + cas.latex(expr) + "\n$$")
	#print(cas.hierarchy(expr))
	'''


	#
	# extinction directional derivative term -----------------------------
	#
	'''
	omegaL = cas.tensor("", rank=1, dimension=3)
	omegaL.setComponent(0, cas.mul(omega_x, L))
	omegaL.setComponent(1, cas.mul(omega_y, L))
	omegaL.setComponent(2, cas.mul(omega_z, L))

	nabla_sigma_t = cas.tensor("", rank=1, dimension=3)
	nabla_sigma_t.setComponent(0, cas.deriv(sigma_t, x.getComponent(0), is_partial = True))
	nabla_sigma_t.setComponent(1, cas.deriv(sigma_t, x.getComponent(1), is_partial = True))
	nabla_sigma_t.setComponent(2, cas.deriv(sigma_t, x.getComponent(2), is_partial = True))

	expr = cas.neg(cas.dot(omegaL, nabla_sigma_t))

	print("\n----------------------------\n")
	print("$$\n" + cas.latex(expr) + "\n$$")
	expr = cas.apply_recursive(expr, cas.ExpandDotProduct())
	expr = cas.integrate(cas.mul( cas.SHBasis(cas.var("l'"), cas.var("m'"), omega, conjugate_complex=True), expr), omega) 
	expr = cas.apply_recursive(expr, cas.CleanupSigns())
	expr = cas.apply_recursive(expr, cas.DistributiveLaw())
	expr = cas.apply_recursive(expr, cas.CleanupSigns())

	#print("\n----------------------------\n")
	#print("expanding the dot product, projecting into SH basis function and integrating over solid angle:\n")
	#print("$$\n" + cas.latex(expr) + "\n$$")

	expr = cas.apply_recursive(expr, cas.SplitIntegrals())
	expr = cas.apply_recursive(expr, cas.CleanupSigns())
	expr = cas.apply_recursive(expr, cas.SHRecursiveRelation())
	expr = cas.apply_recursive(expr, cas.DistributiveLaw())
	expr = cas.apply_recursive(expr, cas.CleanupSigns())
	expr = cas.apply_recursive(expr, cas.SplitIntegrals())
	expr = cas.apply_recursive(expr, cas.CleanupSigns())
	expr = cas.apply_recursive(expr, cas.Factorize())
	expr = cas.apply_recursive(expr, cas.Substitute(L, L_expanded))
	expr = cas.apply_recursive(expr, cas.SwitchDomains())
	expr = cas.apply_recursive(expr, cas.SwitchDomains())
	expr = cas.apply_recursive(expr, cas.Factorize())
	expr = cas.apply_recursive(expr, cas.SHOrthogonalityProperty())
	expr = cas.apply_recursive(expr, cas.SummationOverKronecker())

	print("\n----------------------------\n")
	print("Applying recursive relations, replacing $L$ by its SH expansion and rearranging terms:\n")
	print("$$\n" + cas.latex(expr) + "\n$$")
	'''

	#
	# Squared extinction term -----------------------------
	#
	# trivial....


	#
	# Directional derivative scattering term -----------------------------
	#
	'''
	SL_isotropic_expanded = cas.sum( cas.sum( cas.mul( cas.fun( "\\lambda", cas.var("l"), arglevel=-1), cas.SHCoefficient( "f_p", cas.var("l"), cas.num(0), x ), cas.SHCoefficient( "L", cas.var("l"), cas.var("m"), x ), cas.SHBasis(cas.var("l"), cas.var("m"), omega, conjugate_complex=False) ), cas.var('m'), cas.neg(cas.var('l')), cas.var('l') ), cas.var('l'), cas.num(0), cas.infty() )

	nabla_SL = cas.tensor("", rank=1, dimension=3)
	nabla_SL.setComponent(0, cas.deriv(SL_isotropic_expanded, x.getComponent(0), is_partial = True))
	nabla_SL.setComponent(1, cas.deriv(SL_isotropic_expanded, x.getComponent(1), is_partial = True))
	nabla_SL.setComponent(2, cas.deriv(SL_isotropic_expanded, x.getComponent(2), is_partial = True))


	expr = cas.neg( cas.dot(omega, nabla_SL) )
	# project into SH and integrate over solid angle
	expr = cas.integrate(cas.mul( cas.SHBasis(cas.var("l'"), cas.var("m'"), omega, conjugate_complex=True), expr), omega) 

	expr = cas.apply_recursive(expr, cas.CleanupSigns())
	expr = cas.apply_recursive(expr, cas.ExpandDotProduct())
	expr = cas.apply_recursive(expr, cas.DistributiveLaw())
	expr = cas.apply_recursive(expr, cas.SplitIntegrals())
	expr = cas.apply_recursive(expr, cas.CleanupSigns())
	expr = cas.apply_recursive(expr, cas.SwitchDomains())
	expr = cas.apply_recursive(expr, cas.SwitchDomains())
	expr = cas.apply_recursive(expr, cas.SwitchDomains())
	expr = cas.apply_recursive(expr, cas.Factorize())
	expr = cas.apply_recursive(expr, cas.SHRecursiveRelation())
	expr = cas.apply_recursive(expr, cas.DistributiveLaw())
	expr = cas.apply_recursive(expr, cas.SplitIntegrals())
	expr = cas.apply_recursive(expr, cas.CleanupSigns())
	expr = cas.apply_recursive(expr, cas.Factorize())
	expr = cas.apply_recursive(expr, cas.DistributiveLaw())
	expr = cas.apply_recursive(expr, cas.SHOrthogonalityProperty())
	expr = cas.apply_recursive(expr, cas.DistributiveLaw())
	expr = cas.apply_recursive(expr, cas.CleanupSigns())
	expr = cas.apply_recursive(expr, cas.SplitSums())
	expr = cas.apply_recursive(expr, cas.CleanupSigns())
	expr = cas.apply_recursive(expr, cas.SummationOverKronecker())
	expr = cas.apply_recursive(expr, cas.SplitDerivatives())
	expr = cas.apply_recursive(expr, cas.CleanupSigns())
	expr = cas.apply_recursive(expr, cas.Factorize())

	print("\n----------------------------\n")
	print("$$\n" + cas.latex(expr) + "\n$$")
	'''


	#
	# Extinction scattering term -----------------------------
	#
	# kind of simple too


	#
	# Directional derivative source term -----------------------------
	#
	'''
	nabla_Q = cas.tensor("", rank=1, dimension=3)
	nabla_Q.setComponent(0, cas.deriv(Q, x.getComponent(0), is_partial = True))
	nabla_Q.setComponent(1, cas.deriv(Q, x.getComponent(1), is_partial = True))
	nabla_Q.setComponent(2, cas.deriv(Q, x.getComponent(2), is_partial = True))

	expr = cas.neg(cas.dot(omega, nabla_Q))
	print("\n----------------------------\n")
	print("$$\n" + cas.latex(expr) + "\n$$")

	expr = cas.apply_recursive(expr, cas.ExpandDotProduct())
	expr = cas.apply_recursive(expr, cas.CleanupSigns())
	# project into SH and integrate over solid angle
	expr = cas.integrate(cas.mul( cas.SHBasis(cas.var("l'"), cas.var("m'"), omega, conjugate_complex=True), expr), omega) 
	expr = cas.apply_recursive(expr, cas.DistributiveLaw())
	expr = cas.apply_recursive(expr, cas.CleanupSigns())
	expr = cas.apply_recursive(expr, cas.SplitIntegrals())
	expr = cas.apply_recursive(expr, cas.CleanupSigns())
	expr = cas.apply_recursive(expr, cas.SHRecursiveRelation())
	expr = cas.apply_recursive(expr, cas.DistributiveLaw())
	expr = cas.apply_recursive(expr, cas.SplitIntegrals())
	expr = cas.apply_recursive(expr, cas.CleanupSigns())
	expr = cas.apply_recursive(expr, cas.Factorize())
	expr = cas.apply_recursive(expr, cas.Substitute(Q, Q_expanded))
	expr = cas.apply_recursive(expr, cas.SwitchDomains())
	expr = cas.apply_recursive(expr, cas.SwitchDomains())
	expr = cas.apply_recursive(expr, cas.SwitchDomains())
	expr = cas.apply_recursive(expr, cas.Factorize())
	expr = cas.apply_recursive(expr, cas.SHOrthogonalityProperty())
	expr = cas.apply_recursive(expr, cas.SummationOverKronecker())

	print("\n----------------------------\n")
	print("$$\n" + cas.latex(expr) + "\n$$")
	'''

	
	#
	# Extinction source term -----------------------------
	#
	'''
	expr = cas.mul(sigma_t, Q)
	print("\n----------------------------\n")
	print("$$\n" + cas.latex(expr) + "\n$$")

	expr = cas.integrate(cas.mul( cas.SHBasis(cas.var("l'"), cas.var("m'"), omega, conjugate_complex=True), expr), omega) 
	expr = cas.apply_recursive(expr, cas.Substitute(Q, Q_expanded))
	expr = cas.apply_recursive(expr, cas.SwitchDomains())
	expr = cas.apply_recursive(expr, cas.SwitchDomains())
	expr = cas.apply_recursive(expr, cas.Factorize())
	expr = cas.apply_recursive(expr, cas.SHOrthogonalityProperty())
	expr = cas.apply_recursive(expr, cas.SummationOverKronecker())

	print("\n----------------------------\n")
	print("$$\n" + cas.latex(expr) + "\n$$")
	'''
