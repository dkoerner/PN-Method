import cas




if __name__ == "__main__":

	#
	# extinction derivative term -------------------------------------------
	#
	omega = cas.tensor("\\omega", rank=1, dimension=3)
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

	omegaL = cas.tensor("", rank=1, dimension=3)
	omegaL.setComponent(0, cas.mul(omega.getComponent(0), L))
	omegaL.setComponent(1, cas.mul(omega.getComponent(1), L))
	omegaL.setComponent(2, cas.mul(omega.getComponent(2), L))

	nablasigma_t = cas.tensor("", rank=1, dimension=3)
	nablasigma_t.setComponent(0, cas.mul(cas.nabla().getComponent(0), sigma_t))
	nablasigma_t.setComponent(1, cas.mul(cas.nabla().getComponent(1), sigma_t))
	nablasigma_t.setComponent(2, cas.mul(cas.nabla().getComponent(2), sigma_t))


	# extinction directional derivative term -----------------------------
	'''
	expr = cas.neg(cas.dot(omegaL, nablasigma_t))

	print("\n----------------------------\n")
	print("$$\n" + cas.latex(expr) + "\n$$")

	expr = cas.apply_recursive(expr, cas.ExpandDotProduct())
	expr = cas.integrate(cas.mul( cas.SHBasis(cas.var("l'"), cas.var("m'"), omega, conjugate_complex=True), expr), omega) 
	expr = cas.apply_recursive(expr, cas.CleanupSigns())
	expr = cas.apply_recursive(expr, cas.DistributiveLaw())
	expr = cas.apply_recursive(expr, cas.CleanupSigns())

	print("\n----------------------------\n")
	print("expanding the dot product, projecting into SH basis function and integrating over solid angle:\n")
	print("$$\n" + cas.latex(expr) + "\n$$")

	expr = cas.apply_recursive(expr, cas.SplitIntegrals())
	expr = cas.apply_recursive(expr, cas.CleanupSigns())
	expr = cas.apply_recursive(expr, cas.SHRecursiveRelation())
	expr = cas.apply_recursive(expr, cas.DistributiveLaw())
	expr = cas.apply_recursive(expr, cas.CleanupSigns())
	expr = cas.apply_recursive(expr, cas.SplitIntegrals())
	expr = cas.apply_recursive(expr, cas.CleanupSigns())
	expr = cas.apply_recursive(expr, cas.Substitute(L, L_expanded))
	expr = cas.apply_recursive(expr, cas.Factorize())
	expr = cas.apply_recursive(expr, cas.SwitchDomains())
	expr = cas.apply_recursive(expr, cas.SwitchDomains())
	expr = cas.apply_recursive(expr, cas.Factorize())
	expr = cas.apply_recursive(expr, cas.SHOrthogonalityProperty())
	expr = cas.apply_recursive(expr, cas.SummationOverKronecker())
	
	print("\n----------------------------\n")
	print("Applying recursive relations, replacing $L$ by its SH expansion and rearranging terms:\n")
	print("$$\n" + cas.latex(expr) + "\n$$")
	'''

	# source directional derivative term -----------------------------
	'''
	expr = cas.neg(cas.mul(cas.dot(omega, cas.nabla()), Q))
	print("\n----------------------------\n")
	print("$$\n" + cas.latex(expr) + "\n$$")

	expr = cas.apply_recursive(expr, cas.ExpandDotProduct())
	expr = cas.apply_recursive(expr, cas.DistributiveLaw())
	expr = cas.apply_recursive(expr, cas.CleanupSigns())
	expr = cas.integrate(cas.mul( cas.SHBasis(cas.var("l'"), cas.var("m'"), omega, conjugate_complex=True), expr), omega) 
	expr = cas.apply_recursive(expr, cas.DistributiveLaw())
	expr = cas.apply_recursive(expr, cas.CleanupSigns())
	expr = cas.apply_recursive(expr, cas.SplitIntegrals())
	expr = cas.apply_recursive(expr, cas.CleanupSigns())
	expr = cas.apply_recursive(expr, cas.SHRecursiveRelation())
	expr = cas.apply_recursive(expr, cas.DistributiveLaw())
	expr = cas.apply_recursive(expr, cas.SplitIntegrals())
	expr = cas.apply_recursive(expr, cas.CleanupSigns())
	expr = cas.apply_recursive(expr, cas.Substitute(Q, Q_expanded))
	expr = cas.apply_recursive(expr, cas.SwitchDomains())
	expr = cas.apply_recursive(expr, cas.SwitchDomains())
	expr = cas.apply_recursive(expr, cas.Factorize())
	expr = cas.apply_recursive(expr, cas.SHOrthogonalityProperty())
	expr = cas.apply_recursive(expr, cas.SummationOverKronecker())


	print("\n----------------------------\n")
	print("$$\n" + cas.latex(expr) + "\n$$")
	'''
	


	# directional derivative scattering term ---------------
	SL_isotropic_expanded = cas.sum( cas.sum( cas.mul( cas.fun( "\\lambda", cas.var("l"), arglevel=-1), cas.SHCoefficient( "f_p", cas.var("l"), cas.num(0), x ), cas.SHCoefficient( "L", cas.var("l"), cas.var("m"), x ), cas.SHBasis(cas.var("l"), cas.var("m"), omega, conjugate_complex=False) ), cas.var('m'), cas.neg(cas.var('l')), cas.var('l') ), cas.var('l'), cas.num(0), cas.infty() )

	expr = cas.neg( cas.dot(omega, cas.nabla_new(SL_isotropic_expanded)) )
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
	#print("$$\n" + cas.hierarchy(expr) + "\n$$")

	'''

	print("\n----------------------------\n")
	print("$$\n" + cas.latex(expr) + "\n$$")


	expr = cas.apply_recursive(expr, cas.DistributiveLaw())

	print("\n----------------------------\n")
	print("$$\n" + cas.latex(expr) + "\n$$")
	'''

	# extinction scattering term -----------------------------



	# directional derivative source term --------------------------

	# extinction source term -----------------------------
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
