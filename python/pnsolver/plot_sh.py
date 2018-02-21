import numpy as np
import matplotlib.pyplot as plt
import renderer

#renderer.sh_basis( 0, 0, 0.5, 0.0 )
#renderer.test()
#exit(1)

def sh_delta( order, theta, phi ):
	delta_theta = 0.5*np.pi
	delta_phi = 0.0
	sum = 0.0
	for l in range(0,order+1):
		for m in range(-l,l+1):
			g = renderer.sh_basis_conj( l,m, delta_theta, delta_phi )
			sum += renderer.sh_basis( l,m, theta, phi )*g
	return sum



renderer.sh_init()

print( renderer.sh_basis( 0, 0, 0.0, 0.0 ) )
print( np.sqrt(1.0/(4.0*np.pi)) )
exit(1)

order = 1

theta = np.arange(0, np.pi, 0.01)

ax = plt.subplot(111, projection='polar')

r = np.zeros(theta.shape)
for i in range(theta.shape[0]):
	r[i] = sh_delta( order, theta[i], 0.0 )
ax.plot(theta, r)
#ax.set_rmax(2)
#ax.set_rticks([0.5, 1, 1.5, 2])  # less radial ticks
#ax.set_rlabel_position(-22.5)  # get radial labels away from plotted line
#ax.grid(True)

#ax.set_title("A line plot on a polar axis", va='bottom')
plt.show()