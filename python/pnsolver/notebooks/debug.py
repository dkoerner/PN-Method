import os,sys,inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir) 

import util
import stencil
import pnsolver
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import scipy.io






#filename = "C:/projects/epfl/epfl17/python/pnsolver/results/studies/checkerboard_sopn_p1_not13.mat"
filename = "C:/projects/epfl/epfl17/python/notebooks/sopn/solution_checkerboard_blur10.0_term1.mat"
data = util.load_pn_system(filename)
staggered = False
pni = stencil.PNInfo2D( data["order"], staggered )
domain = pnsolver.Domain( np.array([7.0, 7.0]), data["resolution"], np.array([0.0, 0.0]) )

x = data["x"]
x_complex = pni.to_complex(x)



debug_working_filename = "C:/projects/epfl/epfl17/python/debug_working.mat"
debug_working_data = scipy.io.loadmat(debug_working_filename)
# for some reason exporting a shape=(144700,) results in shape=(1, 144700)
debug_working_x_real = debug_working_data["x_real"].T 
debug_working_x_complex = debug_working_data["x_complex"].T

#util.compare_matrices(debug_working_x_real, x, "working", "ours")
util.compare_matrices(debug_working_x_complex, x_complex, "working", "ours")

print(debug_working_x_complex[6406, 0])
print(x_complex[6406, 0])



'''
# contruct radiance field from coefficients
L = pnsolver.SHEXP( data["order"] )
for index in range(data["numCoeffs"]):
    (l,m) = pni.lm_index(index)
    sh_index = util.sh_index(l,m)
    offset = pni.getOffset(index)*0.5
    
    u0 = util.extract_coefficient_field( x_complex, data["resolution"], data["numCoeffs"], index )
    L.setCoefficientField( l, m, pnsolver.VoxelGrid( u0, domain, offset ) )
    #print(u0.dtype)
    

#pWS = np.array([3.5, 2.8])
pWS = np.array([2.8, 3.5])
L.eval2(pWS)

'''