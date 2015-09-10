import os
#os.chdir(os.path.expanduser("~/Desktop/reflectivity/dynXRD/"))
import pyasf
import reflectivity
import sympy as sp
import pylab as pl
import numpy as np


# Substrate
structsub = pyasf.unit_cell("9008845") #GaAs
Sub=reflectivity.Substrate(structsub)

# Strained layers
struct1=pyasf.unit_cell("9008845") #GaAs
tvector=np.array([0, 3, 4.5, 5.5, 6])*1e4
strain=np.array([20, 15, 10, 5])*1e-4
layer1=reflectivity.Strained_Layer(struct1, tvector, a=strain)

# Geometry of substrate and layers
psi=0
v_perp=sp.Matrix([0,0,1])
Sub.calc_orientation_from_angle(psi, v_perp)
layer1.calc_orientation_from_angle(psi, v_perp)

# Crystal sample    
crystal=reflectivity.Sample(Sub, layer1)

# Miller indices
R = 0,0,4
crystal.set_Miller(R)

# Reflectivity
Energy=8048
crystal.calc_g0_gH(Energy)
thBragg= Sub.thetaBragg
angle=pl.linspace(0.9975, 1.0008,501)*thBragg
XR=crystal.calc_reflectivity(angle, Energy)

# Values and Plot
crystal.print_values(angle, Energy)
pl.plot((angle-thBragg), abs(XR)**2)

#pl.plot((angle-thBragg)*3600*180/np.pi,abs(XR)**2)
# pl.plot((angle-thBragg)*3600*180/np.pi,abs(Sub.XR)**2)
# pl.plot((angle-thBragg)*3600*180/np.pi,abs(layer1.XTvec[0])**2)
# pl.plot((angle-thBragg)*3600*180/np.pi,abs(layer1.XTvec[1])**2)
# pl.plot((angle-thBragg)*3600*180/np.pi,abs(layer1.XTvec[2])**2)
# pl.plot((angle-thBragg)*3600*180/np.pi,abs(layer1.XTvec[3])**2)

#pl.plot(angle-thBragg,abs(XRs)**2)
#pl.plot(angle-thBragg,abs(XRl)**2)
#pl.yscale('log')

pl.show()