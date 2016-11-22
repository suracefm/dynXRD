# import os
# os.chdir(os.path.expanduser("/home/federica/Documenti/gitrep/dynXRD/"))
#import xraylib
#import pyasf
import reflectivity
import sympy as sp
import matplotlib.pyplot as pl
import numpy as np
#import aux_pyasf as auxfunc
#import aux_xraylib as auxfunc

#data = pl.loadtxt("MgO_100nm_002.dat")
#data[:,0] = pl.radians(data[:,0])
##pl.ion()

R = 1,1,1
thickness=1000
Energy=10000

#cryst = reflectivity.crystal("cif/Ge7101739.cif")
cryst = reflectivity.crystal("Ge")
Sub=reflectivity.Substrate(cryst)
v_par=sp.Matrix([1,0,0])
v_perp=sp.Matrix([0,0,1])
Sub.calc_orientation(v_par, v_perp)
layer1=reflectivity.Epitaxial_Layer(cryst, thickness)
layer1.calc_orientation(v_par, v_perp)
crystal=reflectivity.Sample(Sub, layer1)
crystal.set_Miller(R)
crystal.calc_g0_gH(Energy)

angle=np.linspace(0.98, 1.02,501)*float(layer1.calc_Bragg_angle(Energy).subs(layer1.structure.lattice_par_val).evalf())

crystal.calc_reflectivity(angle, Energy)
layer1.calc_amplitudes(angle, Energy)
Sub.calc_amplitudes(angle, Energy)

XRl = layer1.XR
XRs = Sub.XR
XT = layer1.XT

#pl.plot(*data.T)
#pl.plot(angle,abs(XT)**2)

pl.plot(angle,abs(XRl)**2)
#pl.plot(angle,1 - abs(XT)**2 - abs(XRl)**2)
pl.show()

#print(Sub.structure.M, Sub.xrlstructure.M_matrix)
