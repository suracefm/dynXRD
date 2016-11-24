# import os
# os.chdir(os.path.expanduser("~/Desktop/reflectivity/dynXRD/"))
# import pyasf
import sympy as sp

import pylab as pl

from dynXRD import reflectivity
import aux_xraylib as auxfunc
from dynXRD.tests.crystals import MgO


data = pl.loadtxt("test1.dat")
data[:,0] = pl.radians(data[:,0])
#pl.ion()

R = 0,0,2
thickness=1000
Energy=10000

#struct = pyasf.unit_cell("1011169")
# struct = pyasf.unit_cell("1011117") #MgO
struct= MgO.structure_MgO()
cryst = auxfunc.crystal(struct)
Sub= reflectivity.Substrate(cryst)
v_par=sp.Matrix([1,0,0])
v_perp=sp.Matrix([0,0,1])
Sub.calc_orientation(v_par, v_perp)
layer1= reflectivity.Epitaxial_Layer(cryst, thickness)
layer1.calc_orientation(v_par, v_perp)
crystal= reflectivity.Sample(Sub, layer1)
crystal.set_Miller(R)
crystal.calc_g0_gH(Energy)

angle=pl.linspace(0.98, 1.02,501)*float(layer1.calc_Bragg_angle(Energy).subs(layer1.structure.lattice_par_val).evalf())

crystal.calc_reflectivity(angle, Energy)
layer1.calc_amplitudes(angle, Energy)
Sub.calc_amplitudes(angle, Energy)

XRl = layer1.XR
XRs = Sub.XR
XT = layer1.XT

pl.plot(*data.T)
#pl.plot(angle,abs(XT)**2)
pl.plot(angle,abs(XRl)**2)
#pl.plot(angle,1 - abs(XT)**2 - abs(XRl)**2)
pl.show()