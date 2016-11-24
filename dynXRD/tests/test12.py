# import os
# os.chdir(os.path.expanduser("/home/federica/Documenti/gitrep/dynXRD/"))
#import pyasf
import sympy as sp

import pylab as pl

from dynXRD import reflectivity




#import matplotlib.pyplot as pl
#import xraylib
from dynXRD.tests.crystals import MgO

data = pl.loadtxt("test12.dat")
#data[:,0] = pl.radians(data[:,0])
#pl.ion()

R = 0,0,2
thickness=1000
Energy=10000

#cryst1 = reflectivity.crystal("cif/MgO1000053.cif") #MgO
struct1= MgO.structure_MgO()
cryst1 = auxfunc.crystal(struct1)
#cryst2 = reflectivity.crystal("cif/MgO1000053.cif") #MgO
struct2= MgO.structure_MgO()
cryst2 = auxfunc.crystal(struct2)

# cryst1 = reflectivity.crystal("cif/Si2104737.cif")
# cryst2 = reflectivity.crystal("cif/Si2104737.cif")
# cryst1 = reflectivity.crystal("Si")
# cryst2 = reflectivity.crystal("Si")
cryst2.lattice_par_val[cryst2.a] *= 1.01 #strained layer
Sub= reflectivity.Substrate(cryst1)
v_par=sp.Matrix([1,0,0])
v_perp=sp.Matrix([0,0,1])
Sub.calc_orientation(v_par, v_perp)

layer1= reflectivity.Epitaxial_Layer(cryst2, thickness)
layer1.calc_orientation(v_par, v_perp)
crystal= reflectivity.Sample(Sub, layer1)
crystal.set_Miller(R)
crystal.calc_g0_gH(Energy)
thBragg= float(Sub.calc_Bragg_angle(Energy).subs(Sub.structure.lattice_par_val).evalf())
angle=pl.linspace(0.975, 1.006,501)*thBragg

# XRl = layer1.calc_reflection_amplitude(angle, Energy)
#XRs = Sub.calc_reflection_amplitude(angle, Energy)
#
# XT = layer1.calc_transmission_amplitude(angle, Energy)

XR=crystal.calc_reflectivity(angle, Energy)

crystal.print_values(angle, Energy)

pl.plot(*data.T, label='GID_sl', color='red')
#pl.plot(angle-thBragg,abs(XT)**2)
pl.plot(pl.degrees(angle-thBragg),abs(XR)**2, label='dynXRD', color='black')
#pl.plot(pl.degrees(angle-thBragg),abs(Sub.XR)**2, label='dynXRD', color='blue')
#pl.plot(angle-thBragg,1 - abs(XT)**2 - abs(XRl)**2)
pl.yscale('log')
pl.xlabel('Angle (degrees)')
pl.ylabel('Reflectivity')
pl.xlim(-0.25,0.1)
pl.rc('font', size=18)
pl.legend(loc="upper left", prop={'size':19})

#pl.savefig('pics/test12.eps')
pl.show()
