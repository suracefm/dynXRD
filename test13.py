# import os
# os.chdir(os.path.expanduser("~/Desktop/reflectivity/dynXRD/"))
# import pyasf
import reflectivity
import sympy as sp
import pylab as pl
import aux_xraylib as auxfunc
import MgO

data = pl.loadtxt("test13.dat")
data[:,0] = pl.radians(data[:,0])
#pl.ion()

R = 0,0,2
thickness=1000
Energy=10000

#struct = pyasf.unit_cell("1011117")
#struct1 = pyasf.unit_cell("1011117") #MgO
struct1=MgO.structure_MgO()
cryst1 = auxfunc.crystal(struct1)
#cryst2 = reflectivity.crystal("cif/MgO1000053.cif") #MgO
struct2=MgO.structure_MgO()
cryst2 = auxfunc.crystal(struct2)
cryst2.lattice_par_val[cryst2.a] *= 1.01 #strained layer
Sub=reflectivity.Substrate(cryst1)
v_par=sp.Matrix([3,-2,-1])
v_perp=sp.Matrix([1,1,1])
Sub.calc_orientation(v_par, v_perp)
#Sub.set_Miller(R)
layer1=reflectivity.Epitaxial_Layer(cryst2, thickness)
layer1.calc_orientation(v_par, v_perp)
crystal=reflectivity.Sample(Sub, layer1)
crystal.set_Miller(R)
#crystal.calc_layer_Miller()
crystal.calc_g0_gH(Energy)
thBragg= float(Sub.calc_Bragg_angle(Energy).subs(Sub.structure.lattice_par_val).evalf())
angle=pl.linspace(0.99, 1.01,501)*thBragg

#XRl = layer1.calc_reflection_amplitude(angle, Energy)
#XRs = Sub.calc_reflection_amplitude(angle, Energy)
# XT = layer1.calc_transmission_amplitude(angle, Energy)
XR=crystal.calc_reflectivity(angle, Energy)
crystal.print_values(angle, Energy)

pl.plot(data[:,0], data[:,1])
# pl.plot(angle-thBragg,abs(XT)**2)
# pl.plot(angle-thBragg,abs(XRl)**2)
# pl.plot(angle-thBragg,1 - abs(XT)**2 - abs(XRl)**2)
pl.plot(angle-thBragg,abs(XR)**2)
#pl.plot(angle-thBragg,abs(XRs)**2)
#pl.plot(angle-thBragg,abs(XRl)**2)
pl.yscale('log')
pl.show()
