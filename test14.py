import os
os.chdir(os.path.expanduser("~/Desktop/reflectivity/"))
import pyasf
import reflectivity
import sympy as sp
import pylab as pl

data = pl.loadtxt("/afs/desy.de/user/s/suracefm/Desktop/reflectivity/test14.dat")
data[:,0] = pl.radians(data[:,0])
#pl.ion()

R = 0,0,2
thickness=1000
Energy=10000

#struct = pyasf.unit_cell("1011117") 
struct1 = pyasf.unit_cell("1507756") #BaTiO3
struct2 = pyasf.unit_cell("1507756")
struct2.subs[struct2.c] *= 1.01 #strained layer
Sub=reflectivity.Substrate(struct1)
#v_par=sp.Matrix([3,-2,-1])
psi=sp.pi/8
v_perp=sp.Matrix([1,1,3])
Sub.calc_orientation_from_angle(psi, v_perp)
#Sub.set_Miller(R)
layer1=reflectivity.Epitaxial_Layer(struct2, thickness)
layer1.calc_orientation_from_angle(psi, v_perp)
crystal=reflectivity.Sample(Sub, layer1)
#crystal.calc_layer_Miller()
crystal.set_Miller(R)
crystal.calc_g0_gH(Energy)
thBragg= float(Sub.calc_Bragg_angle(Energy).subs(Sub.structure.subs).evalf())
angle=pl.linspace(0.993, 1.007,501)*thBragg

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