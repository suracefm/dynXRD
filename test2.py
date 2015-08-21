import os
os.chdir(os.path.expanduser("~/Desktop/reflectivity/"))
import pyasf
import reflectivity
import sympy as sp
import pylab as pl

data = pl.loadtxt("/afs/desy.de/user/s/suracefm/Desktop/reflectivity/test2.dat")
data[:,0] = pl.radians(data[:,0])
#pl.ion()

R = 0,0,2
thickness=1000
Energy=10000

#struct = pyasf.unit_cell("1011169")
struct = pyasf.unit_cell("1011117") #MgO
Sub=reflectivity.Substrate(struct)
v_par=sp.Matrix([1,-1,0])
v_perp=sp.Matrix([1,1,1])
Sub.calc_orientation(v_par, v_perp)
Sub.set_Miller(R)
layer1=reflectivity.Epitaxial_Layer(struct, thickness)
layer1.calc_orientation(v_par, v_perp)
crystal=reflectivity.Sample(Sub, layer1)
crystal.calc_layer_Miller()
crystal.calc_g0_gH(Energy)
thBragg= float(layer1.calc_Bragg_angle(Energy).subs(layer1.structure.subs).evalf())
angle=pl.linspace(0.993, 1.007,501)*thBragg

XRl = layer1.calc_reflection_amplitude(angle, Energy)
XRs = Sub.calc_reflection_amplitude(angle, Energy)
XT = layer1.calc_transmission_amplitude(angle, Energy)

crystal.print_values(angle, Energy)

pl.plot(data[:,0], data[:,1])
pl.plot(angle-thBragg,abs(XT)**2-1)
pl.plot(angle-thBragg,abs(XRl)**2)
pl.plot(angle-thBragg,1 - abs(XT)**2 - abs(XRl)**2)
pl.yscale('log')
pl.show()