import os
os.chdir(os.path.expanduser("~/Desktop/reflectivity/"))
import pyasf
import reflectivity
import sympy as sp
import pylab as pl

data = pl.loadtxt("/afs/desy.de/user/s/suracefm/Desktop/reflectivity/test11.dat")
data[:,0] = pl.radians(data[:,0])
#pl.ion()

R = 3,0,0
thickness=1000
Energy=10000

#struct = pyasf.unit_cell("1521772")
struct = pyasf.unit_cell("/afs/desy.de/user/s/suracefm/Desktop/reflectivity/cif/LiNbO3_28294.cif") #Li Nb O3
Sub=reflectivity.Substrate(struct)
v_par=sp.Matrix([-1,1,4])#1,0,2
v_perp=sp.Matrix([1,1,0]) #1,2,0
Sub.calc_orientation(v_par, v_perp)

layer1=reflectivity.Epitaxial_Layer(struct, thickness)
layer1.calc_orientation(v_par, v_perp)
crystal=reflectivity.Sample(Sub, layer1)
crystal.set_Miller(R)
#crystal.calc_layer_Miller()
crystal.calc_g0_gH(Energy)
thBragg= float(layer1.calc_Bragg_angle(Energy).subs(layer1.structure.subs).evalf())
angle=pl.linspace(0.995, 1.005,501)*thBragg

XRl = layer1.calc_reflection_amplitude(angle, Energy)
XRs = Sub.calc_reflection_amplitude(angle, Energy)
XT = layer1.calc_transmission_amplitude(angle, Energy)

crystal.print_values(angle, Energy)

pl.plot(data[:,0], data[:,1])
pl.plot(angle-thBragg,abs(XT)**2)
pl.plot(angle-thBragg,abs(XRl)**2)
pl.plot(angle-thBragg,1 - abs(XT)**2 - abs(XRl)**2)
pl.yscale('log')
pl.show()