import os
os.chdir(os.path.expanduser("~/Desktop/reflectivity/dynXRD/"))
import pyasf
import reflectivity
import sympy as sp
import pylab as pl

data = pl.loadtxt("MgO_100nm_002.dat")
data[:,0] = pl.radians(data[:,0])
#pl.ion()

R = 0,0,2
thickness=1000
Energy=10000

#struct = pyasf.unit_cell("1011169")
struct = pyasf.unit_cell("1011117") #MgO
Sub=reflectivity.Substrate(struct)
v_par=sp.Matrix([1,0,0])
v_perp=sp.Matrix([0,0,1])
Sub.calc_orientation(v_par, v_perp)
Sub.set_Miller(R)
layer1=reflectivity.Epitaxial_Layer(struct, thickness)
layer1.calc_orientation(v_par, v_perp)
crystal=reflectivity.Sample(Sub, layer1)
crystal.calc_layer_Miller()
crystal.calc_g0_gH(Energy)

angle=pl.linspace(0.98, 1.02,501)*float(layer1.calc_Bragg_angle(Energy).subs(layer1.structure.subs).evalf())

XRl = layer1.calc_reflection_amplitude(angle, Energy)
XRs = Sub.calc_reflection_amplitude(angle, Energy)

XT = layer1.calc_transmission_amplitude(angle, Energy)

pl.plot(*data.T)
pl.plot(angle,abs(XT)**2)
pl.plot(angle,abs(XRl)**2)
pl.plot(angle,1 - abs(XT)**2 - abs(XRl)**2)
pl.show()