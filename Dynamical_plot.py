import os
import pyasf
import reflectivity
import sympy as sp
import pylab as pl

#data = pl.loadtxt("test2.dat")
#data[:,0] = pl.radians(data[:,0])
#pl.ion()

R = 0,0,2
thickness=8000
Energy=10000

#struct = pyasf.unit_cell("1011169")
struct = pyasf.unit_cell("1011117") #MgO
Sub=reflectivity.Substrate(struct)
v_par=sp.Matrix([1,0,0])
v_perp=sp.Matrix([0,0,1])
Sub.calc_orientation(v_par, v_perp)
layer1=reflectivity.Epitaxial_Layer(struct, thickness)
layer1.calc_orientation(v_par, v_perp)
layer2=reflectivity.Epitaxial_Layer(struct, 15000)
layer2.calc_orientation(v_par, v_perp)
layer3=reflectivity.Epitaxial_Layer(struct, 30000)
layer3.calc_orientation(v_par, v_perp)
crystal=reflectivity.Sample(Sub, layer1, layer2, layer3)
crystal.set_Miller(R)
crystal.calc_g0_gH(Energy)
thBragg= float(layer1.calc_Bragg_angle(Energy).subs(layer1.structure.subs).evalf())
angle=pl.linspace(0.9995, 1.0005,501)*thBragg

crystal.calc_reflectivity(angle, Energy)
layer1.calc_amplitudes(angle, Energy)
layer2.calc_amplitudes(angle, Energy)
layer3.calc_amplitudes(angle, Energy)
Sub.calc_amplitudes(angle, Energy)

XR1 = layer1.XR
XR2 = layer2.XR
XR3 = layer3.XR
XRs = Sub.XR
#XT = layer1.XT

crystal.print_values(angle, Energy)

#pl.plot(data[:,0], data[:,1], label='GID_sl', color='red')
# pl.plot(angle-thBragg,abs(XT)**2-1)
pl.plot(pl.degrees(angle-thBragg),abs(XRs)**2, label='$\infty$', color='black')
pl.plot(pl.degrees(angle-thBragg),abs(XR1)**2, label='$500$ nm', color='red')
pl.plot(pl.degrees(angle-thBragg),abs(XR2)**2, label='$1000$ nm', color='blue')
pl.plot(pl.degrees(angle-thBragg),abs(XR3)**2, label='$2000$ nm', color='green')
# pl.plot(angle-thBragg,1 - abs(XT)**2 - abs(XRl)**2)

pl.xlabel('Angle (degrees)')
pl.ylabel('Reflectivity')
pl.xlim(-0.005,0.008)
#pl.yscale('log')
pl.legend(loc="upper left", title="Thickness", prop={'size':11})

pl.savefig('pics/dynamical.eps')
pl.show()