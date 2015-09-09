import os
os.chdir(os.path.expanduser("~/Desktop/reflectivity/dynXRD/"))
import pyasf
import reflectivity
import sympy as sp
import pylab as pl

data = pl.loadtxt("test12.dat")
#data[:,0] = pl.radians(data[:,0])
#pl.ion()

R = 0,0,2
thickness=1000
Energy=10000

#struct = pyasf.unit_cell("1011169")
struct1 = pyasf.unit_cell("1011117") #MgO
struct2 = pyasf.unit_cell("1011117")
struct2.subs[struct2.a] *= 1.01 #strained layer
Sub=reflectivity.Substrate(struct1)
v_par=sp.Matrix([1,0,0])
v_perp=sp.Matrix([0,0,1])
Sub.calc_orientation(v_par, v_perp)

layer1=reflectivity.Epitaxial_Layer(struct2, thickness)
layer1.calc_orientation(v_par, v_perp)
crystal=reflectivity.Sample(Sub, layer1)
crystal.set_Miller(R)
crystal.calc_g0_gH(Energy)
thBragg= float(Sub.calc_Bragg_angle(Energy).subs(Sub.structure.subs).evalf())
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

pl.savefig('pics/test12.eps')
pl.show()