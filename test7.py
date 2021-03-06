import os
os.chdir(os.path.expanduser("~/Desktop/reflectivity/dynXRD/"))
import pyasf
import reflectivity
import sympy as sp
import pylab as pl

data = pl.loadtxt("test7.dat")
#data[:,0] = pl.radians(data[:,0])
#pl.ion()

R = 0,0,2
thickness=1000
Energy=10000

#struct = pyasf.unit_cell("1011117") 
struct = pyasf.unit_cell("cif/MgO_52026.cif") #MgO
Sub=reflectivity.Substrate(struct)
v_par=sp.Matrix([3,-2,-1])
v_perp=sp.Matrix([1,1,1])
Sub.calc_orientation(v_par, v_perp)
layer1=reflectivity.Epitaxial_Layer(struct, thickness)
layer1.calc_orientation(v_par, v_perp)
crystal=reflectivity.Sample(Sub, layer1)
crystal.set_Miller(R)
crystal.calc_g0_gH(Energy)
thBragg= float(layer1.calc_Bragg_angle(Energy).subs(layer1.structure.subs).evalf())
angle=pl.linspace(0.988, 1.012,501)*thBragg

crystal.calc_reflectivity(angle, Energy)
layer1.calc_amplitudes(angle, Energy)
Sub.calc_amplitudes(angle, Energy)

XRl = layer1.XR
XRs = Sub.XR
XT = layer1.XT


crystal.print_values(angle, Energy)

pl.plot(data[:,0], data[:,1], label='GID_sl', color='red')
# pl.plot(angle-thBragg,abs(XT)**2-1)
pl.plot(pl.degrees(angle-thBragg),abs(XRl)**2, label='dynXRD', color='black')
# pl.plot(angle-thBragg,1 - abs(XT)**2 - abs(XRl)**2)

pl.xlabel('Angle (degrees)')
pl.ylabel('Reflectivity')
pl.xlim(-0.1,0.1)
pl.yscale('log')
pl.rc('font', size=18)
pl.legend(loc="upper left", prop={'size':19})

pl.savefig('pics/test7.eps')
pl.show()