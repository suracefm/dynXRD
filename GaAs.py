import os
os.chdir(os.path.expanduser("~/Desktop/reflectivity/dynXRD/"))
import pyasf
import reflectivity
import sympy as sp
import pylab as pl
import numpy as np


R = 0,0,4
thickness=np.array([3, 1.5, 1, 0.5])*1e4
strain=np.array([20, 15, 10, 5])*1e-4
Energy=8048


structsub = pyasf.unit_cell("9008845") #GaAs
Sub=reflectivity.Substrate(structsub)

Layers=[]
for i in range(4):
    layer=reflectivity.Epitaxial_Layer(pyasf.unit_cell("9008845"), thickness[i])
    Layers.append(layer)
    Layers[i].structure.subs[Layers[i].structure.a]*=(1+strain[i])


psi=0
v_perp=sp.Matrix([0,0,1])
Sub.calc_orientation_from_angle(psi, v_perp)


for i in range(4):
    Layers[i].calc_orientation_from_angle(psi, v_perp)

crystal=reflectivity.Sample(Sub, *Layers)
crystal.set_Miller(R)
#crystal.calc_layer_Miller()
crystal.calc_g0_gH(Energy)
#thBragg= float(Sub.calc_Bragg_angle(Energy).subs(Sub.structure.subs).evalf())
thBragg=Sub.thetaBragg
angle=pl.linspace(0.9975, 1.0008,501)*thBragg

#XRl = layer1.calc_reflection_amplitude(angle, Energy)
#XRs = Sub.calc_reflection_amplitude(angle, Energy)
# XT = layer1.calc_transmission_amplitude(angle, Energy)
XR=crystal.calc_reflectivity(angle, Energy)
crystal.print_values(angle, Energy)

#pl.plot(data[:,0], data[:,1])
# pl.plot(angle-thBragg,abs(XT)**2)
# pl.plot(angle-thBragg,abs(XRl)**2)
# pl.plot(angle-thBragg,1 - abs(XT)**2 - abs(XRl)**2)


pl.plot((angle-thBragg)*3600*180/np.pi,100*abs(XR)**2, color='black')
# pl.plot((angle-thBragg)*3600*180/np.pi,abs(Sub.XR)**2)
# pl.plot((angle-thBragg)*3600*180/np.pi,abs(Layers[0].XT)**2)
# pl.plot((angle-thBragg)*3600*180/np.pi,abs(Layers[1].XT)**2)
# pl.plot((angle-thBragg)*3600*180/np.pi,abs(Layers[2].XT)**2)
# pl.plot((angle-thBragg)*3600*180/np.pi,abs(Layers[3].XT)**2)


#pl.plot(angle-thBragg,abs(XRs)**2)
#pl.plot(angle-thBragg,abs(XRl)**2)
#pl.yscale('log')
pl.ylim(0,100)
pl.xlabel('Angle (sec)')
pl.ylabel('Reflectivity (%)')
pl.rc('font', size=17)

pl.savefig('pics/GaAs.eps')
pl.show()