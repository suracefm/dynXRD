import os
os.chdir(os.path.expanduser("~/Desktop/reflectivity/"))
import pyasf
import reflectivity
import sympy as sp
import pylab as pl
import numpy as np


R = 0,0,4
#thickness=np.array([3, 1.5, 1, 0.5])*1e4
strain=np.array([20, 15, 10, 5])*1e-4
Energy=8048


structsub = pyasf.unit_cell("9008845") #GaAs
Sub=reflectivity.Substrate(structsub)
struct1=pyasf.unit_cell("9008845")
tvector=np.array([0, 3, 4.5, 5.5, 6])*1e4
layer1=reflectivity.Strained_Layer(struct1, tvector, a=strain)

# Layers=[]
# for i in range(4):
#     layer=reflectivity.Epitaxial_Layer(pyasf.unit_cell("9008845"), thickness[i])
#     Layers.append(layer)
#     Layers[i].structure.subs[Layers[i].structure.a]*=(1+strain[i])


psi=0
v_perp=sp.Matrix([0,0,1])
Sub.calc_orientation_from_angle(psi, v_perp)
layer1.calc_orientation_from_angle(psi, v_perp)


# for i in range(4):
#     Layers[i].calc_orientation_from_angle(psi, v_perp)

crystal=reflectivity.Sample(Sub, layer1)
crystal.set_Miller(R)
#crystal.calc_layer_Miller()
crystal.calc_g0_gH(Energy)
thBragg= float(Sub.calc_Bragg_angle(Energy).subs(Sub.structure.subs).evalf())
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

pl.plot((angle-thBragg)*3600*180/np.pi,abs(XR)**2)
# pl.plot((angle-thBragg)*3600*180/np.pi,abs(Sub.XR)**2)
# pl.plot((angle-thBragg)*3600*180/np.pi,abs(layer1.XTvec[0])**2)
# pl.plot((angle-thBragg)*3600*180/np.pi,abs(layer1.XTvec[1])**2)
# pl.plot((angle-thBragg)*3600*180/np.pi,abs(layer1.XTvec[2])**2)
# pl.plot((angle-thBragg)*3600*180/np.pi,abs(layer1.XTvec[3])**2)

#pl.plot(angle-thBragg,abs(XRs)**2)
#pl.plot(angle-thBragg,abs(XRl)**2)
#pl.yscale('log')

pl.show()