import os
os.chdir(os.path.expanduser("~/Desktop/reflectivity/dynXRD/"))
import pyasf
import reflectivity
import sympy as sp
import pylab as pl
import numpy as np


data = pl.loadtxt("sto_m28_ez_8_myt1_norm.dat")

R = 0,0,2
xmax=25
strainmax=0.0106
t=.62*1e4
x=np.linspace(0, xmax, 101)
strain=(strainmax*np.exp(-x)*(1+x))[::-1]
strain=np.delete(strain, 0)
tvector=(x*t)
#pl.plot(tvector[1:], strain[::-1])
#pl.show()

Energy=8000


structsub = pyasf.unit_cell("1507756") #BaTiO3
structsub.subs[structsub.a] = 3.9064
structsub.subs[structsub.c] = 3.9064
Sub=reflectivity.Substrate(structsub)
struct1=pyasf.unit_cell("1507756")
struct1.subs[struct1.a] = 3.9064
struct1.subs[struct1.c] = 3.9064
layer1=reflectivity.Strained_Layer(struct1, tvector, c=strain)



psi=0
v_perp=sp.Matrix([0,0,1])
Sub.calc_orientation_from_angle(psi, v_perp)
layer1.calc_orientation_from_angle(psi, v_perp)

crystal=reflectivity.Sample(Sub, layer1)
crystal.set_Miller(R)
crystal.calc_g0_gH(Energy)
thBragg= float(Sub.calc_Bragg_angle(Energy).subs(Sub.structure.subs).evalf())
angle=pl.linspace(data[:,0][0]*np.pi/180, data[:,0][-1]*np.pi/180,len(data[:,0]))


XR=crystal.calc_reflectivity(angle, Energy)
crystal.print_values(angle, Energy)



#pl.plot(angle*180/np.pi,abs(layer1.XR)**2)

#pl.plot(angle*180/np.pi,abs(Sub.XR)**2)

pl.plot(data[:,0],abs(XR)**2)

pl.plot(data[:,0], data[:,1]*9.5*1e-4)

pl.yscale('log')

pl.show()