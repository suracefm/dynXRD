import os
os.chdir(os.path.expanduser("~/Desktop/reflectivity/dynXRD"))
import pyasf
import reflectivity
import sympy as sp
import pylab as pl
import numpy as np
import lmfit
import rexs


strainmax=0.011423522689800361
straindepth=10333.260937547546
scale=1.9996547379290535
thetaoffset=-0.068374122229889256
m=0.00037979104292197959




data = rexs.tools.loaddat("sto_m28_ez_112-404.dat")
feff = pl.loadtxt("s82_f1f2.dat")

R = 0,0,2

def my_strain_profile(strainmax, straindepth, xmax=25, npoints=501):
    x=np.linspace(0, xmax, 501)
    strain=(strainmax*np.exp(-x)*(1+x))[::-1]
    strain=np.delete(strain, 0)
    return x*straindepth, strain




#pl.plot(tvector[1:], strain[::-1])
#pl.show()

fields=data.fields[1::12]
Energy=np.array([float(en) for en in fields])
Q=data["Q"]


structsub = pyasf.unit_cell("1512124") #BaTiO3
structsub.get_tensor_symmetry()
structsub.build_unit_cell()
structsub.feed_feff("ti", *feff.T)

structsub.subs[structsub.a] = 3.9064
Sub=reflectivity.Substrate(structsub)

struct1=pyasf.unit_cell("1525437")

struct1.add_atom("Sr", [0,0,0])
struct1.AU_positions.pop("Ba1")
struct1.AU_formfactors.pop("Ba1")
struct1.U.pop("Ba1")
struct1.get_tensor_symmetry()
struct1.build_unit_cell()
struct1.feed_feff("Ti1", *feff.T)

struct1.subs[struct1.a] = 3.9064
struct1.subs[struct1.c] = 3.9064


tvector, strain = my_strain_profile(strainmax, straindepth)
layer1=reflectivity.Strained_Layer(struct1, tvector, c=strain)



psi=0
v_perp=sp.Matrix([0,0,1])
Sub.calc_orientation_from_angle(psi, v_perp)
layer1.calc_orientation_from_angle(psi, v_perp)

crystal=reflectivity.Sample(Sub, layer1)
crystal.set_Miller(R)

XR_matrix = []
angle_matrix=[]

for en in Energy:
    crystal.calc_g0_gH(en)
    #thBragg= float(Sub.calc_Bragg_angle(en).subs(Sub.structure.subs).evalf())
    angle=np.arcsin(Q*12398/(4*np.pi*en))
    angle_matrix.append(angle)
    XR=crystal.calc_reflectivity(angle+thetaoffset*np.pi/180, en)
    #crystal.print_values(angle, en)
    XR_matrix.append(XR)
    pl.plot(angle*180/np.pi+thetaoffset, abs(XR)**2)
    pl.plot(angle*180/np.pi+thetaoffset, (scale-m*en)*data[str(en)])
    pl.yscale('log')
    pl.show()

XR_matrix = np.array(XR_matrix)
angle_matrix = np.array(angle_matrix)




#pl.plot(angle*180/np.pi,abs(layer1.XR)**2)

#pl.plot(angle*180/np.pi,abs(Sub.XR)**2)

#pl.plot(data[:,0],abs(XR)**2)

#pl.plot(data[:,0]+thetaoffset, data[:,1]*scale)

#pl.yscale('log')

#pl.show()

# 
# 
# def residual(params, Energy, data, Q):
#     strainmax = params['strainmax'].value
#     straindepth = params['straindepth'].value
#     scale = params['scale'].value
#     m = params['m'].value
#     thetaoffset = params['thetaoffset'].value
#     
#     
#     
#     c = layer1.strain.keys()[0]
#     tvector, strain = my_strain_profile(strainmax, straindepth)
#     layer1.strain[c] = strain
#     layer1.thickness=tvector
#     res=[]
#     for en in Energy:
#         crystal.calc_g0_gH(en)
#         #thBragg= float(Sub.calc_Bragg_angle(en).subs(Sub.structure.subs).evalf())
#         angle=np.arcsin(Q*12398/(4*np.pi*en))
#         XR=crystal.calc_reflectivity(angle+thetaoffset*np.pi/180, en)
#         err=(scale-m*en)*data[str(en)]-abs(XR)**2
#         res.append(err)
#     res=np.concatenate(res)
#     print (res**2).sum()
#     return res
# 
params = lmfit.Parameters()
params.add('strainmax', value=strainmax, min=0.004, max=0.02)
params.add('straindepth', value=straindepth)
params.add('scale', value=scale)
params.add('thetaoffset', value=thetaoffset, max=0.0001, min=-0.1)
params.add('m', value=m)
# 
# residual(params, Energy, data, Q)
# 
# out = lmfit.minimize(residual, params, args=(Energy, data, Q))
print params