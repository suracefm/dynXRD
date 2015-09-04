import os
os.chdir(os.path.expanduser("~/Desktop/reflectivity/dynXRD/"))
import pyasf
import reflectivity
import sympy as sp
import pylab as pl
import numpy as np
import lmfit


strainmax=0.010338974326850121
straindepth=5925.9959587231851
scale=0.0007561892024690888
thetaoffset=0.0072608882909874772

# strainmax=0.010214814368263415 #gaussian
# straindepth=6666.7569700496733
# scale=0.00043979920675453233
# thetaoffset=0.0083924870611093993

# strainmax=0.010277496963932327
# straindepth=6637.1261811733111
# scale=0.00085304866780904886
# thetaoffset=0.0084925890283824579 

data = pl.loadtxt("sto_m28_ez_8_myt1_norm.dat")

R = 0,0,2

def my_strain_profile(strainmax, straindepth, xmax=25, npoints=501):
    x=np.linspace(0, xmax, 501)
    strain=(strainmax*np.exp(-x)*(1+x))[::-1]
    strain=np.delete(strain, 0)
    return x*straindepth, strain

# def my_strain_profile(strainmax, straindepth, xmax=25, npoints=501): #gaussian
#     x=np.linspace(0, xmax, 501)
#     strain=(strainmax*np.exp(-x**2))[::-1]
#     strain=np.delete(strain, 0)
#     return x*straindepth, strain


# def my_strain_profile(strainmax, straindepth, xmax=25, npoints=501): #lorentzian
#     x=np.linspace(0, xmax, 501)
#     strain=(strainmax/(1+x**2))[::-1]
#     strain=np.delete(strain, 0)
#     return x*straindepth, strain


#pl.plot(tvector[1:], strain[::-1])
#pl.show()

Energy=8000


structsub = pyasf.unit_cell("1512124") #BaTiO3
structsub.get_tensor_symmetry()
structsub.build_unit_cell()

structsub.subs[structsub.a] = 3.9064
Sub=reflectivity.Substrate(structsub)

struct1=pyasf.unit_cell("1525437")

struct1.add_atom("Sr", [0,0,0])
struct1.AU_positions.pop("Ba1")
struct1.AU_formfactors.pop("Ba1")
struct1.U.pop("Ba1")
struct1.get_tensor_symmetry()
struct1.build_unit_cell()

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
crystal.calc_g0_gH(Energy)
thBragg= float(Sub.calc_Bragg_angle(Energy).subs(Sub.structure.subs).evalf())
angle=pl.linspace(data[:,0][0]*np.pi/180, data[:,0][-1]*np.pi/180,len(data[:,0]))


XR=crystal.calc_reflectivity(angle, Energy)
crystal.print_values(angle, Energy)



#pl.plot(angle*180/np.pi,abs(layer1.XR)**2)

#pl.plot(angle*180/np.pi,abs(Sub.XR)**2)

pl.plot(data[:,0],abs(XR)**2, color='black', label='calculated')

pl.plot(data[:,0]+thetaoffset, data[:,1]*scale, color='red', label='measured')

pl.yscale('log')
pl.rc('font', size=11)
pl.legend(loc="upper left")
pl.savefig('pics/SrTiO3_2.eps')

pl.show()


def residual(params, angle, data):
    strainmax = params['strainmax'].value
    straindepth = params['straindepth'].value
    scale = params['scale'].value
    thetaoffset = params['thetaoffset'].value
    
    
    
    c = layer1.strain.keys()[0]
    tvector, strain = my_strain_profile(strainmax, straindepth)
    layer1.strain[c] = strain
    layer1.thickness=tvector
    XR=crystal.calc_reflectivity(angle+thetaoffset*np.pi/180, Energy)
    res=scale*data-abs(XR)**2
    print (res**2).sum()
    return res

params = lmfit.Parameters()
params.add('strainmax', value=strainmax, min=0.004, max=0.02)
params.add('straindepth', value=straindepth)
params.add('scale', value=scale, min=1e-4, max=5e-3)
params.add('thetaoffset', value=thetaoffset, max=0.02)

#out = lmfit.minimize(residual, params, args=(angle, data[:,1]))
