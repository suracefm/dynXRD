#
# example of using xraylib to get crystal data
#
# see: 
#       http://dx.doi.org/10.1016/j.sab.2011.09.011   (paper) 
#       https://github.com/tschoonj/xraylib/          (code) 
#       http://lvserver.ugent.be/xraylib-web          (web interface, but crystals not included!)
#


def structure_MgO():
    cryst = {} # cryst.copy()

    cryst['name']   = b"MgO"
    cryst['a']      =  4.217
    cryst['b']      =  4.217
    cryst['c']      =  4.217
    cryst['alpha']  = 90.0
    cryst['beta']   = 90.0
    cryst['gamma']  = 90.0
    cryst['n_atom'] =  8
    cryst['atom']   = [{'x': 0.0,  'y': 0.0,  'z': 0.0,  'Zatom': 12, 'fraction': 1.0, 'label':'Mg'},
                       {'x': 0.0,  'y': 0.5,  'z': 0.5,  'Zatom': 12, 'fraction': 1.0, 'label':'Mg'},
                       {'x': 0.5,  'y': 0.0,  'z': 0.5,  'Zatom': 12, 'fraction': 1.0, 'label':'Mg'},
                       {'x': 0.5,  'y': 0.5,  'z': 0.0,  'Zatom': 12, 'fraction': 1.0, 'label':'Mg'},
                       {'x': 0.5,  'y': 0.5,  'z': 0.5,  'Zatom': 8 , 'fraction': 1.0, 'label':'O'},
                       {'x': 0.5,  'y': 0.0,  'z': 0.0,  'Zatom': 8 , 'fraction': 1.0, 'label':'O'},
                       {'x': 0.0,  'y': 0.5,  'z': 0.0,  'Zatom': 8 , 'fraction': 1.0, 'label':'O'},
                       {'x': 0.0,  'y': 0.0,  'z': 0.5,  'Zatom': 8 , 'fraction': 1.0, 'label':'O'}
                       ]
    cryst['volume'] = 74.99128631299997

    return cryst


# #
# #  import block
# #
# import xraylib
# import numpy as np
# import scipy.constants.codata
#
#
#
#
#
#
# #
# # get crystal data for silicon crystal
# #
# #cryst = xraylib.Crystal_GetCrystal("Si")
# cryst = structure_MgO()
#
# # print some info
# print ("  Unit cell dimensions [A] are %f %f %f" % (cryst['a'],cryst['b'],cryst['c']))
# print ("  Unit cell angles are %f %f %f" % (cryst['alpha'],cryst['beta'],cryst['gamma']))
# print ("  Unit cell volume [A] is %f" % (cryst['volume']) )
#
# #
# # define miller indices and compute dSpacing
# #
# hh = 1
# kk = 1
# ll = 1
# debyeWaller = 1.0
# rel_angle = 1.0  # ratio of (incident angle)/(bragg angle) -> we work at Bragg angle
#
# dspacing = xraylib.Crystal_dSpacing(cryst,hh,kk,ll )
# print("dspacing: %f A \n"%dspacing)
#
# #
# # define energy and get Bragg angle
# #
# ener = 12.398 # keV
# braggAngle = xraylib.Bragg_angle(cryst,ener,hh,kk,ll )
# print("Bragg angle: %f degrees \n"%(braggAngle*180/np.pi))
#
#
# #
# # get the structure factor (at a given energy)
# #
# f0 = xraylib.Crystal_F_H_StructureFactor(cryst, ener, 0, 0, 0, debyeWaller, 1.0)
# fH = xraylib.Crystal_F_H_StructureFactor(cryst, ener, hh, kk, ll, debyeWaller, 1.0)
# print("f0: (%f , %f) \n"%(f0.real,f0.imag))
# print("fH: (%f , %f) \n"%(fH.real,fH.imag))
#
# #
# # convert structure factor in chi (or psi) = - classical_e_radius wavelength^2 fH /(pi volume)
# #
# codata = scipy.constants.codata.physical_constants
# codata_c, tmp1, tmp2 = codata["speed of light in vacuum"]
# codata_h, tmp1, tmp2 = codata["Planck constant"]
# codata_ec, tmp1, tmp2 = codata["elementary charge"]
# codata_r, tmp1, tmp2 = codata["classical electron radius"]
#
# ev2meter = codata_h*codata_c/codata_ec
# wavelength = ev2meter/(ener*1e3)
# print("Photon energy: %f keV \n"%ener)
# print("Photon wavelength: %f A \n"%(1e10*wavelength))
#
# volume = cryst['volume'] *1e-10*1e-10*1e-10 # volume of silicon unit cell in m^3
# cte = - codata_r * wavelength*wavelength/(np.pi * volume)
#
# chi0 = cte*f0
# chiH = cte*fH
#
# print("chi0: (%e , %e) \n"%(chi0.real,chi0.imag))
# print("chiH: (%e , %e) \n"%(chiH.real,chiH.imag))



