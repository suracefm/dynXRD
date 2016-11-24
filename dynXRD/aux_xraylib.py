"""Auxiliary functions for dynXRD"""
import xraylib
import math
import numpy as np
import sympy as sp
from sympy.utilities import lambdify
import types
from operator import itemgetter


def get_crystal(element):
    return xraylib.Crystal_GetCrystal(element)

class crystal(object):
    """ Class containing symbols for lattice parameters and rotation matrix """
    def __init__(self, cryst):
        self.structure=cryst
        lattice_par_sym={}

        a=sp.symbols('a')
        lattice_par_sym['a']=a
        if cryst['b']==cryst['a']:
            lattice_par_sym['b']=a
        else:
            b=sp.symbols('b')
            lattice_par_sym['b']=b
        if cryst['c']==cryst['a']:
            lattice_par_sym['c']=a
        elif cryst['c']==cryst['b']:
            lattice_par_sym['c']=b
        else:
            c=sp.symbols('c')
            lattice_par_sym['c']=c

        alpha=sp.symbols('alpha')
        lattice_par_sym['alpha']=alpha
        if cryst['beta']==cryst['alpha']:
            lattice_par_sym['beta']=alpha
        else:
            beta=sp.symbols('beta')
            lattice_par_sym['alpha']=beta
        if cryst['gamma']==cryst['alpha']:
            lattice_par_sym['gamma']=alpha
        elif cryst['gamma']==cryst['beta']:
            lattice_par_sym['gamma']=beta
        else:
            gamma=sp.symbols('gamma')
            lattice_par_sym['gamma']=gamma

        lattice_par_val={key:cryst[str(key)] for key in lattice_par_sym.values()}

        angles={'alpha', 'beta', 'gamma'}
        angles_sym=list(set([lattice_par_sym[x] for x in angles]))
        special_angles={90.0:sp.pi/2, 120.0:2*sp.pi/3, 60.0:sp.pi/3, 45.0:sp.pi/4}

        for angle in angles_sym:
            if cryst[str(angle)] in special_angles.keys():
                for x in lattice_par_sym.keys():
                    lattice_par_sym[x]=lattice_par_sym[x].subs(angle, special_angles[cryst[str(angle)]])

        for angle in angles_sym:
            lattice_par_val[angle]=lattice_par_val[angle]*np.pi/180

        #symbols
        self.a=lattice_par_sym['a']
        self.b=lattice_par_sym['b']
        self.c=lattice_par_sym['c']
        self.alpha=lattice_par_sym['alpha']
        self.beta=lattice_par_sym['beta']
        self.gamma=lattice_par_sym['gamma']

        self.lattice_par_sym=lattice_par_sym # strings -> symbols (angles in radians)
        self.lattice_par_val=lattice_par_val # symbols -> values (angles in radians)
        self.get_M_matrix()
        self.volume=sp.sqrt(-self.a**2*self.b**2*self.c**2*sp.cos(self.alpha)**2 + 2*self.a**2*self.b**2*self.c**2*\
        sp.cos(self.alpha)*sp.cos(self.beta)*sp.cos(self.gamma) - self.a**2*self.b**2*self.c**2*sp.cos(self.beta)**2 - \
        self.a**2*self.b**2*self.c**2*sp.cos(self.gamma)**2 + self.a**2*self.b**2*self.c**2)

    def get_M_matrix(self):
        """ Calculate M matrix, i.e. matrix changing bases from primitive lattice
        vectors to x, y, z (with x pointing along a)"""
        a, b, c, alpha, beta, gamma = itemgetter('a', 'b', 'c', 'alpha', 'beta', 'gamma')(self.lattice_par_sym)
        M=sp.Matrix([[a, b*sp.cos(gamma), c*sp.cos(beta)],\
         [0, b*sp.sin(gamma), c*(sp.cos(alpha) - sp.cos(beta)*sp.cos(gamma))/sp.sin(gamma)],\
         [0, 0, c*sp.sqrt(-(sp.cos(alpha)*sp.cos(beta) - sp.cos(gamma))**2\
         /(sp.sin(alpha)**2*sp.sin(beta)**2) + 1)*sp.sin(alpha)*sp.sin(beta)/sp.sin(gamma)]])
        self.M_matrix=M
        return M

    def FH_structure_factor(self, Energy, Miller):
        #debye_temp_factor=1.0
        #rel_angle=1.0
        return xraylib.Crystal_F_H_StructureFactor(self.structure, Energy, Miller[0], Miller[1], Miller[2], 1.0, 1.0)

def makefunc(expr, mathmodule = "numpy", dummify=False):
    symbols = list(expr.atoms(sp.Symbol))
    symbols.sort(key=str)
    func = lambdify(symbols, expr, mathmodule, dummify=dummify)
    func.kw = symbols
    func.expr = expr
    func.kwstr = map(lambda x: x.name, symbols)
    func.dictcall = types.MethodType(dictcall, func)
    func.__doc__ = str(expr)
    return func

dictcall = lambda self, d: self.__call__(*[d.get(k, d.get(k.name, k)) for k in self.kw])
