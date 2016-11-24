"""Auxiliary functions for dynXRD"""
import pyasf
import math
import numpy as np
import sympy as sp


def get_crystal(element):
    return pyasf.unit_cell(element)

class crystal(object):
    """ Class containing symbols for lattice parameters and rotation matrix """
    def __init__(self, cryst):
        self.structure=cryst
        self.a=cryst.a
        self.b=cryst.b
        self.c=cryst.c
        self.alpha=cryst.alpha
        self.beta=cryst.beta
        self.gamma=cryst.gamma

        #self.lattice_par_sym=lattice_par_sym # strings -> symbols (angles in radians)
        self.lattice_par_val=cryst.subs # symbols -> values (angles in radians)
        self.M_matrix=cryst.M
        self.volume=cryst.V

    def FH_structure_factor(self, Energy, Miller):
        #debye_temp_factor=1.0
        #rel_angle=1.0
        FH=self.structure.DAFS(Energy, Miller)
        if not isinstance(FH, sp.numbers.Zero):
            FH = FH[0]
        return FH

def makefunc(expr, mathmodule = "numpy"):
    return pyasf.makefunc(expr, mathmodule=mathmodule)

dictcall = lambda self, d: self.__call__(*[d.get(k, d.get(k.name, k)) for k in self.kw])
