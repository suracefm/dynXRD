import xraylib
import math
import numpy as np
import sympy as sp
import auxfunc


xraylib.XRayInit()
xraylib.SetErrorMessages(0)

Miller = 1,1,1
Energy=10000

structure=auxfunc.get_crystal("Si")
#structure = xraylib.Crystal_GetCrystal("Graphite")
crystal=auxfunc.crystal(structure)
Msym=crystal.get_M_matrix()
Mnum=Msym.subs(crystal.lattice_par_val).evalf()
FH=crystal.FH_structure_factor(Energy, Miller)
print(crystal.lattice_par_sym, crystal.lattice_par_val, Msym, Mnum)
