import sys
sys.path.append("/home/federica/Documenti/gitrep/dynXRD/")
import pyasf
import reflectivity
import sympy as sp
#import matplotlib


#struct1 = pyasf.unit_cell("cif/2104737.cif") #Si cubic
struct2 = pyasf.unit_cell("1000270") #triclinic
#struct3 = pyasf.unit_cell("1000065") #graphite hexagonal
print(struct2.V)
