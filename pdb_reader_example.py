#Bruno Iochins Grisci
#JUNE/2017

import sys
import math
from pdb_reader import PDB_reader

pdb_file = sys.argv[1]
pdb = PDB_reader(pdb_file)

naa = pdb.get_number_amino_acids()
print(naa)

ca = pdb.get_ca_info()
print(ca)

n = pdb.get_N_info()
print(n)

c = pdb.get_C_info()
print(c)

print(len(ca), len(n), len(c))

angles = pdb.get_angles()
for a in xrange(0, naa*2, 2):
    print(round(math.degrees(angles[a]),2), round(math.degrees(angles[a+1]),2))
    
#pis = [math.radians(90.0)]*len(angles)    
pis = [math.radians(360.0),math.radians(176.63),
       math.radians(148.48),math.radians(-21.96),
       math.radians(114.02),math.radians(29.89),
       math.radians(-88.0),math.radians(-38.16),
       math.radians(-74.24),math.radians(360.0)]

print("example start")
print(len(pdb.get_all_pos()))
for p in pdb.get_all_pos():
    print(p)

pdb.rotate_to(pis)

print("########")
angles = pdb.get_angles()
for a in xrange(0, naa*2, 2):
    print(round(math.degrees(angles[a]),2), round(math.degrees(angles[a+1]),2))

print("example end")
print(len(pdb.get_all_pos()))
for p in pdb.get_all_pos():
    print(p)


    
pdb.write_pdb("test.pdb")
