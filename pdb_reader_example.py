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
    print(math.degrees(angles[a]), math.degrees(angles[a+1]))
    
pis = [math.radians(45.0)]*len(angles)    
pdb.rotate_to(pis)

print("########")
angles = pdb.get_angles()
for a in xrange(0, naa*2, 2):
    print(math.degrees(angles[a]), math.degrees(angles[a+1]))
    
pdb.write_pdb("test.pdb")
