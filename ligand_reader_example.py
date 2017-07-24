#Bruno Iochins Grisci
#JULY/2017

import sys
import math
from ligand_reader import LIGAND_reader

pdb_file = sys.argv[1]
pdb = LIGAND_reader(pdb_file)

print(pdb.get_sequence())
print(pdb.get_atoms())
print(pdb.get_amino_acids())

naa = pdb.get_number_amino_acids()
print(naa)

ca = pdb.get_ca_info()
print(ca)

n = pdb.get_N_info()
print(n)

c = pdb.get_C_info()
print(c)

print(len(ca), len(n), len(c))

#angles = pdb.get_angles()
#for a in xrange(0, naa*2, 2):
#    print(round(math.degrees(angles[a]),2), round(math.degrees(angles[a+1]),2))
