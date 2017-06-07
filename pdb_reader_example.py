#Bruno Iochins Grisci
#JUNE/2017

import sys
from pdb_reader import PDB_reader

pdb_file = sys.argv[1]
pdb = PDB_reader(pdb_file)

ca = pdb.get_ca_info()
print(ca)

n = pdb.get_N_info()
print(n)

c = pdb.get_C_info()
print(c)

print(len(ca), len(n), len(c))
