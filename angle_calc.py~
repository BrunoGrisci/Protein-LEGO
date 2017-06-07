#Bruno Iochins Grisci
#JUNE/2017

import sys
import numpy as np
import math
from pdb_reader import PDB_reader

def calc_angles(pos1, pos2, pos3, pos4):
    pos1 = np.array(pos1)
    pos2 = np.array(pos2)
    pos3 = np.array(pos3)
    pos4 = np.array(pos4)
    
    bond12 = pos2 - pos1
    bond23 = pos3 - pos2
    bond34 = pos4 - pos3
    
    
    
    cos_angle = np.inner(bond12, bond34) / (math.sqrt(np.inner(bond12, bond12)) * math.sqrt(np.inner(bond34, bond34)))
    angle = np.degrees(np.arccos(cos_angle))
    return round(angle, 2)

pdb_file = sys.argv[1]
pdb = PDB_reader(pdb_file)

ca = pdb.get_ca_info()
n = pdb.get_N_info()
c = pdb.get_C_info()

print(len(ca), len(n), len(c))

print("A_A   PHI    PSI")
print("----------------")

for aa in xrange(len(ca)):
    name   = ca[aa][1]
    
    if aa > 0:
        pre_c_pos = c[aa-1][2]
    n_pos  =  n[aa][2]
    ca_pos = ca[aa][2]
    c_pos  =  c[aa][2]
    if aa < len(ca) - 1:
        nex_n_pos = n[aa+1][2]
    
    if aa > 0: 
        phi = calc_angles(pre_c_pos, n_pos, ca_pos, c_pos)
    else:
        phi = 360.0
    
    if aa < len(ca) - 1: 
        psi = calc_angles(n_pos, ca_pos, c_pos, nex_n_pos)
    else:
        psi = 360.0
    print(name + "    " + str(phi) + "    " + str(psi))

'''
for amino_acid in zip(ca, n, c):
    name = amino_acid[0][1]
    ca_pos = amino_acid[0][2]
    n_pos  = amino_acid[1][2]
    c_pos  = amino_acid[2][2]
    
    phi = calc_angles(n_pos, ca_pos, c_pos)
    
    print(name + "    " + str(phi) + "    " + str(psi))'''
    
     
     
