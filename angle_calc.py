#Bruno Iochins Grisci
#JUNE/2017

import sys
import numpy as np
import math
from pdb_reader import PDB_reader

def calc_angles(pos1, pos2, pos3):
    pos1 = np.array(pos1)
    pos2 = np.array(pos2)
    pos3 = np.array(pos3)
    
    bond12 = pos2 - pos1
    bond23 = pos3 - pos2
    
    cos_phi = np.inner(pos1, pos2) / (math.sqrt(np.inner(pos1, pos1)) * math.sqrt(np.inner(pos2, pos2)))
    cos_psi = np.inner(pos2, pos3) / (math.sqrt(np.inner(pos2, pos2)) * math.sqrt(np.inner(pos3, pos3)))
    
    #print(cos_phi, cos_psi)
    
    phi = np.degrees(np.arccos(cos_phi))
    psi = np.degrees(np.arccos(cos_psi))
    
    return phi, psi

pdb_file = sys.argv[1]
pdb = PDB_reader(pdb_file)

ca = pdb.get_ca_info()
n = pdb.get_N_info()
c = pdb.get_C_info()

print(len(ca), len(n), len(c))

print("A_A   PHI    PSI")
print("----------------")

'''for aa in xrange(len(ca)):
    name   = ca[aa][1]
    ca_pos = ca[aa][2]
    n_pos  =  n[aa][2]
    c_pos  =  c[aa][2]'''

for amino_acid in zip(ca, n, c):
    name = amino_acid[0][1]
    ca_pos = amino_acid[0][2]
    n_pos  = amino_acid[1][2]
    c_pos  = amino_acid[2][2]
    
    phi, psi = calc_angles(n_pos, ca_pos, c_pos)
    
    print(name + "    " + str(phi) + "    " + str(psi))
    
     
     
