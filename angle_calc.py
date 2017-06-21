#Bruno Iochins Grisci
#JUNE/2017

import sys
import numpy as np
import math
import matplotlib.pyplot as plt

from pdb_reader import PDB_reader

def calc_angles(pos1, pos2, pos3, pos4):
    #https://math.stackexchange.com/questions/47059/how-do-i-calculate-a-dihedral-angle-given-cartesian-coordinates

    pos1 = np.array(pos1)
    pos2 = np.array(pos2)
    pos3 = np.array(pos3)
    pos4 = np.array(pos4)
    
    bond12 = pos2 - pos1
    bond23 = pos3 - pos2
    bond34 = pos4 - pos3
    
    normal1 = np.cross(bond12, bond23) / np.linalg.norm(np.cross(bond12, bond23))  
    normal2 = np.cross(bond23, bond34) / np.linalg.norm(np.cross(bond23, bond34))
    m1 = np.cross(normal1, bond23 / np.linalg.norm(bond23))
    
    x = np.inner(normal1, normal2)
    y = np.inner(m1, normal2)
    angle = math.atan2(y, x)
    angle = np.degrees(angle)
    return round(-angle, 2)

pdb_file = sys.argv[1]
pdb = PDB_reader(pdb_file)

ca = pdb.get_ca_info()
n = pdb.get_N_info()
c = pdb.get_C_info()

print(len(ca), len(n), len(c))

print("A_A   PHI     PSI        OMEGA")
print("------------------------------")

phis = []
psis = []
omegas = []

angle_file = open(pdb_file.replace('.pdb', '_angles.txt'), 'w')

for aa in xrange(len(ca)):
    name   = ca[aa][1]
    
    if aa > 0:
        pre_c_pos = c[aa-1][2]
    n_pos  =  n[aa][2]
    ca_pos = ca[aa][2]
    c_pos  =  c[aa][2]
    if aa < len(ca) - 1:
        nex_n_pos  =  n[aa+1][2]
        nex_ca_pos = ca[aa+1][2]
    
    if aa > 0: 
        phi = calc_angles(pre_c_pos, n_pos, ca_pos, c_pos)
    else:
        phi = 360.0
    
    if aa < len(ca) - 1: 
        psi = calc_angles(n_pos, ca_pos, c_pos, nex_n_pos)
        omega = calc_angles(ca_pos, c_pos, nex_n_pos, nex_ca_pos)
    else:
        psi = 360.0
        omega = 360.0
    
    phis.append(phi)
    psis.append(psi)
    omegas.append(omega)
    line = "{:3s}  {:7.2f}  {:7.2f}  {:7.2f}".format(name, phi, psi, omega)
    print(line)
    angle_file.write(line + '\n')
    
angle_file.close()

plt.plot(phis, psis, 'ro')
plt.axis([-180, 180, -180, 180])
plt.grid(True)
plt.title('Ramachandran Map ' + pdb_file)
plt.xlabel('PHI')
plt.ylabel('PSI')
#plt.show()
#plt.savefig(pdb_file.replace('.pdb', '.png'), bbox_inches='tight')
fig = plt.gcf()
fig.set_size_inches(10, 10)
fig.savefig(pdb_file.replace('.pdb', '.png'), dpi=100, bbox_inches='tight')
    
     
     
