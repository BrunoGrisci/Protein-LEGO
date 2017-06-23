#Bruno Iochins Grisci
#JUNE/2017

import sys
import numpy as np
import math
import matplotlib.pyplot as plt

from pdb_reader import PDB_reader

pdb_file = sys.argv[1]
pdb = PDB_reader(pdb_file)

print("A_A   PHI     PSI        OMEGA")
print("------------------------------")

names = pdb.get_sequence()
phis = pdb.get_angles()[::2]
psis = pdb.get_angles()[1::2]
omegas = pdb.get_omegas()

angle_file = open(pdb_file.replace('.pdb', '_angles.txt'), 'w')

for a in xrange(len(names)):
    line = "{:3s}  {:7.2f}  {:7.2f}  {:7.2f}".format(names[a], math.degrees(phis[a]), math.degrees(psis[a]), math.degrees(omegas[a]))
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
    
     
     
