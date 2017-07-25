#Bruno Iochins Grisci
#JULY/2017

import sys
import math
import random

from ligand_reader import LIGAND_reader

pdb_file = sys.argv[1]
pdb = LIGAND_reader(pdb_file)

print(pdb.get_sequence())
print(pdb.get_atoms())
print(pdb.get_amino_acids())
print(pdb.get_all_pos())

naa = pdb.get_number_amino_acids()
print(naa)

'''pdb.translate([4, 4, 4])
pdb.rotate([math.pi/2.0, math.pi/2.0, math.pi/2.0])
pdb.rotate_to([math.pi/2.0]*10)'''

pdb.translate([ random.uniform(-50, 50) for i in range(3) ])
pdb.rotate([ random.uniform(-math.pi, math.pi) for i in range(3) ])   
pdb.rotate_to([ random.uniform(-math.pi, math.pi) for i in range(10) ])

print(pdb.get_all_pos())

pdb.write_pdb(pdb_file.replace(".pdb", "-TEST.pdb"))






