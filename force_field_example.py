#Bruno Iochins Grisci
#JUNE/2017

import sys

from force_field import Force_field
from pdb_reader import PDB_reader

bonded_file = sys.argv[1]
nonbonded_file = sys.argv[2]
aminoacids_file = sys.argv[3]
pdb_file = sys.argv[4]

pdb = PDB_reader(pdb_file)
ff = Force_field(bonded_file, nonbonded_file, aminoacids_file)
ff.insert_PDB(pdb)
print(ff.non_bonded())
