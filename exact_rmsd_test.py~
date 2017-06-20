#Bruno Iochins Grisci
#JUNE/2017
#Usage: python exact_rmsd_test.py reference_file mobile_file atoms_type

import math
import numpy as np
import copy
import sys
import rmsd

from pdb_reader import PDB_reader

def calc_exact_rmsd(ref_atoms, mob_atoms):
    ra = np.array(copy.deepcopy(ref_atoms))
    ma = np.array(copy.deepcopy(mob_atoms))
    ra -= rmsd.centroid(ra)
    ma -= rmsd.centroid(ma)
    return rmsd.kabsch_rmsd(ra, ma)
    
reference_file = sys.argv[1]
mobile_file = sys.argv[2]
atoms_type = sys.argv[3]

pdb_ref = PDB_reader(reference_file)
pdb_mob = PDB_reader(mobile_file)

#pdb_mob.move_to(pdb_ref.get_all_pos()[0])
#pdb_mob.write_pdb("testN.pdb")

print(len(pdb_ref.get_atoms()), len(pdb_mob.get_atoms()))

pdb_mob.match_atoms(pdb_ref.get_atoms(), pdb_ref.get_amino_acids_index())
pdb_mob.remove_nones()
pdb_ref.match_atoms(pdb_mob.get_atoms(), pdb_mob.get_amino_acids_index())
pdb_ref.remove_nones()

print(len(pdb_ref.get_atoms()), len(pdb_mob.get_atoms()))

for atom in zip(pdb_ref.get_atoms(), pdb_mob.get_atoms()):
    print(atom)

RMSD = None
if atoms_type == "all":
    RMSD = calc_exact_rmsd(pdb_ref.get_all_pos(), pdb_mob.get_all_pos())  
elif atoms_type == "backbone":
    RMSD = calc_exact_rmsd(pdb_ref.get_backbone_pos(), pdb_mob.get_backbone_pos())  
elif atoms_type == "ca":
    RMSD = calc_exact_rmsd(pdb_ref.get_ca_pos(), pdb_mob.get_ca_pos())    
print(RMSD)
