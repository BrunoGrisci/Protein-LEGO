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

pdb_ref = PDB_reader(reference_file)
pdb_mob = PDB_reader(mobile_file)

oa = pdb_mob.get_atoms()

print(len(pdb_ref.get_atoms()), len(pdb_mob.get_atoms()))

pdb_mob.match_atoms(pdb_ref.get_atoms(), pdb_ref.get_amino_acids_index())
pdb_mob.remove_nones()
pdb_ref.match_atoms(pdb_mob.get_atoms(), pdb_mob.get_amino_acids_index())
pdb_ref.remove_nones()

print(len(pdb_ref.get_atoms()), len(pdb_mob.get_atoms()))

for atom in zip(pdb_ref.get_atoms(), pdb_mob.get_atoms()):
    print(atom)

rmsd = calc_exact_rmsd(pdb_ref.get_all_pos(), pdb_mob.get_all_pos())    
print(rmsd)
