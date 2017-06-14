#Bruno Iochins Grisci
#JUNE/2017
#Usage: python protein_folder.py reference_file mobile_file

import math
import numpy as np
import copy
import sys
import rmsd

from PSO import PSO
from pdb_reader import PDB_reader

def calc_exact_rmsd(ref_atoms, mob_atoms):
    ra = np.array(copy.deepcopy(ref_atoms))
    ma = np.array(copy.deepcopy(mob_atoms))
    ra -= rmsd.centroid(ra)
    ma -= rmsd.centroid(ma)
    return rmsd.kabsch_rmsd(ra, ma)

def evaluator(solutions, pdb_ref, pdbs_mob):
    scores = []
    for solution in zip(solutions, pdbs_mob):
        solution[1].rotate_to(solution[0]):
        score = calc_exact_rmsd(pdb_ref.get_all_pos(), solution[1].get_all_pos())
        scores.append(score)
    return scores

reference_file = sys.argv[1]
mobile_file = sys.argv[2]

pdb_ref = PDB_reader(reference_file)
    
pop_size = 200    
dim = 2 * pdb_ref.get_number_amino_acids()
min_iterations = 1000    

mobs = []
for i in xrange(pop_size):
    mobs.append(PDB_reader(mobile_file))
for mob in mobs:
    mob.match_atoms(pdb_ref.get_atoms(), pdb_ref.get_amino_acids())
    mob.remove(nones)    
pdb_ref.match_atoms(mobs[0].get_atoms(), mobs[0].get_amino_acids())
pdb_ref.remove_nones()

pso = PSO(swarm_size=pop_size, dimensions=dim, lower_bounds=0.0, upper_bounds=2.0*math.pi, minimization=True)

#Start the regular loop 
for i in xrange(min_iterations):
    locations = pso.get_locations()
    scores = evaluator(locations, pdb_ref, mobs)
    pso.run_step(scores)
print(pso.get_best_location())
print(pso.get_best_score())
print("Finished run")
    
print(pso.get_best_location())
print(pso.get_best_score())
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
