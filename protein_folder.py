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

def calc_exact_rmsd(solutions, ref_atoms, mob_atoms):
    ra = np.array(copy.deepcopy(ref_atoms))
    ma = np.array(copy.deepcopy(mob_atoms))
    ra -= rmsd.centroid(ra)
    ma -= rmsd.centroid(ma)
    return rmsd.kabsch_rmsd(ra, ma)

def evaluator(solutions, ref_atoms, mob_atoms):
    scores = []
    for solution in solutions:
        score = rmsd(solution, ref_atoms, mob_atoms)
        scores.append(score)
    return scores

reference_file = sys.argv[1]
mobile_file = sys.argv[2]

pop_size = 200
dim = 6
min_iterations = 1000
initial_lr = 1.0
lr = initial_lr
    
pdb_ref = PDB_reader(reference_file)
pdb_mob = PDB_reader(mobile_file)
pdb_mob.match_atoms(pdb_ref.get_atoms(), pdb_ref.get_amino_acids())

latest_best_scores = []
    
pso = PSO(swarm_size=pop_size, dimensions=dim, lower_bounds=0.0, upper_bounds=2.0*math.pi, movement_step=initial_lr, minimization=True)

#Start with regular loop 
for i in xrange(min_iterations):
    locations = pso.get_locations()
    scores = evaluator(locations, pdb_ref.get_all_pos(), pdb_mob.get_all_pos())
    pso.run_step(scores, initial_lr)
print(pso.get_best_location())
print(pso.get_best_score())
print("Finished run")
    
print(pso.get_best_location())
print(pso.get_best_score())
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
