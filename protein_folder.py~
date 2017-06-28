#Bruno Iochins Grisci
#JUNE/2017
#Usage: python protein_folder.py reference_file mobile_file

import math
import numpy as np
import copy
import sys
import time
import matplotlib.pyplot as plt
import rmsd

from PSO import PSO
from pdb_reader import PDB_reader

def calc_exact_rmsd(ref_atoms, mob_atoms):
    ra = np.array(copy.deepcopy(ref_atoms))
    ma = np.array(copy.deepcopy(mob_atoms))
    ra -= rmsd.centroid(ra)
    ma -= rmsd.centroid(ma)
    return rmsd.kabsch_rmsd(ra, ma)

def evaluator(solutions, pdb_ref, pdb_mob):
    scores = []
    for solution in solutions:
        pdb_mob.rotate_to([2.0*math.pi] + solution + [2.0*math.pi])
        score = calc_exact_rmsd(pdb_ref.get_all_pos(), pdb_mob.get_all_pos())
        scores.append(score)
    return scores

reference_file = sys.argv[1]
mobile_file = sys.argv[2]

pdb_ref = PDB_reader(reference_file)
pdb_mob = PDB_reader(mobile_file)

pdb_mob.match_atoms(pdb_ref.get_atoms(), pdb_ref.get_amino_acids_index())
pdb_mob.remove_nones()    
pdb_ref.match_atoms(pdb_mob.get_atoms(), pdb_mob.get_amino_acids_index())
pdb_ref.remove_nones()
   
pdb_mob.rotate_omegas()
pdb_mob.set_peptide_bond_angles()
    
pop_size = 200    
dim = 2 * pdb_ref.get_number_amino_acids() - 2
min_iterations = 300

start_time = time.time()

pso = PSO(swarm_size=pop_size, dimensions=dim, lower_bounds=-math.pi, upper_bounds=math.pi, minimization=True)
#Start the regular loop 
scores_over_time = []
for i in xrange(min_iterations):
    locations = pso.get_locations()
    scores = evaluator(locations, pdb_ref, pdb_mob)
    pso.run_step(scores)
    print(pso.get_best_score())
    scores_over_time.append(pso.get_best_score())
elapsed_time = time.time() - start_time
print("Elapsed time: " + str(elapsed_time))
print([2.0*math.pi] + pso.get_best_location() + [2.0*math.pi])
location_in_degrees = []
for loc in [2.0*math.pi] + pso.get_best_location() + [2.0*math.pi]:
    location_in_degrees.append(math.degrees(loc))
print(location_in_degrees)
print(pso.get_best_score())
print("Finished run")
   
pdb_mob.rotate_to([2.0*math.pi] + pso.get_best_location() + [2.0*math.pi])
pdb_mob.write_pdb(reference_file.replace(".pdb", "-F.pdb"))

print("###")
angles = pdb_mob.get_angles()
for a in xrange(0, pdb_mob.get_number_amino_acids()*2, 2):
    print(round(math.degrees(angles[a]),2), round(math.degrees(angles[a+1]),2))

plt.plot(scores_over_time)
plt.title('RMSD over iterations: ' + reference_file.replace(".pdb", "-F"))
plt.xlabel('Iterations')
plt.ylabel('RMSD')
fig = plt.gcf()
#fig.set_size_inches(10, 10)
fig.savefig(reference_file.replace(".pdb", "-F_rmsd.png"), dpi=100, bbox_inches='tight') 
