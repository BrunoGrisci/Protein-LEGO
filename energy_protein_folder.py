#Bruno Iochins Grisci
#JUNE/2017
#Usage: python protein_folder.py reference_file mobile_file bonded_file nonbonded_file aminoacids_file

import math
import numpy as np
import copy
import sys
import time
import matplotlib.pyplot as plt
import rmsd

from PSO import PSO
from pdb_reader import PDB_reader
from force_field import Force_field

def calc_exact_rmsd(ref_atoms, mob_atoms):
    ra = np.array(copy.deepcopy(ref_atoms))
    ma = np.array(copy.deepcopy(mob_atoms))
    ra -= rmsd.centroid(ra)
    ma -= rmsd.centroid(ma)
    return rmsd.kabsch_rmsd(ra, ma)

def evaluator(solutions, force_field, pdb_ref, pdb_mob):
    scores = []
    for solution in solutions:
        pdb_mob.rotate_to([2.0*math.pi] + solution + [2.0*math.pi])
        force_field.insert_PDB(pdb_mob)
        score = force_field.non_bonded()
        scores.append(score)
    return scores

def save_log(i, dimensions, solutions, scrs):
    header = 'id'
    for d in xrange(dimensions):
        header = header + ',d' + str(d)
    header = header + ',score\n'
    log = open('log/' + str(i) + '.csv', 'w')
    log.write(header)
    pid = 1
    for p in zip(solutions, scrs):
        data = p[0] + [p[1]]
        line = str(pid)
        for d in data:
            line = line + ',' + str(d)
        line = line + '\n'
        log.write(line)
        pid += 1
    log.close()    

reference_file = sys.argv[1]
mobile_file = sys.argv[2]
bonded_file = sys.argv[3]
nonbonded_file = sys.argv[4]
aminoacids_file = sys.argv[5]

ff = Force_field(bonded_file, nonbonded_file, aminoacids_file)

pdb_ref = PDB_reader(reference_file)
pdb_mob = PDB_reader(mobile_file)

pdb_mob.match_atoms(pdb_ref.get_atoms(), pdb_ref.get_amino_acids_index())
pdb_mob.remove_nones()    
pdb_ref.match_atoms(pdb_mob.get_atoms(), pdb_mob.get_amino_acids_index())
pdb_ref.remove_nones()
   
pdb_mob.rotate_omegas()
pdb_mob.set_peptide_bond_angles()

ff.insert_PDB(pdb_ref)
e_ref = ff.non_bonded()
print('Reference energy: ' + str(e_ref))
    
pop_size = 60
dim = 2 * pdb_ref.get_number_amino_acids() - 2
min_iterations = 1000   

start_time = time.time()

pso = PSO(swarm_size=pop_size, dimensions=dim, lower_bounds=-math.pi, upper_bounds=math.pi, minimization=True)
#Start the regular loop 
scores_over_time = []
rmsds_over_time = []
for i in xrange(min_iterations):
    locations = pso.get_locations()
    scores = evaluator(locations, ff, pdb_ref, pdb_mob)
    
    if i%1 == 0:
        dlocations = [[math.degrees(x) for x in group] for group in locations]
        save_log(i, dim, dlocations, scores)
    
    pso.run_step(scores) 
    scores_over_time.append(pso.get_best_score())
    pdb_mob.rotate_to([2.0*math.pi] + pso.get_best_location() + [2.0*math.pi])
    best_rmsd = calc_exact_rmsd(pdb_ref.get_all_pos(), pdb_mob.get_all_pos())
    rmsds_over_time.append(best_rmsd)
    print(i+1, pso.get_best_score(), best_rmsd, pso.get_best_score() - e_ref)
    
    if (i+1)%20 == 0:
        print("###############" + str(i+1) + "##############")   
        elapsed_time = time.time() - start_time
        print("Elapsed time: " + str(elapsed_time))
        print([2.0*math.pi] + pso.get_best_location() + [2.0*math.pi])
        location_in_degrees = []
        for loc in [2.0*math.pi] + pso.get_best_location() + [2.0*math.pi]:
            location_in_degrees.append(math.degrees(loc))
        print(location_in_degrees)
        print(pso.get_best_score())
           
        pdb_mob.rotate_to([2.0*math.pi] + pso.get_best_location() + [2.0*math.pi])
        #pdb_mob.write_pdb(reference_file.replace(".pdb", "-F" + str(i+1) + ".pdb"))

        print("###")
        angles = pdb_mob.get_angles()
        for a in xrange(0, pdb_mob.get_number_amino_acids()*2, 2):
            print(round(math.degrees(angles[a]),2), round(math.degrees(angles[a+1]),2))

        '''fig, ax1 = plt.subplots()
        ax1.plot(scores_over_time, 'b-')
        ax1.set_title('Energy over iterations: ' + reference_file.replace(".pdb", "-F"))
        ax1.set_xlabel('Iterations')
        # Make the y-axis label, ticks and tick labels match the line color.
        ax1.set_ylabel('Energy', color='b')
        ax1.tick_params('y', colors='b')
        ax2 = ax1.twinx()
        ax2.plot(rmsds_over_time, 'r-')
        ax2.set_ylabel('RMSD', color='r')
        ax2.tick_params('y', colors='r')
        fig = plt.gcf()
        print("###############################")
        fig.savefig(reference_file.replace(".pdb", "-F" + str(i+1) + "_energy.png"), dpi=100, bbox_inches='tight')'''

print("Finished run")
