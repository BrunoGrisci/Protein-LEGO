#Bruno Iochins Grisci
#JULY/2017

import sys
import math
import random
import numpy as np
import copy
import time
import matplotlib.pyplot as plt

from pdb_reader import PDB_reader
from ligand_reader import LIGAND_reader
from PSO import PSO

def save_log(i, dimensions, solutions, scrs, bsol, bscr):
    header = 'id'
    for d in xrange(dimensions):
        header = header + ',d' + str(d)
    header = header + ',score\n'
    log = open('log/' + str(i) + '.csv', 'w')
    log.write(header)
    pid = 1
    for p in zip(solutions, scrs):
        line = str(pid)
        for d in p[0][0:3]:
            line = line + ',' + str(d)
        for d in p[0][3:]:
            line = line + ',' + str(math.degrees(d))  
        line = line + ',' + str(p[1])          
        line = line + '\n'
        log.write(line)
        pid += 1
    line = 'best'
    for d in bsol[0:3]:
        line = line + ',' + str(d)
    for d in bsol[3:]:
        line = line + ',' + str(math.degrees(d))
    line = line + ',' + str(bscr)    
    log.write(line)            
    log.close() 

def calc_rmsd(ref_atoms, mob_atoms):
    #not aligned
    distance_sum = 0.0
    for coord in zip(ref_atoms, mob_atoms):
        distance_sum += (coord[0][0] - coord[1][0])**2
        distance_sum += (coord[0][1] - coord[1][1])**2 
        distance_sum += (coord[0][2] - coord[1][2])**2  
    score = math.sqrt(distance_sum / float(len(ref_atoms)))
    return score

def evaluator(solutions, pdb_pro, pdb_lig, pdb_ref):
    scores = []
    for solution in solutions:
        pdb_lig1 = copy.deepcopy(pdb_lig)
        pdb_lig1.translate(solution[0:3])
        pdb_lig1.rotate(solution[3:6])
        pdb_lig1.rotate_to(solution[6:])
        score = calc_rmsd(pdb_ref.get_all_pos(), pdb_lig1.get_all_pos())
        scores.append(score)
    return scores

protein_file = sys.argv[1]
ligand_file = sys.argv[2]

pdb_pro = PDB_reader(protein_file)
pdb_lig = LIGAND_reader(ligand_file)
pdb_ref = LIGAND_reader(ligand_file)

pdb_lig.translate([ random.uniform(-50, 50) for i in range(3) ])
pdb_lig.rotate([ random.uniform(-math.pi, math.pi) for i in range(3) ])   
pdb_lig.rotate_to([ random.uniform(-math.pi, math.pi) for i in range(10) ])

e_ref = 0.0
print('Reference energy: ' + str(e_ref))
    
pop_size = 60
dim = 3 + 3 + 10
min_iterations = 1000   
lb = [-50.0]*3 + [-math.pi]*3 + [-math.pi]*10
ub = [50.0]*3 + [math.pi]*3 + [math.pi]*10

start_time = time.time()

pso = PSO(swarm_size=pop_size, dimensions=dim, lower_bounds=lb, upper_bounds=ub, minimization=True)
#Start the regular loop 
scores_over_time = []
rmsds_over_time = []
for i in xrange(min_iterations):
    locations = pso.get_locations()
    scores = evaluator(locations, pdb_pro, pdb_lig, pdb_ref)
    
    #if i > 0 and i%1 == 0:
    #    save_log(i, dim, locations, scores, pso.get_best_location(), pso.get_best_score())
    
    pso.run_step(scores) 
    scores_over_time.append(pso.get_best_score())
    pdb_lig1 = copy.deepcopy(pdb_lig)
    pdb_lig1.translate(pso.get_best_location()[0:3])
    pdb_lig1.rotate(pso.get_best_location()[3:6])
    pdb_lig1.rotate_to(pso.get_best_location()[6:])
    best_rmsd = calc_rmsd(pdb_ref.get_all_pos(), pdb_lig1.get_all_pos())
    rmsds_over_time.append(best_rmsd)
    print(i+1, pso.get_best_score(), best_rmsd, pso.get_best_score() - e_ref)
    
    '''if (i+1)%20 == 0:
        print("###############" + str(i+1) + "##############")   
        elapsed_time = time.time() - start_time
        print("Elapsed time: " + str(elapsed_time))
        print(pso.get_best_location())
        location_in_degrees = []
        for loc in pso.get_best_location()[3:-1]:
            location_in_degrees.append(math.degrees(loc))
        print(pso.get_best_location()[0:3] + location_in_degrees)
        print(pso.get_best_score())
        
        pdb_lig1 = copy.deepcopy(pdb_lig)
        pdb_lig1.translate(pso.get_best_location()[0:3])
        pdb_lig1.rotate(pso.get_best_location()[3:6])
        pdb_lig1.rotate_to(pso.get_best_location()[6:])
        pdb_lig1.write_pdb(ligand_file.replace(".pdb", "-F" + str(i+1) + ".pdb"))

        fig, ax1 = plt.subplots()
        ax1.plot(scores_over_time, 'b-')
        ax1.set_title('Energy over iterations: ' + ligand_file.replace(".pdb", "-F"))
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
        fig.savefig(ligand_file.replace(".pdb", "-F" + str(i+1) + "_energy.png"), dpi=100, bbox_inches='tight')'''

print("Finished run")
