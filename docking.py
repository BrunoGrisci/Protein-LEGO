#Bruno Iochins Grisci
#JULY/2017

import sys
import math

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
        data = p[0] + [p[1]]
        line = str(pid)
        for d in data:
            line = line + ',' + str(d)
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

def evaluator(solutions, pdb_pro, pdb_lig):
    scores = []
    for solution in solutions:
        pdb_lig.translate(solution[0], solution[1], solution[2])
        pdb_lig.rotate(solution[3], solution[4], solution[5])
        score = 0
        scores.append(score)
    return scores

protein_file = sys.argv[1]
ligand_file = sys.argv[2]

pdb_pro = PDB_reader(protein_file)
pdb_lig = LIGAND_reader(ligand_file)

   
