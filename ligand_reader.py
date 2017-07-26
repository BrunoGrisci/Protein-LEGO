#Bruno Iochins Grisci
#JULY/2017
#This code gets the hetatoms positions from a .pdb file

import copy
import numpy as np
import math

from pdb_reader import PDB_reader

class LIGAND_reader(PDB_reader):
    ATOM_TAG = "HETATM"
    END_TAG = "END"
    
    def write_pdbqt(self, file_name):
        j = 0
        f = open(file_name, "r+")        
        lines = f.readlines()
        f.seek(0)
        f.truncate()
        for line in lines:
            if self.ATOM_TAG in line:
                i = line[12:16].replace(" ", "")
                i = int(i[1:]) - 1
                at = line[77:80].replace(" ", "").replace('\n', '')
                line = "{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}".format(self.ATOM_TAG, j+1, self.atoms_full[i], " ", self.amino_acids[i], " ", self.amino_acids_number[i], " ", self.atoms_pos[i][0], self.atoms_pos[i][1], self.atoms_pos[i][2], self.more_stuff[i][0], self.more_stuff[i][1], " ", at)
                line = line + "\n"
                j += 1
            f.write(line)                 
        f.close()   
    
    def rotate_to(self, angles): 
        rotable_bounds = [( 7,  1,  8,  9), 
                          ( 2,  3, 16, 17),
                          ( 3,  4, 23, 24),
                          ( 6,  7, 33, 34),
                          ( 1,  8,  9, 10),
                          ( 3, 16, 17, 18),
                          ( 4, 23, 24, 25),
                          (23, 24, 25, 26),
                          ( 7, 33, 34, 35),
                          (33, 34, 35, 36)]
        rotable_atoms = [[ 9, 14],
                         [17, 22],
                         [24, 30],
                         [34, 40],
                         [10, 14],
                         [18, 22],
                         [25, 30],
                         [26, 30],
                         [35, 40],
                         [36, 40]]
    
        for i in xrange(len(rotable_bounds)):
            current_angle = self.calc_angles(self.atoms_pos[rotable_bounds[i][0]-1], self.atoms_pos[rotable_bounds[i][1]-1], self.atoms_pos[rotable_bounds[i][2]-1], self.atoms_pos[rotable_bounds[i][3]-1])   
            dtheta = self.angle_diff(current_angle, angles[i])
            for a in xrange(rotable_atoms[i][0]-1, rotable_atoms[i][1]):
                self.atoms_pos[a] = self.rotate_atom_around_bond(dtheta, self.atoms_pos[a], self.atoms_pos[rotable_bounds[i][1]-1], self.atoms_pos[rotable_bounds[i][2]-1])     
