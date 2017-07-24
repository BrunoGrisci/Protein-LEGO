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
