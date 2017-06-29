#Bruno Iochins Grisci
#JUNE/2017

import sys
import numpy as np
import math

from aa_reader import AA_reader
from pdb_reader import PDB_reader

def write_pdb(file_name, amino_acids):
    pdb = open(file_name, "w")
    for aa in amino_acids:
        for atom in aa.get_atoms():
            line = "{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}".format(atom[0], atom[1], atom[10], " ", atom[3], " ", atom[4], " ", atom[5], atom[6], atom[7], atom[8], atom[9], " ", " ")
            line = line + "\n"
            pdb.write(line)
    pdb.write("TER")
    pdb.close()

AMINO_ACIDS = {'A':'amino_acids/alanine.pdb', 
               'R':'amino_acids/arginine.pdb', 
               'N':'amino_acids/asparagine.pdb', 
               'D':'amino_acids/aspartic_acid.pdb', 
               'ASX':'B', 
               'C':'amino_acids/cysteine.pdb', 
               'E':'amino_acids/glutamic_acid.pdb', 
               'Q':'amino_acids/glutamine.pdb', 
               'GLX':'Z', 
               'G':'amino_acids/glycine.pdb', 
               'H':'amino_acids/histidine.pdb', 
               'I':'amino_acids/isoleucine.pdb', 
               'L':'amino_acids/leucine.pdb', 
               'K':'amino_acids/lysine.pdb', 
               'M':'amino_acids/methionine.pdb', 
               'F':'amino_acids/phenalalanine.pdb', 
               'P':'amino_acids/proline.pdb', 
               'S':'amino_acids/serine.pdb', 
               'T':'amino_acids/threonine.pdb', 
               'W':'amino_acids/tryptophan.pdb', 
               'Y':'amino_acids/tyrosine.pdb', 
               'V':'amino_acids/valine.pdb'}

PPL = 1.32 #peptide bound length = 1.32A
               
aa_sequence = sys.argv[1]

amino_acids = []

for aa in aa_sequence: 
    amino_acids.append(AA_reader(AMINO_ACIDS[aa]))

for aa in amino_acids:
    aa.send_origin()   

print("###################################################")

aai = 1
if len(amino_acids) > 1:
    OC_pos = amino_acids[0].get_pos_OC()
    C_pos = amino_acids[0].get_pos_C()
    amino_acids[0].remove_OH()
    amino_acids[0].set_amino_acid_index(aai)

if len(amino_acids) > 2:
    for aa in amino_acids[1:-1]:
        aa.remove_H()
        distance_vector = np.array(OC_pos) - np.array(C_pos)
        distance = np.linalg.norm(distance_vector)
        r = PPL/distance   #peptide bound length = 1.32A    
        aa.relocate_N(list(np.array(C_pos) + r * distance_vector))
        C_pos = aa.get_pos_C()
        OC_pos = aa.get_pos_OC()
        aa.remove_OH()
        aai += 1
        aa.set_amino_acid_index(aai)

if len(amino_acids) > 1:
    amino_acids[-1].remove_H()
    amino_acids[-1].relocate_N(list(np.array(C_pos) + r * distance_vector))
    amino_acids[-1].set_amino_acid_index(aai+1)

index = 1
for aa in amino_acids:
    index = aa.set_atoms_index(index)

for aa in amino_acids:    
    for atom in aa.get_atoms():
        print(atom)

write_pdb("data/" + aa_sequence + ".pdb", amino_acids)

pdb = PDB_reader("data/" + aa_sequence + ".pdb")
pdb.set_peptide_bond_angles()
pdb.set_NH_angles()
pdb.rotate_omegas()
pis = [math.radians(180.0)]*(len(aa_sequence)*2)
pdb.rotate_to(pis)
pdb.write_pdb("data/" + aa_sequence + ".pdb")
