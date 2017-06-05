#Bruno Iochins Grisci
#JUNE/2017

import sys
from pdb_reader import PDB_reader

def write_pdb(file_name, amino_acids):
    pdb = open(file_name, "w")
    for aa in amino_acids:
        for atom in aa.get_atoms():
            line = ""
            for info in atom:
                line = line + str(info) + "    "
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
               
aa_sequence = sys.argv[1]

amino_acids = []

for aa in aa_sequence: 
    amino_acids.append(PDB_reader(AMINO_ACIDS[aa]))

for aa in amino_acids:
    aa.send_origin()
for aa in amino_acids:    
    for atom in aa.get_atoms():
        print(atom)    

print("###################################################")

aai = 1
OC_pos = amino_acids[0].get_pos_OC()
amino_acids[0].remove_OH()
amino_acids[0].set_amino_acid_index(aai)

for aa in amino_acids[1:-1]:
    aa.remove_H()
    aa.relocate_N(OC_pos)
    OC_pos = aa.get_pos_OC()
    aa.remove_OH()
    aai += 1
    aa.set_amino_acid_index(aai)

amino_acids[-1].remove_H()
amino_acids[-1].relocate_N(OC_pos)
amino_acids[-1].set_amino_acid_index(aai+1)

for aa in amino_acids:    
    for atom in aa.get_atoms():
        print(atom)

write_pdb(aa_sequence + ".pdb", amino_acids)
