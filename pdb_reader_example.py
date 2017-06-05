#Bruno Iochins Grisci
#MAY/2017
#This code is a usage example of the PDB_reader class

from pdb_reader import PDB_reader

alanine = PDB_reader("amino_acids/alanine.pdb")
print(len(alanine.get_atoms()))
for a in alanine.get_atoms():
    print(a)
   
alanine.send_origin()
print(len(alanine.get_atoms())) 
for a in alanine.get_atoms():
    print(a)

ala_pos = alanine.get_pos_OC()
print(ala_pos)

alanine.remove_H()
print(len(alanine.get_atoms())) 
for a in alanine.get_atoms():
    print(a)
    
alanine.remove_OH()
print(len(alanine.get_atoms())) 
for a in alanine.get_atoms():
    print(a)
    
alanine.relocate_N(ala_pos)
print(len(alanine.get_atoms())) 
for a in alanine.get_atoms():
    print(a)
    
alanine.set_amino_acid_index(5)
print(len(alanine.get_atoms())) 
for a in alanine.get_atoms():
    print(a)
