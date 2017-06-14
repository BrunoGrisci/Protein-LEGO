#Bruno Iochins Grisci
#MAY/2017
#This code gets the atoms positions from a .pdb file

import copy

class PDB_reader:

    ATOM_TAG = "ATOM"
    END_TAG = "TER"
    ALPHA_CARBON = "CA"
    NITROGEN = "N"
    CARBON = "C"
    BACKBONE_ATOMS = ("N", "CA", "C", "O")
    
    def __init__(self, file_name):
        self.file_name = file_name
        self.atoms = []
        self.atoms_pos = []
        self.amino_acids = []
        self.amino_acids_number = []
        self.read_pdb()
        
    def read_pdb(self):
        stop = False
        pdb = open(self.file_name, "r")
        
        while not stop:
            line = pdb.readline()
            if not line:
                stop = True
            else:
                line = line.split()
                if line[0] == self.END_TAG:
                    stop = True
                elif line[0] == self.ATOM_TAG:
                    atom = line[2]
                    aa = line[3]
                    if self.is_number(line[4]):
                        amino_acid = int(line[4])
                    elif self.is_number(line[5]):
                        amino_acid = int(line[5])
                    pos_init = 0
                    for i in xrange(len(line)):
                        if "." in line[i]:
                            pos_init = i
                            break
                    pos  = map(float, line[pos_init:pos_init+3])
                    self.atoms.append(atom)
                    self.atoms_pos.append(pos)
                    self.amino_acids.append(aa)
                    self.amino_acids_number.append(amino_acid)                  
        pdb.close()
    
    def get_atoms(self):
        return copy.deepcopy(self.atoms)
       
    def get_amino_acids(self):
        return copy.deepcopy(self.amino_acids)       
        
    def get_amino_acids_index(self):
        return copy.deepcopy(self.amino_acids_number)
        
    def get_all_pos(self):
        return copy.deepcopy(self.atoms_pos)
    
    def get_backbone_pos(self):
        backbone_pos = []
        for a in xrange(len(self.atoms)):
            if self.atoms[a] in self.BACKBONE_ATOMS:
                backbone_pos.append(self.atoms_pos[a])
        return copy.deepcopy(backbone_pos) 
        
    def get_ca_pos(self):
        ca_pos = []
        for a in xrange(len(self.atoms)):
            if self.atoms[a] == self.ALPHA_CARBON:
                ca_pos.append(self.atoms_pos[a])
        return copy.deepcopy(ca_pos)
        
    def get_ca_info(self):
        ca_info = []
        for a in zip(self.atoms, self.amino_acids, self.atoms_pos):
            if a[0] == self.ALPHA_CARBON:
                ca_info.append(a)
        return copy.deepcopy(ca_info)
    
    def get_N_info(self):
        n_info = []
        for a in zip(self.atoms, self.amino_acids, self.atoms_pos):
            if a[0] == self.NITROGEN:
                n_info.append(a)
        return copy.deepcopy(n_info)
                
    def get_C_info(self):
        c_info = []
        for a in zip(self.atoms, self.amino_acids, self.atoms_pos):
            if a[0] == self.CARBON:
                c_info.append(a)
        return copy.deepcopy(c_info) 
        
    def match_atoms(self, ref_atoms, ref_amino_acids):
        aux_atoms = []
        aux_pos = []
        aux_aa = []
        aux_aan = []
        for a in xrange(len(ref_atoms)):  
            #print(ref_atoms[a], ref_amino_acids[a])  
            if ref_atoms[a] == None:
                pass   
            elif (ref_atoms[a], ref_amino_acids[a]) in zip(self.atoms, self.amino_acids_number):
                i = zip(self.atoms, self.amino_acids_number).index((ref_atoms[a], ref_amino_acids[a]))
                aux_atoms.append(self.atoms.pop(i))
                aux_pos.append(self.atoms_pos.pop(i))
                aux_aa.append(self.amino_acids_number.pop(i))
                aux_aan.append(self.amino_acids.pop(i))
            elif (ref_atoms[a][1:]+ref_atoms[a][0], ref_amino_acids[a]) in zip(self.atoms, self.amino_acids_number):
                i = zip(self.atoms, self.amino_acids_number).index((ref_atoms[a][1:]+ref_atoms[a][0], ref_amino_acids[a]))
                aux_atoms.append(self.atoms.pop(i))
                aux_pos.append(self.atoms_pos.pop(i))
                aux_aa.append(self.amino_acids_number.pop(i))
                aux_aan.append(self.amino_acids.pop(i))
            elif (ref_atoms[a][-1]+ref_atoms[a][0:-1], ref_amino_acids[a]) in zip(self.atoms, self.amino_acids_number):
                i = zip(self.atoms, self.amino_acids_number).index((ref_atoms[a][-1]+ref_atoms[a][0:-1], ref_amino_acids[a]))
                aux_atoms.append(self.atoms.pop(i))
                aux_pos.append(self.atoms_pos.pop(i))
                aux_aa.append(self.amino_acids_number.pop(i))
                aux_aan.append(self.amino_acids.pop(i))
            elif (ref_atoms[a] == "OXT" and ("OC", ref_amino_acids[a]) in zip(self.atoms, self.amino_acids_number)):
                    i = zip(self.atoms, self.amino_acids_number).index(("OC", ref_amino_acids[a]))
                    aux_atoms.append(self.atoms.pop(i))
                    aux_pos.append(self.atoms_pos.pop(i))
                    aux_aa.append(self.amino_acids_number.pop(i))
                    aux_aan.append(self.amino_acids.pop(i))
            elif (ref_atoms[a] == "OC" and ("OXT", ref_amino_acids[a]) in zip(self.atoms, self.amino_acids_number)):
                    i = zip(self.atoms, self.amino_acids_number).index(("OXT", ref_amino_acids[a]))
                    aux_atoms.append(self.atoms.pop(i))
                    aux_pos.append(self.atoms_pos.pop(i))
                    aux_aa.append(self.amino_acids_number.pop(i))
                    aux_aan.append(self.amino_acids.pop(i))
            elif (ref_atoms[a] == "H" and ("1H", ref_amino_acids[a]) in zip(self.atoms, self.amino_acids_number)):
                    i = zip(self.atoms, self.amino_acids_number).index(("1H", ref_amino_acids[a]))
                    aux_atoms.append(self.atoms.pop(i))
                    aux_pos.append(self.atoms_pos.pop(i))
                    aux_aa.append(self.amino_acids_number.pop(i))
                    aux_aan.append(self.amino_acids.pop(i))
            elif (ref_atoms[a] == "1H" and ("H", ref_amino_acids[a]) in zip(self.atoms, self.amino_acids_number)):
                    i = zip(self.atoms, self.amino_acids_number).index(("H", ref_amino_acids[a]))
                    aux_atoms.append(self.atoms.pop(i))
                    aux_pos.append(self.atoms_pos.pop(i))
                    aux_aa.append(self.amino_acids_number.pop(i))
                    aux_aan.append(self.amino_acids.pop(i))
            elif (len(ref_atoms[a]) == 3 and "H" in ref_atoms[a] and "3" in ref_atoms[a]):
                ra = ref_atoms[a].replace("3", "1")
                ra = ra[-1]+ra[0:-1]
                print("Changing: ")
                print(ref_atoms[a], ra)
                if (ra, ref_amino_acids[a]) in zip(self.atoms, self.amino_acids_number):
                    i = zip(self.atoms, self.amino_acids_number).index((ra, ref_amino_acids[a]))
                    aux_atoms.append(self.atoms.pop(i))
                    aux_pos.append(self.atoms_pos.pop(i))
                    aux_aa.append(self.amino_acids_number.pop(i))
                    aux_aan.append(self.amino_acids.pop(i))
            elif (len(ref_atoms[a]) == 3 and "H" in ref_atoms[a] and "1" in ref_atoms[a]):
                ra = ref_atoms[a].replace("1", "3")
                ra = ra[1:]+ra[0]
                print("Changing: ")
                print(ref_atoms[a], ra)
                if (ra, ref_amino_acids[a]) in zip(self.atoms, self.amino_acids_number):
                    i = zip(self.atoms, self.amino_acids_number).index((ra, ref_amino_acids[a]))
                    aux_atoms.append(self.atoms.pop(i))
                    aux_pos.append(self.atoms_pos.pop(i))
                    aux_aa.append(self.amino_acids_number.pop(i))
                    aux_aan.append(self.amino_acids.pop(i))
            else:
                aux_atoms.append(None)
                aux_pos.append(None)
                aux_aa.append(None)
                aux_aan.append(None)
                print("Not found: ")
                print(ref_atoms[a], ref_amino_acids[a])
        #print(zip(self.atoms, self.amino_acids_number))
        self.atoms = copy.deepcopy(aux_atoms)
        self.atoms_pos = copy.deepcopy(aux_pos)
        self.amino_acids_number = copy.deepcopy(aux_aa)
        self.amino_acids = copy.deepcopy(aux_aan)
        
    def remove_nones(self):
        self.atoms = [x for x in self.atoms if x is not None] 
        self.atoms_pos = [x for x in self.atoms_pos if x is not None] 
        self.amino_acids_number = [x for x in self.amino_acids_number if x is not None]  
        self.amino_acids = [x for x in self.amino_acids if x is not None]   
        
    def is_number(self, s):
        try:
            int(s)
            return True
        except ValueError:
            return False
                
                
       
    
