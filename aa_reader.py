#Bruno Iochins Grisci
#JUNE/2017
#This code gets the atoms positions from a .pdb file

import copy

class AA_reader:

    ATOM_TAG = "ATOM"
    END_TAG = "TER"
    ALPHA_CARBON = "CA"
    NITROGEN = "N"
    BACKBONE_ATOMS = ("N", "CA", "C", "O")
    H_ATOMS = ("2H", "H")
    OH_ATOMS = ("OC", "HOC", "HC")
    WATER_ATOMS = ("2H", "H", "OC", "HOC", "HC")
    
    IATOMTAG = 0
    IATOMINDEX = 1
    IATOM = 2
    IAMINOACID = 3
    IAMINOACIDINDEX = 4
    IX = 5
    IY = 6
    IZ = 7
    IA = 8
    IB = 9
    
    def __init__(self, file_name):
        self.file_name = file_name
        self.atoms = []
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
                    line[self.IATOMINDEX] = int(line[self.IATOMINDEX])
                    line[self.IAMINOACIDINDEX] = int(line[self.IAMINOACIDINDEX])
                    line[self.IX] = float(line[self.IX])
                    line[self.IY] = float(line[self.IY])
                    line[self.IZ] = float(line[self.IZ])
                    line[self.IA] = float(line[self.IA])
                    line[self.IB] = float(line[self.IB])
                    self.atoms.append(line)                 
        pdb.close()
    
    def get_atoms(self):
        return copy.deepcopy(self.atoms)
        
    def send_origin(self):
        aci = 0
        for atom in self.atoms:
            if atom[self.IATOM] == self.NITROGEN:
                aci = self.atoms.index(atom)
        tx = self.atoms[aci][self.IX] * -1.0
        ty = self.atoms[aci][self.IY] * -1.0
        tz = self.atoms[aci][self.IZ] * -1.0
        
        for atom in self.atoms:
            atom[self.IX] += tx
            atom[self.IY] += ty
            atom[self.IZ] += tz
            
    def get_pos_OC(self):
        ohi = 0
        for atom in self.atoms:
            if atom[self.IATOM] == "OC":
                ohi = self.atoms.index(atom)
        return self.atoms[ohi][self.IX], self.atoms[ohi][self.IY], self.atoms[ohi][self.IZ]  
            
    def remove_H(self):
        hi = 0
        for atom in self.atoms:
            if atom[self.IATOM] in self.H_ATOMS:
                hi = self.atoms.index(atom)
        self.atoms.pop(hi)

    def remove_OH(self):         
        ohi = []
        for atom in self.atoms:
            if atom[self.IATOM] in self.OH_ATOMS:
                ohi.append(self.atoms.index(atom))
        for i in reversed(ohi):
            self.atoms.pop(i) 
            
    def relocate_N(self, (x, y, z)):
        ni = 0
        for atom in self.atoms:
            if atom[self.IATOM] == self.NITROGEN:
                ni = self.atoms.index(atom)
        tx = x - self.atoms[ni][self.IX]
        ty = y - self.atoms[ni][self.IY]
        tz = z - self.atoms[ni][self.IZ]
        for atom in self.atoms:
            atom[self.IX] += tx
            atom[self.IY] += ty
            atom[self.IZ] += tz
            
    def set_amino_acid_index(self, index):
        for atom in self.atoms:
            atom[self.IAMINOACIDINDEX] = index
                
                
       
    
