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
    OH_ATOMS = ("OC", "HOC", "HC", "HO")
    WATER_ATOMS = ("2H", "H", "OC", "HOC", "HC", "HO")
    
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
    IFULLATOM = 10
    
    def __init__(self, file_name):
        self.file_name = file_name
        self.atoms = []
        self.full_atoms = []
        self.read_pdb()
        
    def read_pdb(self):
        stop = False
        pdb = open(self.file_name, "r")
        
        while not stop:
            original_line = pdb.readline()
            if not original_line:
                stop = True
            else:
                dic = {}
                line = original_line.split()
                if line[0] == self.END_TAG:
                    stop = True
                elif line[0] == self.ATOM_TAG:     
                    dic[self.IATOMTAG] = original_line[0:4]        
                    dic[self.IATOMINDEX] = int(original_line[6:11])
                    dic[self.IATOM] = original_line[12:16].replace(" ", "")
                    dic[self.IFULLATOM] = original_line[12:16]
                    dic[self.IAMINOACID] = original_line[17:20].replace(" ", "")
                    dic[self.IAMINOACIDINDEX] = int(original_line[22:26])
                    dic[self.IX] = float(original_line[30:38])
                    dic[self.IY] = float(original_line[38:46])
                    dic[self.IZ] = float(original_line[46:54])
                    dic[self.IA] = float(original_line[54:60])
                    dic[self.IB] = float(original_line[60:66])
                    self.atoms.append(dic)                
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
        ohi = None
        for atom in self.atoms:
            if atom[self.IATOM] == "OC":
                ohi = self.atoms.index(atom)
        return self.atoms[ohi][self.IX], self.atoms[ohi][self.IY], self.atoms[ohi][self.IZ]  
 
    def get_pos_AC(self):
        aci = None
        for atom in self.atoms:
            if atom[self.IATOM] == self.ALPHA_CARBON:
                aci = self.atoms.index(atom)
        return self.atoms[aci][self.IX], self.atoms[aci][self.IY], self.atoms[aci][self.IZ]  

    def get_pos_C(self):
        ci = None
        for atom in self.atoms:
            if atom[self.IATOM] == "C":
                ci = self.atoms.index(atom)
        return self.atoms[ci][self.IX], self.atoms[ci][self.IY], self.atoms[ci][self.IZ]  
            
    def remove_H(self):
        hi = None
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
        ni = None
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
            
    def set_atoms_index(self, index):
        i = index
        for atom in self.atoms:
            atom[self.IATOMINDEX] = i
            i += 1
        return i
                
                
       
    
