#Bruno Iochins Grisci
#MAY/2017
#This code gets the atoms positions from a .pdb file

import copy
import numpy as np
import math

class PDB_reader:

    ATOM_TAG = "ATOM"
    END_TAG = "TER"
    ALPHA_CARBON = "CA"
    NITROGEN = "N"
    CARBON = "C"
    BACKBONE_ATOMS = ("N", "CA", "C", "O")
    OC_ATOMS = ("C", "O", "OC", "HOC", "HC", "HO")
    NH_ATOMS = ("N", "1H", "H1", "2H", "H2")
    
    def __init__(self, file_name):
        self.file_name = file_name
        self.atoms = []
        self.atoms_pos = []
        self.amino_acids = []
        self.amino_acids_number = []
        self.more_stuff = []
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
                    extra = map(float, line[pos_init+3:pos_init+5])
                    self.atoms.append(atom)
                    self.atoms_pos.append(pos)
                    self.amino_acids.append(aa)
                    self.amino_acids_number.append(amino_acid)    
                    self.more_stuff.append(extra)              
        pdb.close()
    
    def get_atoms(self):
        return copy.deepcopy(self.atoms)
       
    def get_amino_acids(self):
        return copy.deepcopy(self.amino_acids)       
        
    def get_amino_acids_index(self):
        return copy.deepcopy(self.amino_acids_number)
        
    def get_number_amino_acids(self):
        return len(set(self.amino_acids_number))
        
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
        aux_ms = []
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
                aux_ms.append(self.more_stuff.pop(i))
            elif (ref_atoms[a][1:]+ref_atoms[a][0], ref_amino_acids[a]) in zip(self.atoms, self.amino_acids_number):
                i = zip(self.atoms, self.amino_acids_number).index((ref_atoms[a][1:]+ref_atoms[a][0], ref_amino_acids[a]))
                aux_atoms.append(self.atoms.pop(i))
                aux_pos.append(self.atoms_pos.pop(i))
                aux_aa.append(self.amino_acids_number.pop(i))
                aux_aan.append(self.amino_acids.pop(i))
                aux_ms.append(self.more_stuff.pop(i))
            elif (ref_atoms[a][-1]+ref_atoms[a][0:-1], ref_amino_acids[a]) in zip(self.atoms, self.amino_acids_number):
                i = zip(self.atoms, self.amino_acids_number).index((ref_atoms[a][-1]+ref_atoms[a][0:-1], ref_amino_acids[a]))
                aux_atoms.append(self.atoms.pop(i))
                aux_pos.append(self.atoms_pos.pop(i))
                aux_aa.append(self.amino_acids_number.pop(i))
                aux_aan.append(self.amino_acids.pop(i))
                aux_ms.append(self.more_stuff.pop(i))
            elif (ref_atoms[a] == "OXT" and ("OC", ref_amino_acids[a]) in zip(self.atoms, self.amino_acids_number)):
                    i = zip(self.atoms, self.amino_acids_number).index(("OC", ref_amino_acids[a]))
                    aux_atoms.append(self.atoms.pop(i))
                    aux_pos.append(self.atoms_pos.pop(i))
                    aux_aa.append(self.amino_acids_number.pop(i))
                    aux_aan.append(self.amino_acids.pop(i))
                    aux_ms.append(self.more_stuff.pop(i))
            elif (ref_atoms[a] == "OC" and ("OXT", ref_amino_acids[a]) in zip(self.atoms, self.amino_acids_number)):
                    i = zip(self.atoms, self.amino_acids_number).index(("OXT", ref_amino_acids[a]))
                    aux_atoms.append(self.atoms.pop(i))
                    aux_pos.append(self.atoms_pos.pop(i))
                    aux_aa.append(self.amino_acids_number.pop(i))
                    aux_aan.append(self.amino_acids.pop(i))
                    aux_ms.append(self.more_stuff.pop(i))
            elif (ref_atoms[a] == "H" and ("1H", ref_amino_acids[a]) in zip(self.atoms, self.amino_acids_number)):
                    i = zip(self.atoms, self.amino_acids_number).index(("1H", ref_amino_acids[a]))
                    aux_atoms.append(self.atoms.pop(i))
                    aux_pos.append(self.atoms_pos.pop(i))
                    aux_aa.append(self.amino_acids_number.pop(i))
                    aux_aan.append(self.amino_acids.pop(i))
                    aux_ms.append(self.more_stuff.pop(i))
            elif (ref_atoms[a] == "1H" and ("H", ref_amino_acids[a]) in zip(self.atoms, self.amino_acids_number)):
                    i = zip(self.atoms, self.amino_acids_number).index(("H", ref_amino_acids[a]))
                    aux_atoms.append(self.atoms.pop(i))
                    aux_pos.append(self.atoms_pos.pop(i))
                    aux_aa.append(self.amino_acids_number.pop(i))
                    aux_aan.append(self.amino_acids.pop(i))
                    aux_ms.append(self.more_stuff.pop(i))
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
                    aux_ms.append(self.more_stuff.pop(i))
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
                    aux_ms.append(self.more_stuff.pop(i))
            else:
                aux_atoms.append(None)
                aux_pos.append(None)
                aux_aa.append(None)
                aux_aan.append(None)
                aux_ms.append(None)
                print("Not found: ")
                print(ref_atoms[a], ref_amino_acids[a])
        #print(zip(self.atoms, self.amino_acids_number))
        self.atoms = copy.deepcopy(aux_atoms)
        self.atoms_pos = copy.deepcopy(aux_pos)
        self.amino_acids_number = copy.deepcopy(aux_aa)
        self.amino_acids = copy.deepcopy(aux_aan)
        self.more_stuff = copy.deepcopy(aux_ms)
        
    def remove_nones(self):
        self.atoms = [x for x in self.atoms if x is not None] 
        self.atoms_pos = [x for x in self.atoms_pos if x is not None] 
        self.amino_acids_number = [x for x in self.amino_acids_number if x is not None]  
        self.amino_acids = [x for x in self.amino_acids if x is not None]  
        self.more_stuff = [x for x in self.more_stuff if x is not None] 
        
    def calc_angles(self, pos1, pos2, pos3, pos4):
        #https://math.stackexchange.com/questions/47059/how-do-i-calculate-a-dihedral-angle-given-cartesian-coordinates
        pos1 = np.array(pos1)
        pos2 = np.array(pos2)
        pos3 = np.array(pos3)
        pos4 = np.array(pos4)
        bond12 = pos2 - pos1
        bond23 = pos3 - pos2
        bond34 = pos4 - pos3
        normal1 = np.cross(bond12, bond23) / np.linalg.norm(np.cross(bond12, bond23))  
        normal2 = np.cross(bond23, bond34) / np.linalg.norm(np.cross(bond23, bond34))
        m1 = np.cross(normal1, bond23 / np.linalg.norm(bond23))
        x = np.inner(normal1, normal2)
        y = np.inner(m1, normal2)
        angle = math.atan2(y, x)
        #angle = np.degrees(angle)
        return -angle

    def get_angles(self):
        angles = []
    
        ca = self.get_ca_info()
        n  = self.get_N_info()
        c  = self.get_C_info() 
        
        for aa in xrange(len(ca)):
            name   = ca[aa][1]
            if aa > 0:
                pre_c_pos = c[aa-1][2]
            n_pos  =  n[aa][2]
            ca_pos = ca[aa][2]
            c_pos  =  c[aa][2]
            if aa < len(ca) - 1:
                nex_n_pos = n[aa+1][2]       
            if aa > 0: 
                phi = self.calc_angles(pre_c_pos, n_pos, ca_pos, c_pos)
            else:
                phi = 2 * math.pi   
            if aa < len(ca) - 1: 
                psi = self.calc_angles(n_pos, ca_pos, c_pos, nex_n_pos)
            else:
                psi = 2 * math.pi 
            angles.append(phi)
            angles.append(psi)
        return angles

    def rotate(self, angle):
        rot = self.rotaxis2m(angle, self.atoms_pos[20])
        self.atoms_pos = np.matrix.tolist(np.matrix(self.atoms_pos) * rot.transpose())
        
    def rotate_to(self, angles):
        n_aa = self.get_number_amino_acids()
        for i in xrange(n_aa):     
            #ROTATE PHI
            n_i = zip(self.atoms, self.amino_acids_number).index(("N", i + min(self.amino_acids_number)))   
            ca_i = zip(self.atoms, self.amino_acids_number).index(("CA", i + min(self.amino_acids_number)))
            current_angles = self.get_angles()
            dphi = math.atan2(math.sin(angles[2*i] - current_angles[2*i]), math.cos(angles[2*i] - current_angles[2*i]))
            n_pos = self.atoms_pos[n_i]
            ca_pos = self.atoms_pos[ca_i]                
            ia = 0
            for atom in zip(self.atoms, self.amino_acids_number):
                if i > 0 and (atom[1] > i + min(self.amino_acids_number) or (atom[1] == i + min(self.amino_acids_number) and (atom[0] in self.OC_ATOMS))): 
                    self.atoms_pos[ia] = self.rotate_atom_around_bond(dphi, self.atoms_pos[ia], n_pos, ca_pos)   
                ia += 1        
            #ROTATE PSI    
            c_i  = zip(self.atoms, self.amino_acids_number).index(("C",  i + min(self.amino_acids_number)))  
            ca_i = zip(self.atoms, self.amino_acids_number).index(("CA", i + min(self.amino_acids_number)))
            current_angles = self.get_angles()
            dpsi = math.atan2(math.sin(angles[2*i+1] - current_angles[2*i+1]), math.cos(angles[2*i+1] - current_angles[2*i+1]))              
            c_pos = self.atoms_pos[c_i] 
            ca_pos = self.atoms_pos[ca_i]
            ia = 0
            for atom in zip(self.atoms, self.amino_acids_number):
                if atom[1] > i + min(self.amino_acids_number) and i + min(self.amino_acids_number) < max(self.amino_acids_number): 
                    self.atoms_pos[ia] = self.rotate_atom_around_bond(dpsi, self.atoms_pos[ia], ca_pos, c_pos)          
                ia += 1  
            #self.write_pdb("t" + str(i) + ".pdb")
            
    def normalize(self, v):
        norm = np.linalg.norm(v)
        if norm == 0: 
           return v
        return v/norm  
    
    def rotate_atom_around_bond(self, theta, atom_pos, bond_start, bond_end):
        #https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula
        v = np.array(atom_pos) - np.array(bond_start)
        k = np.array(bond_end) - np.array(bond_start)
        k = self.normalize(k)
        rot_pos = v * np.cos(theta) + (np.cross(k, v)) * np.sin(theta) + k * (np.dot(k,v)) * (1.0 - np.cos(theta))
        return list(rot_pos + np.array(bond_start))
        
    def is_number(self, s):
        try:
            int(s)
            return True
        except ValueError:
            return False
                
    def write_pdb(self, file_name):
        pdb = open(file_name, "w")
        for i in xrange(len(self.atoms)):
            line = "{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}".format(self.ATOM_TAG, i+1, self.atoms[i], " ", self.amino_acids[i], " ", self.amino_acids_number[i], " ", self.atoms_pos[i][0], self.atoms_pos[i][1], self.atoms_pos[i][2], self.more_stuff[i][0], self.more_stuff[i][1], " ", " ")
            line = line + "\n"
            pdb.write(line)
        pdb.write(self.END_TAG)
        pdb.close()                
                
       
    
