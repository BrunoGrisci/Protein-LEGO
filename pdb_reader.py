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
    NH_ATOMS = ("N", "H", "1H", "H1", "2H", "H2", "3H", "H3")
    NHC_ATOMS = ("N", "H", "1H", "H1", "2H", "H2", "3H", "H3", "CA")
    
    def __init__(self, file_name):
        self.file_name = file_name
        self.sequence = []
        self.atoms = []
        self.atoms_full = []
        self.atoms_pos = []
        self.amino_acids = []
        self.amino_acids_number = []
        self.more_stuff = []
        self.read_pdb()
        
    def read_pdb(self):
    #https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/pdbintro.html
        stop = False
        pdb = open(self.file_name, "r")       
        while not stop:
            original_line = pdb.readline()
            if not original_line:
                stop = True
            else:
                line = original_line.split()
                if line[0] == self.END_TAG:
                    stop = True
                elif line[0] == self.ATOM_TAG:                    
                    atom_full  = original_line[12:16]
                    atom       = atom_full.replace(" ", "")
                    aa         = original_line[17:20].replace(" ", "")
                    amino_acid = int(original_line[22:26])
                    pos        = map(float, [original_line[30:38], original_line[38:46], original_line[46:54]])
                    extra      = map(float, [original_line[54:60], original_line[60:66]])
                    
                    self.atoms.append(atom)
                    self.atoms_full.append(atom_full)
                    self.atoms_pos.append(pos)
                    self.amino_acids.append(aa)
                    self.amino_acids_number.append(amino_acid)    
                    self.more_stuff.append(extra)              
        pdb.close()
        seq = list(set(zip(self.amino_acids, self.amino_acids_number)))
        seq.sort(key=lambda tup: tup[1])
        for aa in seq:
            self.sequence.append(aa[0])

    def get_sequence(self):
        return copy.deepcopy(self.sequence)
        
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

    def move_to_origin(self):
        first_pos = np.matrix([self.atoms_pos[0]] * len(self.atoms_pos))
        pos = np.matrix(self.atoms_pos)
        self.atoms_pos = np.matrix.tolist(pos - first_pos)
        
    def move_to(self, target):
        self.move_to_origin()
        movement = np.matrix([target] * len(self.atoms_pos))
        pos = np.matrix(self.atoms_pos)
        self.atoms_pos = np.matrix.tolist(pos + movement)

    def translate(self, mov):
        pos = np.matrix(self.atoms_pos)
        movement = np.matrix([mov] * len(self.atoms_pos))
        self.atoms_pos = np.matrix.tolist(pos + movement)
    
    def rotate(self, rot):   
        rotX = np.matrix([[1.0, 0.0, 0.0], [0.0, math.cos(rot[0]), -math.sin(rot[0])], [0.0, math.sin(rot[0]), math.cos(rot[0])]])
        rotY = np.matrix([[math.cos(rot[1]), 0.0, math.sin(rot[1])], [0.0, 1.0, 0.0], [-math.sin(rot[1]), 0.0, math.cos(rot[1])]])
        rotZ = np.matrix([[math.cos(rot[2]), -math.sin(rot[2]), 0.0], [math.sin(rot[2]), math.cos(rot[2]), 0.0], [0.0, 0.0, 1.0]])
        rotXYZ = rotZ * rotY * rotX            
        pos = np.matrix(self.atoms_pos) 
        self.atoms_pos = np.matrix.tolist(pos * rotXYZ.transpose())       
        
    def match_atoms(self, ref_atoms, ref_amino_acids):
        aux_atoms = []
        aux_atomsf = []
        aux_pos = []
        aux_aa = []
        aux_aan = []
        aux_ms = []
        for a in xrange(len(ref_atoms)):
            i = -1    
            if ref_atoms[a] == None:
                pass   
            elif (ref_atoms[a], ref_amino_acids[a]) in zip(self.atoms, self.amino_acids_number):
                i = zip(self.atoms, self.amino_acids_number).index((ref_atoms[a], ref_amino_acids[a]))
            elif (ref_atoms[a][1:]+ref_atoms[a][0], ref_amino_acids[a]) in zip(self.atoms, self.amino_acids_number):
                i = zip(self.atoms, self.amino_acids_number).index((ref_atoms[a][1:]+ref_atoms[a][0], ref_amino_acids[a]))
            elif (ref_atoms[a][-1]+ref_atoms[a][0:-1], ref_amino_acids[a]) in zip(self.atoms, self.amino_acids_number):
                i = zip(self.atoms, self.amino_acids_number).index((ref_atoms[a][-1]+ref_atoms[a][0:-1], ref_amino_acids[a]))
            elif (ref_atoms[a] == "OXT" and ("OC", ref_amino_acids[a]) in zip(self.atoms, self.amino_acids_number)):
                i = zip(self.atoms, self.amino_acids_number).index(("OC", ref_amino_acids[a]))
            elif (ref_atoms[a] == "OC" and ("OXT", ref_amino_acids[a]) in zip(self.atoms, self.amino_acids_number)):
                i = zip(self.atoms, self.amino_acids_number).index(("OXT", ref_amino_acids[a]))
            elif (ref_atoms[a] == "H" and ("1H", ref_amino_acids[a]) in zip(self.atoms, self.amino_acids_number)):
                i = zip(self.atoms, self.amino_acids_number).index(("1H", ref_amino_acids[a]))
            elif (ref_atoms[a] == "1H" and ("H", ref_amino_acids[a]) in zip(self.atoms, self.amino_acids_number)):
                i = zip(self.atoms, self.amino_acids_number).index(("H", ref_amino_acids[a]))
            elif (len(ref_atoms[a]) == 3 and "H" in ref_atoms[a] and "3" in ref_atoms[a]):
                ra = ref_atoms[a].replace("3", "1")
                ra = ra[-1]+ra[0:-1]
                if (ra, ref_amino_acids[a]) in zip(self.atoms, self.amino_acids_number):
                    i = zip(self.atoms, self.amino_acids_number).index((ra, ref_amino_acids[a]))
            elif (len(ref_atoms[a]) == 3 and "H" in ref_atoms[a] and "1" in ref_atoms[a]):
                ra = ref_atoms[a].replace("1", "3")
                ra = ra[1:]+ra[0]
                if (ra, ref_amino_acids[a]) in zip(self.atoms, self.amino_acids_number):
                    i = zip(self.atoms, self.amino_acids_number).index((ra, ref_amino_acids[a]))
            else:
                aux_atoms.append(None)
                aux_atomsf.append(None)
                aux_pos.append(None)
                aux_aa.append(None)
                aux_aan.append(None)
                aux_ms.append(None)
            if i >= 0:
                aux_atoms.append(self.atoms.pop(i))
                aux_atomsf.append(self.atoms_full.pop(i))
                aux_pos.append(self.atoms_pos.pop(i))
                aux_aa.append(self.amino_acids_number.pop(i))
                aux_aan.append(self.amino_acids.pop(i))
                aux_ms.append(self.more_stuff.pop(i))                
        self.atoms = copy.deepcopy(aux_atoms)
        self.atoms_full = copy.deepcopy(aux_atomsf)
        self.atoms_pos = copy.deepcopy(aux_pos)
        self.amino_acids_number = copy.deepcopy(aux_aa)
        self.amino_acids = copy.deepcopy(aux_aan)
        self.more_stuff = copy.deepcopy(aux_ms)
        
    def remove_nones(self):
        self.atoms = [x for x in self.atoms if x is not None]
        self.atoms_full = [x for x in self.atoms_full if x is not None] 
        self.atoms_pos = [x for x in self.atoms_pos if x is not None] 
        self.amino_acids_number = [x for x in self.amino_acids_number if x is not None]  
        self.amino_acids = [x for x in self.amino_acids if x is not None]  
        self.more_stuff = [x for x in self.more_stuff if x is not None]

    def angle_diff(self, angle_source, angle_target):
        a = angle_target - angle_source
        if a > math.pi:
            a -= math.pi * 2.0
        if a < -math.pi:
            a += math.pi * 2.0
        return a
        
    def calc_angle_3(self, pos1, posC, pos2):
        #https://stackoverflow.com/questions/19729831/angle-between-3-points-in-3d-space
        '''In pseudo-code, the vector BA (call it v1) is:
        v1 = {A.x - B.x, A.y - B.y, A.z - B.z}
        Similarly the vector BC (call it v2) is:
        v2 = {C.x - B.x, C.y - B.y, C.z - B.z}
        The dot product of v1 and v2 is a function of the cosine of the angle between them (it's scaled by the product of their magnitudes). So first normalize v1 and v2:
        v1mag = sqrt(v1.x * v1.x + v1.y * v1.y + v1.z * v1.z)
        v1norm = {v1.x / v1mag, v1.y / v1mag, v1.z / v1mag}
        v2mag = sqrt(v2.x * v2.x + v2.y * v2.y + v2.z * v2.z)
        v2norm = {v2.x / v2mag, v2.y / v2mag, v2.z / v2mag}
        Then calculate the dot product:
        res = v1norm.x * v2norm.x + v1norm.y * v2norm.y + v1norm.z * v2norm.z
        And finally, recover the angle:
        angle = acos(res)'''
        pos1 = np.array(pos1)
        posC = np.array(posC)
        pos2 = np.array(pos2)
        bond1C = self.normalize(pos1 - posC)
        bond2C = self.normalize(pos2 - posC)
        dp = np.dot(bond1C, bond2C)
        angle = np.arccos(dp)
        return angle  

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
        return -angle

    def get_peptide_bond_angles(self):
        angles = []
        ca = self.get_ca_info()
        n  = self.get_N_info()
        c  = self.get_C_info()
        for aa in xrange(len(ca)):
            if aa < len(ca) - 1:
                name       = ca[aa][1]
                c_pos      =  c[aa][2]                
                nex_n_pos  =  n[aa+1][2]
                nex_ca_pos = ca[aa+1][2]
                alpha = self.calc_angle_3(c_pos, nex_n_pos, nex_ca_pos)
                angles.append(alpha)
            else:
                angles.append(2.0*math.pi)
        return angles 

    def get_peptide_bond_H_angles(self):
        angles = []
        ca = self.get_ca_info()
        n  = self.get_N_info()
        c  = self.get_C_info()
        for aa in xrange(len(ca)):
            if aa < len(ca) - 1:
                name       = ca[aa][1]
                c_pos      =  c[aa][2]                
                nex_n_pos  =  n[aa+1][2]
                nex_ca_pos = ca[aa+1][2]
                alpha = self.calc_angle_3(c_pos, nex_n_pos, nex_ca_pos)
                angles.append(alpha)
            else:
                angles.append(2.0*math.pi)
        return angles 
        
    def set_peptide_bond_angles(self, angles=[]):
        n_aa = self.get_number_amino_acids()
        if len(angles) == 0:
            angles = [math.radians(120.0)]*n_aa
        for i in xrange(n_aa):
            if i + min(self.amino_acids_number) < max(self.amino_acids_number):
                #ROTATE ALPHA
                c_i   = zip(self.atoms, self.amino_acids_number).index(("C",  i + min(self.amino_acids_number)))   # C from aminoacid i
                nn_i  = zip(self.atoms, self.amino_acids_number).index(("N",  i+1 + min(self.amino_acids_number))) # N from aminoacid i+1
                nca_i = zip(self.atoms, self.amino_acids_number).index(("CA", i+1 + min(self.amino_acids_number))) #CA from aminoacid i+1                                   
                c_pos   = self.atoms_pos[c_i]
                nn_pos  = self.atoms_pos[nn_i]
                nca_pos = self.atoms_pos[nca_i]
                current_alpha = self.calc_angle_3(c_pos, nn_pos, nca_pos)
                dalpha = self.angle_diff(current_alpha, angles[i])
                if (current_alpha < dalpha):
                    dalpha = self.angle_diff(current_alpha, -angles[i])
                ia = 0
                for atom in zip(self.atoms, self.amino_acids_number):
                    if (atom[1] > i+1 + min(self.amino_acids_number) or (atom[1] == i+1 + min(self.amino_acids_number) and (atom[0] not in self.NH_ATOMS))): 
                        self.atoms_pos[ia] = self.bend_bonds(dalpha, self.atoms_pos[ia], c_pos, nn_pos, nca_pos)  
                    ia += 1                   

    def set_NH_angles(self, angles=[]):
        n_aa = self.get_number_amino_acids()
        if len(angles) == 0:
            angles = [math.radians(120.0)]*n_aa
        for i in xrange(n_aa):
            if i + min(self.amino_acids_number) < max(self.amino_acids_number):
                #ROTATE ALPHA
                c_i   = zip(self.atoms, self.amino_acids_number).index(("C",  i + min(self.amino_acids_number)))   # C from aminoacid i
                nn_i  = zip(self.atoms, self.amino_acids_number).index(("N",  i+1 + min(self.amino_acids_number))) # N from aminoacid i+1
                nca_i = zip(self.atoms, self.amino_acids_number).index(("CA", i+1 + min(self.amino_acids_number))) #CA from aminoacid i+1
                nh_i  = -1.0
                for a in self.NH_ATOMS:
                    if a != "N":                                                 
                        nh_i = zip(self.atoms, self.amino_acids_number).index((a, i+1+min(self.amino_acids_number))) if (a, i+1+min(self.amino_acids_number)) in zip(self.atoms, self.amino_acids_number) else -1
                        if nh_i >= 0:
                            break                                   
                c_pos   = self.atoms_pos[c_i]
                nn_pos  = self.atoms_pos[nn_i]
                nca_pos = self.atoms_pos[nca_i]
                if nh_i >= 0:
                    nh_pos  = self.atoms_pos[nh_i]
                    current_alphaH = self.calc_angle_3(c_pos, nn_pos, nh_pos)
                    dalphaH = self.angle_diff(current_alphaH, angles[i])
                    if (current_alphaH < dalphaH):
                        dalphaH = self.angle_diff(current_alphaH, -angles[i])
                ia = 0
                for atom in zip(self.atoms, self.amino_acids_number):
                    if nh_i >= 0 and atom[1] == i+1 + min(self.amino_acids_number) and atom[0] in self.NH_ATOMS and atom[0] != "N":
                        self.atoms_pos[ia] = self.bend_bonds(dalphaH, self.atoms_pos[ia], c_pos, nn_pos, nh_pos)  
                    ia += 1

    def bend_bonds(self, theta, atom_pos, pos1, posC, pos2):
        #https://pt.stackoverflow.com/questions/25923/vetores-e-%C3%82ngulos-geometria-molecular
        posC = np.array(posC)
        atom_pos = np.array(atom_pos) - posC
        pos1 = np.array(pos1) - posC
        pos2 = np.array(pos2) - posC
        bond1C = self.normalize(pos1)
        bond2C = self.normalize(pos2)
        ortho = self.normalize(np.cross(bond1C, bond2C))
        c = np.cos(theta)
        s = np.sin(theta)
        t = 1.0 - c
        rot = np.matrix([[c + ortho[0] * ortho[0] * t, ortho[0] * ortho[1] * t - ortho[2] * s, ortho[0] * ortho[2] * t + ortho[1] * s], 
                         [ortho[0] * ortho[1] * t + ortho[2] * s, c + ortho[1] * ortho[1] * t, ortho[1] * ortho[2] * t - ortho[0] * s],
                         [ortho[2] * ortho[0] * t - ortho[1] * s, ortho[2] * ortho[1] * t + ortho[0] * s, c + ortho[2] * ortho[2] * t]])      
        new_pos = np.matrix.tolist(np.matrix(atom_pos) * rot.transpose())[0]               
        new_pos = list(np.array(new_pos) + posC)
        return new_pos
        
    def get_omegas(self):
        angles = []        
        ca = self.get_ca_info()
        n  = self.get_N_info()
        c  = self.get_C_info()       
        for aa in xrange(len(ca)):
            if aa > 0:
                name    = ca[aa][1]
                pca_pos = ca[aa-1][2]
                pc_pos  =  c[aa-1][2]                
                n_pos   =  n[aa][2]
                ca_pos  = ca[aa][2]
                omega   = self.calc_angles(pca_pos, pc_pos, n_pos, ca_pos)
                angles.append(omega)
            else:
                angles.append(2.0*math.pi)    
        return angles        

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

    def rotate_omegas(self, angles=[]):
        n_aa = self.get_number_amino_acids()
        if len(angles) == 0:
            angles = [math.pi]*n_aa
        for i in xrange(n_aa):
            if i > 0:
                #ROTATE OMEGA
                pca_i  = zip(self.atoms, self.amino_acids_number).index(("CA",  i-1 + min(self.amino_acids_number))) #CA from aminoacid i-1
                pc_i   = zip(self.atoms, self.amino_acids_number).index(("C",   i-1 + min(self.amino_acids_number))) #C from aminoacid i-1
                n_i    = zip(self.atoms, self.amino_acids_number).index(("N",   i + min(self.amino_acids_number))) #N from aminoacid i
                ca_i   = zip(self.atoms, self.amino_acids_number).index(("CA",  i + min(self.amino_acids_number))) #CA from aminoacid i
                pca_pos = self.atoms_pos[pca_i]
                pc_pos  = self.atoms_pos[pc_i]
                n_pos   = self.atoms_pos[n_i]
                ca_pos  = self.atoms_pos[ca_i]
                current_omega = self.calc_angles(pca_pos, pc_pos, n_pos, ca_pos)
                domega = self.angle_diff(current_omega, angles[i])
                ia = 0
                for atom in zip(self.atoms, self.amino_acids_number):
                    if (atom[1] > i + min(self.amino_acids_number) or (atom[1] == i + min(self.amino_acids_number) and (atom[0] != "N"))): 
                        self.atoms_pos[ia] = self.rotate_atom_around_bond(domega, self.atoms_pos[ia], pc_pos, n_pos)
                    ia += 1   
        
    def rotate_to(self, angles):
        n_aa = self.get_number_amino_acids()
        for i in xrange(n_aa):     
            #ROTATE PHI
            if i > 0:
                pc_i  = zip(self.atoms, self.amino_acids_number).index(("C",  i-1 + min(self.amino_acids_number)))
                n_i   = zip(self.atoms, self.amino_acids_number).index(("N",  i + min(self.amino_acids_number)))   
                ca_i  = zip(self.atoms, self.amino_acids_number).index(("CA", i + min(self.amino_acids_number)))
                c_i   = zip(self.atoms, self.amino_acids_number).index(("C",  i + min(self.amino_acids_number)))
                pc_pos = self.atoms_pos[pc_i]
                n_pos  = self.atoms_pos[n_i]
                ca_pos = self.atoms_pos[ca_i]
                c_pos  = self.atoms_pos[c_i] 
                current_angle = self.calc_angles(pc_pos, n_pos, ca_pos, c_pos)   
                dphi = self.angle_diff(current_angle, angles[2*i])
                ia = 0
                for atom in zip(self.atoms, self.amino_acids_number):
                    if (atom[1] > i + min(self.amino_acids_number) or (atom[1] == i + min(self.amino_acids_number) and (atom[0] not in self.NHC_ATOMS))): 
                        self.atoms_pos[ia] = self.rotate_atom_around_bond(dphi, self.atoms_pos[ia], n_pos, ca_pos)   
                    ia += 1        
            #ROTATE PSI 
            if i + min(self.amino_acids_number) < max(self.amino_acids_number):
                n_i  = zip(self.atoms, self.amino_acids_number).index(("N",  i + min(self.amino_acids_number)))
                ca_i = zip(self.atoms, self.amino_acids_number).index(("CA", i + min(self.amino_acids_number)))
                c_i  = zip(self.atoms, self.amino_acids_number).index(("C",  i + min(self.amino_acids_number)))
                nn_i = zip(self.atoms, self.amino_acids_number).index(("N",  i+1 + min(self.amino_acids_number)))  
                n_pos  = self.atoms_pos[n_i]
                ca_pos = self.atoms_pos[ca_i]
                c_pos  = self.atoms_pos[c_i]
                nn_pos = self.atoms_pos[nn_i] 
                current_angle = self.calc_angles(n_pos, ca_pos, c_pos, nn_pos)                            
                dpsi = self.angle_diff(current_angle, angles[2*i+1])
                ia = 0
                for atom in zip(self.atoms, self.amino_acids_number):
                    if (atom[1] > i+min(self.amino_acids_number) or (atom[1] == i+min(self.amino_acids_number) and (atom[0]=="O"))): 
                        self.atoms_pos[ia] = self.rotate_atom_around_bond(dpsi, self.atoms_pos[ia], ca_pos, c_pos)         
                    ia += 1  
            
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
                
    def write_pdb(self, file_name):
        pdb = open(file_name, "w")
        for i in xrange(len(self.atoms)):
            line = "{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}".format(self.ATOM_TAG, i+1, self.atoms_full[i], " ", self.amino_acids[i], " ", self.amino_acids_number[i], " ", self.atoms_pos[i][0], self.atoms_pos[i][1], self.atoms_pos[i][2], self.more_stuff[i][0], self.more_stuff[i][1], " ", " ")
            line = line + "\n"
            pdb.write(line)
        pdb.write(self.END_TAG)
        pdb.close()               
                  
