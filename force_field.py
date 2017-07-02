#Bruno Iochins Grisci
#JUNE/2017

import sys
import pprint 
import math

import pdb_reader

class Force_field():
    #ffbonded
    BONDTYPES = {'I':   0,
                 'J':   1,
                 'FUNC':2,
                 'B0':  3,
                 'KB':  4}
    CONSTRAINTTYPES = {}             
    ANGLETYPES = {'I':    0,
                  'J':    1,
                  'K':    2,
                  'FUNC': 3,
                  'TH0':  4,
                  'CTH':  5}
    DIHEDRALTYPES = {'I':    0,
                     'J':    1,
                     'K':    2,
                     'L':    3,
                     'FUNC': 4,
                     'PHASE':5,
                     'KD':   6,
                     'PN':   7}
    ##ffnonbonded             
    ATOMTYPES = {'NAME':   0,
                 'ATNUM':  1,
                 'MASS':   2,
                 'CHARGE': 3,
                 'PTYPE':  4,
                 'SIGMA':  5,
                 'EPSILON':6}
                 
    def __init__(self, file_bonds, file_nonbonds, file_aminoacids):
        self.PDB = None
        self.BONDED = self.read_file(file_bonds)
        self.NONBONDED = self.read_file(file_nonbonds)
        self.AMINOACIDS = self.read_file(file_aminoacids)
        self.LJPARAMETERS = {}
        self.build_LJ_parameters()
        #pp = pprint.PrettyPrinter(indent=4)
        #pp.pprint(self.LJPARAMETERS)
        
    def insert_PDB(self, pdb):
        self.PDB = pdb        
                        
    def read_file(self, file_name):
        dic = {}
        key = ''
        subkey = ''        
        f = open(file_name, 'r')
        for original_line in f:
            if original_line[0] != ';' and original_line.strip() and original_line not in ['\n', '\r\n']:
                    if ' [' in original_line:
                        subkey = original_line[original_line.find('[')+1:original_line.find(']')].replace(' ', '')
                        if dic[key] == []:
                            dic[key] = {}
                        dic[key][subkey] = []
                    elif '[' in original_line:
                        key = original_line[original_line.find('[')+1:original_line.find(']')].replace(' ', '')
                        subkey = ''
                        dic[key] = []                    
                    else:
                        if subkey != '':
                            info = original_line.split(';')[0]
                            dic[key][subkey].append(map(self.to_number, info.split()))
                        else:
                            info = original_line.split(';')[0]
                            dic[key].append(map(self.to_number, info.split()))
        f.close()
        return dic
        
    def to_number(self, n):
        try:
            float(n)
            return float(n)
        except ValueError:
            return n    
            
    def build_LJ_parameters(self):
        cons = math.pow(2.0, (1.0/6.0))
        for i in self.NONBONDED['atomtypes']:
            ai_name  = i[self.ATOMTYPES['NAME']]
            ai_sigma = i[self.ATOMTYPES['SIGMA']]
            ai_eps   = i[self.ATOMTYPES['EPSILON']]
            ai_rmin  = ai_sigma * cons
            dic_B = {}
            for j in self.NONBONDED['atomtypes']:
                aj_name  = j[self.ATOMTYPES['NAME']]
                aj_sigma = j[self.ATOMTYPES['SIGMA']]
                aj_eps   = j[self.ATOMTYPES['EPSILON']]            
                aj_rmin  = aj_sigma * cons
                aij_rmin = (ai_rmin + aj_rmin)/2.0
                aij_eps  = math.sqrt(ai_eps * aj_eps)
                dic_A = {'rmin': aij_rmin, 'epsilon': aij_eps}
                dic_B[aj_name] = dic_A
            self.LJPARAMETERS[ai_name] = dic_B

    def are_bonded(self, atom1, atom2, aa1, aa2, aai1, aai2):
        if abs(aai1-aai2) > 1:
            return False
        elif aai1-aai2 == 1 and atom2 == 'C':
            atom2 = '-C'
            bonds = self.AMINOACIDS[aa1]['bonds']
            return [atom1, atom2] in bonds or [atom2, atom1] in bonds
        elif aai2-aai1 == 1 and atom1 == 'C':
            atom1 = '-C'
            bonds = self.AMINOACIDS[aa2]['bonds']
            return [atom1, atom2] in bonds or [atom2, atom1] in bonds            
        elif abs(aai1-aai2) == 0:
            bonds = self.AMINOACIDS[aa1]['bonds']
            return [atom1, atom2] in bonds or [atom2, atom1] in bonds        
 
    def correct_atom_name(self, atom, amino_acid):
        for a in self.AMINOACIDS[amino_acid]['atoms']:
            if atom == a[0]:
                return atom
            elif atom[1:] + atom[0]      == a[0]:
                return atom[1:]+atom[0]
            elif atom == 'OXT' and 'OC2' == a[0]:
                return 'OC2'
            elif atom == 'OC'  and 'OC1' == a[0]:
                return 'OC1'
            elif atom == 'O'  and 'OC1' == a[0]:
                return 'OC1'                
            elif atom == 'H'   and 'H1'  == a[0]:
                return 'H1'
            elif atom == '1H'   and 'H'  == a[0]:
                return 'H'                
            elif (len(atom) == 3 and "H" in atom and "3" in atom):
                ra = atom.replace("3", "1")                
                if ra == a[0]:
                    return ra
                else:
                    if ra[1:] + ra[0] == a[0]:
                        return ra[1:] + ra[0]

    def non_bonded(self):
        lj = 0.0
        coulomb = 0.0
        e1 = 8.99 * math.pow(10.0, 9.0)
        atoms = self.PDB.get_atoms()
        atoms_pos = self.PDB.get_all_pos()
        amino_acids = self.PDB.get_amino_acids()
        aa_index = self.PDB.get_amino_acids_index()
        aa_0 = min(aa_index)
        aa_N = max(aa_index)
        for i in xrange(len(atoms)):
            for j in xrange(len(atoms)):
                if i > j:
                    aai = amino_acids[i]
                    aaj = amino_acids[j]
                    if aa_index[i] == aa_0:
                        aai = 'N' + aai
                    if aa_index[j] == aa_0:
                        aaj = 'N' + aaj
                    if aa_index[i] == aa_N:
                        aai = 'C' + aai
                    if aa_index[j] == aa_N:
                        aaj = 'C' + aaj                   
                    atomi = self.correct_atom_name(atoms[i], aai)
                    atomj = self.correct_atom_name(atoms[j], aaj)  
                    if(atomi is None):
                        print(atoms[i], atomi, aai)
                    if(atomj is None):
                        print(atoms[j], atomj, aaj)                                                                
                    if not self.are_bonded(atomi, atomj, aai, aaj, aa_index[i], aa_index[j]) and atomi is not None and atomj is not None:  
                        atomi_key = ''
                        qi = 0.0    
                        for a in self.AMINOACIDS[aai]['atoms']:
                            if atomi == a[0]:
                                atomi_key = a[1]
                                qi = a[2]                                                   
                        atomj_key = ''
                        qj = 0.0
                        for a in self.AMINOACIDS[aaj]['atoms']:
                            if atomj == a[0]:
                                atomj_key = a[1]
                                qj = a[2] 
                        if(atomj_key == ''):
                            print(atoms[j], atomj, aaj)                                                                         
                        rij  = math.sqrt((atoms_pos[i][0] - atoms_pos[j][0])**2 + (atoms_pos[i][1] - atoms_pos[j][1])**2 + (atoms_pos[i][2] - atoms_pos[j][2])**2)/10.0 #a to nm                         
                        eps  = self.LJPARAMETERS[atomi_key][atomj_key]['epsilon']  
                        rmin = self.LJPARAMETERS[atomi_key][atomj_key]['rmin']  
                        lj      += eps * (math.pow(rmin/rij, 12.0) - 2.0 * math.pow(rmin/rij, 6.0)) 
                        coulomb += (qi * qj) / (e1 * rij)
        return lj + coulomb            
                       
