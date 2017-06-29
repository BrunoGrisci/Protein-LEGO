#Bruno Iochins Grisci
#JUNE/2017

import sys

class Force_field():

    def __init__(self, file_bonds, file_nonbonds, file_aminoacids):
        self.file_bonds = file_bonds
        self.file_nonbonds = file_nonbonds
        self.file_aminoacids = file_aminoacids
        self.BONDED = {}
        self.NONBONDED = {}
        self.AMINOACIDS = {}
        
        self.read_bonds()
        self.read_nonbonds()
        self.read_aminoacids()
        
    def read_bonds(self):
        fbonds = open(self.file_bonds, 'r')                    
        fbonds.close()
        
    def read_nonbonds(self):
        fnonbonds = open(self.file_nonbonds, 'r')
        fnonbonds.close()
        
    def read_aminoacids(self):
        faminoacids = open(self.file_aminoacids, 'r')
        key = ''
        subkey = ''
        for original_line in faminoacids:
            if original_line[0] != ';' and original_line.strip() and original_line not in ['\n', '\r\n']:
                    if ' [' in original_line:
                        subkey = original_line[original_line.find('[')+1:original_line.find(']')].replace(' ', '')
                        if self.AMINOACIDS[key] == []:
                            self.AMINOACIDS[key] = {}
                        self.AMINOACIDS[key][subkey] = []
                    elif '[' in original_line:
                        key = original_line[original_line.find('[')+1:original_line.find(']')].replace(' ', '')
                        subkey = ''
                        self.AMINOACIDS[key] = []                    
                    else:
                        if subkey != '':
                            self.AMINOACIDS[key][subkey].append(map(self.to_number, original_line.split()))
                        else:
                            self.AMINOACIDS[key].append(map(self.to_number, original_line.split()))
        
        for k in self.AMINOACIDS:
            for kk in self.AMINOACIDS[k]:
                print(k, kk)
        faminoacids.close()
        
    def to_number(self, n):
        try:
            float(n)
            return float(n)
        except ValueError:
            return n
        
