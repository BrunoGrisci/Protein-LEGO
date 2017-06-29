#Bruno Iochins Grisci
#JUNE/2017

import sys
import pprint 

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
        self.BONDED = self.read_file(file_bonds)
        self.NONBONDED = self.read_file(file_nonbonds)
        self.AMINOACIDS = self.read_file(file_aminoacids)
        pp = pprint.PrettyPrinter(indent=2)
        pp.pprint(self.AMINOACIDS)
                        
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
