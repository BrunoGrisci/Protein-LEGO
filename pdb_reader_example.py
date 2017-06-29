#Bruno Iochins Grisci
#JUNE/2017

import sys
import math
from pdb_reader import PDB_reader

pdb_file = sys.argv[1]
pdb = PDB_reader(pdb_file)

naa = pdb.get_number_amino_acids()
print(naa)

ca = pdb.get_ca_info()
print(ca)

n = pdb.get_N_info()
print(n)

c = pdb.get_C_info()
print(c)

print(len(ca), len(n), len(c))

angles = pdb.get_angles()
for a in xrange(0, naa*2, 2):
    print(round(math.degrees(angles[a]),2), round(math.degrees(angles[a+1]),2))
    
pis = [math.radians(180.0)]*len(angles)    
'''pis = [math.radians(360.0),math.radians(176.63),
       math.radians(148.48),math.radians(-21.96),
       math.radians(114.02),math.radians(29.89),
       math.radians(-88.0),math.radians(-38.16),
       math.radians(-74.24),math.radians(360.0)]'''
       
'''pis = [math.radians(360.0),math.radians(-129.17),
       math.radians(67.55),math.radians(172.84),
       math.radians(-59.59),math.radians(-19.44),
       math.radians(-66.36),math.radians(-61.65),
       math.radians(-64.12),math.radians(-17.93),
       math.radians(-78.25),math.radians(-33.68),
       math.radians(-83.98),math.radians(-18.06),
       math.radians(-85.68),math.radians(-40.82),
       math.radians(-75.09),math.radians(-30.75),
       math.radians(-77.64),math.radians(-46.97),
       math.radians(-61.33),math.radians(-27.44),
       math.radians(-60.70),math.radians(-47.46),
       math.radians(-71.15),math.radians(-38.57),
       math.radians(-46.24),math.radians(-50.74),
       math.radians(-69.13),math.radians(-47.41),
       math.radians(-41.89),math.radians(-52.57),
       math.radians(-82.58),math.radians(-23.68),
       math.radians(-53.40),math.radians(-63.44),
       math.radians(-61.23),math.radians(-30.44),
       math.radians(-61.13),math.radians(-32.28),
       math.radians(-80.62),math.radians(-60.14),
       math.radians(-45.89),math.radians(-34.36),
       math.radians(-74.49),math.radians(-47.76),
       math.radians(-83.45),math.radians(10.99),
       math.radians(-134.55),math.radians(360.00)]'''

for a in pdb.get_peptide_bond_angles():
    print(math.degrees(a))

print("example start")
print(len(pdb.get_all_pos()))

pdb.set_peptide_bond_angles()
pdb.set_NH_angles()
pdb.rotate_omegas()
pdb.rotate_to(pis)

print("########")
print("PEPTIDE BOND")
for a in pdb.get_peptide_bond_angles():
    print(round(math.degrees(a),2))

print("########")
print("OMEGA")
for a in pdb.get_omegas():
    print(round(math.degrees(a),2))

print("########")
print("PHI-PSI")
angles = pdb.get_angles()
for a in xrange(0, naa*2, 2):
    print(round(math.degrees(angles[a]),2), round(math.degrees(angles[a+1]),2))

print("example end")
print(len(pdb.get_all_pos()))

pdb.write_pdb("data/test.pdb")
