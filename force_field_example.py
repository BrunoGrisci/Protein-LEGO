#Bruno Iochins Grisci
#JUNE/2017

import sys

from force_field import Force_field

bonded_file = sys.argv[1]
nonbonded_file = sys.argv[2]
aminoacids_file = sys.argv[3]

ff = Force_field(bonded_file, nonbonded_file, aminoacids_file)
