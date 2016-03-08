# (c) 2015-2016 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
# (c) 2015-2016 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
from htmd import *
import numpy as np
import sys

mol = Molecule('filtered.pdb')

#mol.read('traj.xtc', frames=551)

a1=mol.get('resname', 'protein and name CA')
#array(['MET', 'LYS', 'VAL', 'ILE', 'PHE', 'LEU', 'LYS', 'ASP', 'VAL',
#       'LYS', 'GLY', 'MET', 'GLY', 'LYS', 'LYS', 'GLY', 'GLU', 'ILE',
#       'LYS', 'ASN', 'VAL', 'ALA', 'ASP', 'GLY', 'TYR', 'ALA', 'ASN',
#       'ASN', 'PHE', 'LEU', 'PHE', 'LYS', 'GLN', 'GLY', 'LEU', 'ALA',
#       'ILE', 'GLU', 'ALA'], dtype=object)


mol.read( 'traj.xtc' )
mol.wrap()
for f in range(len(mol.coords)):
	mol.frame=f;
	a2=mol.get('resname', 'protein and name CA' )
#array(['MET', 'LYS', 'PHE', 'LEU', 'LYS', 'ASP', 'VAL', 'LYS', 'GLY',
#       'MET', 'GLY', 'LYS', 'LYS', 'GLY', 'GLU', 'ILE', 'LYS', 'ASN',
#       'VAL', 'ALA', 'ASP', 'GLY', 'TYR', 'ALA', 'ASN', 'ASN', 'PHE',
#       'LEU', 'PHE', 'LYS', 'GLN', 'GLY', 'LEU', 'ALA', 'ILE', 'GLU', 'ALA'], dtype=object)

	if( not np.array_equal(a1 , a2 ) ):
		print("ERROR")
		sys.exit(1)

sys.exit(0)
