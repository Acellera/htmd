# (c) 2015-2017 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import htmd
from htmd import *
mol = Molecule('filtered.pdb')
mol.read('I4100R34/I4100R34-SDOERR_BARNA-0-4-RND6632_9.filtered.xtc')
mol.view()
mol.wrap( )
mol.view()
mol.wrap( wrapsel="chain B")
mol.view()
input()
