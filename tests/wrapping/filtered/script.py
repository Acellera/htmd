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
