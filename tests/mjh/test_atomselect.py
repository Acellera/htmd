# (c) 2015-2017 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
from htmd import *

a =Molecule('filtered.pdb')
print( a.coords.shape )
print(a.box.shape)
print(a.box)
print(a.bonds)
sel1 = a.atomselect("protein")

a.read("filtered.pdb")
print( a.coords.shape )
print(a.box.shape)
print(a.box)
print(a.bonds)
sel2 = a.atomselect("protein")

a.read("traj.xtc" )
print( a.coords.shape )
print(a.box.shape)
print(a.box)
print(a.bonds)
sel3 = a.atomselect("protein")

print(sel1)
print(sel2)
print(sel3)

b=Molecule('filtered.pdb')
b.coords = a.coords[:,:,0].reshape((644, 3, 1))

b.write( 'temp.pdb' )
b =Molecule('temp.pdb' )
print(b.coords.shape)
sel4 = b.atomselect("protein")
print(sel4)
