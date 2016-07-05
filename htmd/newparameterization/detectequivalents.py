# (c) 2015-2016 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#

from htmd.molecule.molecule import Molecule

from htmd.molecule.vmdparser import guessAnglesAndDihedrals
import numpy as np
import scipy.sparse.csgraph as sp
import sys

def _my_breadth_first( mol, con, start, depth=1 ):
  con = con.copy()
  N   = mol.coords.shape[0]

  ret=[]

  search=[ ]
  for j in start:
    for i in range(N):
    
      if( con[j,i]==1 and j!=i):
        #con[j,i] = 0
        con[i,j] = 0
        if i not in search:
          search.append( i )
  if len(search):
    ret.append(search)
    n = _my_breadth_first(mol, con, search, depth+1 ) 
    if len(n):
      for i in n: ret.append(i)

 # top levelreturn. canonicalise the ordering of the elements at each level
  if(depth==1):
    rr=[]
    for r in ret:
       for i in range(len(r)):
         r[i] = mol.element[r[i]]
       r.sort()
       rr.append("!")
       for i in r: rr.append(i)
    ret = rr
  return ret

def detectEquivalents( mol ):
  bonds = mol.bonds
  natoms= mol.coords.shape[0]
  conn  = np.zeros( (natoms, natoms), dtype=np.bool )

  # Make a connectivity matrix
  #print(natoms)

  for b in bonds:
    conn[b[0],b[1]] = True
    conn[b[1],b[0]] = True


  paths=[]
  uniq_paths=[]
  for b in range(natoms):
    a = _my_breadth_first( mol, conn, [b] ) 
    paths.append(a)

  for b in paths:
    found=False
    for c in uniq_paths:
      if c==b: 
        found=True
        break  
    if not found:
      uniq_paths.append(b)
#  print("UNIQUE")
  #print(uniq_paths)

  equiv_groups=[]
  for b in uniq_paths:
    e=[]
    i=0
    for c in paths:
      if b==c: e.append(i)
      i=i+1
    equiv_groups.append(e)
 
#  print(equiv_groups)

  equiv_atoms         = list(range(natoms))
  equiv_group_by_atom = list(range(natoms))
  i=0
  for a in equiv_groups:
    if type(a )==int: a=[a]
    for b in a:
      equiv_atoms[b] = a
      equiv_group_by_atom[b]=i
    i=i+1

  return (equiv_groups, equiv_atoms, equiv_group_by_atom)  


if __name__ == "__main__":
  m = Molecule( "ethanol.mol2" )
 # (angles, dihedrals) = guessAnglesAndDihedrals( m.bonds )
  
  (sd,sp) = detectEquivalents(m)
  print(sd)
  print(sp)

