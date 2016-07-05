# (c) 2015-2016 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#

from htmd.molecule.molecule import Molecule
from htmd.molecule.vmdparser import guessAnglesAndDihedrals
import numpy as np
import scipy.sparse.csgraph as sp
import os
import sys

class SoftDihedral:
  def __init__( self, atoms, left=[], right=[], left_1=[], right_1=[], equiv=[] ):
    self.atoms  = atoms
    self.left   = left
    self.right  = right
    self.left_1   = left_1
    self.right_1  = right_1
    self.equivalents=equiv


def detectSoftDihedrals( mol, equivalent_atoms ):
  bonds = mol.bonds
  natoms= mol.coords.shape[0]
  conn  = np.zeros( (natoms, natoms), dtype=np.bool )

  # Make a connectivity matrix
  #print(natoms)
  for b in bonds:
    conn[b[0],b[1]] = True
    conn[b[1],b[0]] = True

  # Iterate over each of the dihedrals, checking to see which partition the graph
  possible_soft=[]
  for d in mol.dihedrals:
    a0 = d[0]
    a1 = d[1]
    a2 = d[2]
    a3 = d[3]
    c=conn.copy()
    c[a1,a2] = c[a2,a1] = False
    left  =  sp.breadth_first_tree( c, a1, directed=False ).indices.flatten()
    right =  sp.breadth_first_tree( c, a2, directed=False ).indices.flatten()
    left = np.unique(left)
    right= np.unique(right)

    c=conn.copy()
    c[a0,a1] = c[a1,a0] = False
    c[a1,a2] = c[a2,a1] = False
    c[a3,a2] = c[a2,a3] = False
    left_1   =  sp.breadth_first_tree( c, a0, directed=False ).indices.flatten()
    right_1  =  sp.breadth_first_tree( c, a3, directed=False ).indices.flatten()
    left_1   = np.unique(left_1)
    right_1  = np.unique(right_1)
#    print("MARK")
#    print(len(left_1))
#    print(len(right_1))

    if not ( a2 in left)  and not (a1 in right): 
      possible_soft.append( SoftDihedral( d, left, right,  left_1, right_1, [] ) ) 

  final_soft=[]
  e = mol.element
  for d in possible_soft:
      a1   = d.atoms[1]
      a2   = d.atoms[2]
      left = d.left
      right= d.right

      # Exclude trivial dihedrals with just one H atom on a side
   #   if (len(left) == 1)   and (e[left[0]] =='H') : continue
   #   if (len(right) == 1)  and (e[right[0]]=='H') : continue

      # Exclude methyls
      if(len(left) == 3 ):
         if e[a1] == 'C' and e[left[0]]=='H' and e[left[1]]=='H' and e[left[2]]=='H': continue
      if(len(right) == 3 ):
         if e[a1] == 'C' and e[right[0]]=='H' and e[right[1]]=='H' and e[right[2]]=='H': continue
      found=False
      # check to see if the torsional pair of atoms are already included in the list.
      for g in final_soft:
        f=g.atoms
        if f[1] == a1 and f[2] == a2: 
          found=True
          break;
        if f[2] == a1 and f[1] == a2: 
          found=True
          break;
      if not found:  
        final_soft.append(d)

#  if equivalent_atoms:
  final_soft = remove_equivalents( mol, final_soft, equivalent_atoms )

  idx=0
  for t in final_soft:
    print("Dihedral %d: %d-%d-%d-%d" % (idx, t.atoms[0], t.atoms[1], t.atoms[2], t.atoms[3] ) )
    if( len(t.equivalents) ):
       print(" Has equivalent dihedrals through symmetry: " )
       for s in t.equivalents:
         print(" Dihedral %d-%d-%d-%d" % (s[0], s[1], s[2], s[3] ) )
    idx=idx+1

  return final_soft

def remove_equivalents( mol, soft, equiv ):
  final_soft = []
  equivalent_atom_groups   = equiv[0] # list of groups of equivalent atoms
  equivalent_atoms         = equiv[1] # list of equivalent atoms, indexed by atom
  equivalent_group_by_atom = equiv[2] # mapping from atom index to equivalent atom group

#  print("MAPPPING")
#  print(equivalent_group_by_atom )
  # for each of the soft dihedrals, remove any which are equivalent to others through symmetry
  # compare only the middle two atoms since we care about not duplicating X-A-B-X and X-A'-B'-X 
  final_soft=[]
 
  for i in range(len(soft)):
    h1 = equivalent_group_by_atom[ soft[i].atoms[0]]
    h2 = equivalent_group_by_atom[ soft[i].atoms[1]]
    h3 = equivalent_group_by_atom[ soft[i].atoms[2]]
    h4 = equivalent_group_by_atom[ soft[i].atoms[3]]
    found = None
    for j in range(len(final_soft)):
        g1 = equivalent_group_by_atom[ final_soft[j].atoms[0] ]
        g2 = equivalent_group_by_atom[ final_soft[j].atoms[1] ]
        g3 = equivalent_group_by_atom[ final_soft[j].atoms[2] ]
        g4 = equivalent_group_by_atom[ final_soft[j].atoms[3] ]
       
        if   h2==g2 and h3==g3:
            found=j
        if   h3==g2 and h2==g3:
            found=j

    if found != None:    
       # check to see which the two -- the one in already in the final list
       # and the one we've just found -- has more "weight"
       # So that we aren't choosing based on the arbitrary graph traversal ordering
       # Weight is the # of atoms in the left_1 + right_1 groups

       already_in_list = final_soft[found]
       candidate       = soft[i]
#       print(candidate.left_1)
#       print(candidate.right_1)
       if( ( len( already_in_list.left_1 ) + len( already_in_list.right_1 ) ) <  \
           ( len( candidate.left_1 ) + len( candidate.right_1 ) ) ):
         final_soft.remove( already_in_list )
         final_soft.append( candidate ) 
#         print("SWAPPING" )
#         print( soft[i] )
       else:
#         print("DISCARDING")
#        print( soft[i] )
         pass
    else:
#       print("ADDING" )
#       print( soft[i] )
       final_soft.append( soft[i] )


  # now for each of the unique soft dihedrals note the dihedrals that also use the same type
  for s in final_soft:
    a1 = s.atoms[0]
    a2 = s.atoms[1]
    a3 = s.atoms[2]
    a4 = s.atoms[3]
#    print("LOOKING FOR EQUIVALENTS FOR SOFT DIHEDRAL %d-%d-%d-%d" % ( a1,a2,a3,a4 ) )
    h1 = equivalent_group_by_atom[ a1 ]
    h2 = equivalent_group_by_atom[ a2 ]
    h3 = equivalent_group_by_atom[ a3 ]
    h4 = equivalent_group_by_atom[ a4 ]

    for d in mol.dihedrals:
       b1 = d[0] 
       b2 = d[1] 
       b3 = d[2] 
       b4 = d[3] 
#       print("   CHECKING DIHEDRAL %d-%d-%d-%d" % ( b1,b2,b3,b4))

       g1 = equivalent_group_by_atom[ b1 ]
       g2 = equivalent_group_by_atom[ b2 ]
       g3 = equivalent_group_by_atom[ b3 ]
       g4 = equivalent_group_by_atom[ b4 ]

       found=False
        
       if( a1==b1 and a2==b2 and a3==b3 and a4==b4 ): found = True 
       if( a4==b1 and a3==b2 and a2==b3 and a1==b4 ): found = True 
       if found==False: # d is not the soft dihedral itself
         found=False
#         print( "    CHECKING %d-%d-%d-%d vs %d-%d-%d-%d" %( h1,h2,h3,h4,g1,g2,g3,g4))
         if( h1==g1 and h2==g2 and h3==g3 and h4==g4 ): found = True 
         if( h4==g1 and h3==g2 and h2==g3 and h1==g4 ): found = True 
         if found==True: # this dihedral shares a type with the soft dihedral
#           print("APPENDING %d %d %d %d" % ( b1,b2,b3,b4 ) )
           s.equivalents.append( [ b1, b2, b3, b4 ] )

  return final_soft;

if __name__ == "__main__":
  m = Molecule( "ethanol.mol2" )
  (angles, dihedrals) = guessAnglesAndDihedrals( m.bonds )
  
  m.angles   = angles
  m.dihedrals= dihedrals
  sd = detectSoftDihedrals(m)
  print("Soft dihedrals to fit:")
  for t in sd:
    s=t[0]
    print("%d-%d-%d-%d" % ( s[0],s[1], s[2], s[3] ) )
    print("%s-%s-%s-%s" % ( m.name[s[0]], m.name[s[1]], m.name[s[2]], m.name[s[3]] ))
