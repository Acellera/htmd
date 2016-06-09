#!/usr/bin/env python

# Check to verify that atomselects return nothing, to match VMD
from htmd import *

for pdb in [ "1gmi.pdb", "1hod.pdb", "1w7b.pdb",  "1zec.pdb", "2dhi.pdb" , "2lzp.pdb", "2x72.pdb" ]:
	print(pdb)
	m = Molecule(pdb)
	for v in [ "P", "P1", "P2", "P3", "P4" ]:
		s = m.atomselect( "segid " + v + " and not protein", indexes=True )
		if(len(s)):
			print(v)
			print(s)
			print(m.name[s])
			print(m.resname[s])

