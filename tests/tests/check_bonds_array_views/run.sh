#!/usr/bin/env python
import sys

from htmd.molecule.molecule import Molecule
from htmd.home import home
import numpy as np
from os import path

ref = Molecule(path.join(home(), 'data', 'metricdistance', 'filtered.pdb'))
ref.read(path.join(home(), 'data', 'metricdistance', 'traj.xtc'))
ref.coords = np.atleast_3d(ref.coords[:, :, 0] )
len1 = len( ref._guessBonds() )
ref.coords = ref.coords[:, :, 0] 
len2 = len( ref._guessBonds() )
ref.coords = np.array( ref.coords, dtype=np.float32 )
len3 = len( ref._guessBonds() )
print(len1)
print(len2)
print(len3)
if(len1!=4562): sys.exit(1)
if(len2!=4562): sys.exit(2)
if(len3!=4562): sys.exit(3)

