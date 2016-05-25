#!/usr/bin/env python

from htmd.molecule.molecule import Molecule
from htmd.home import home
import numpy as np
from os import path

ref = Molecule(path.join(home(), 'data', 'metricdistance', 'filtered.pdb'))
ref.read(path.join(home(), 'data', 'metricdistance', 'traj.xtc'))
ref.coords = np.atleast_3d(ref.coords[:, :, 0])
ref._guessBonds()

