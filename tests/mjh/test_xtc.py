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
from htmd.molecule.coordreaders import XTCread
import os
import sys
import numpy as np

x1 = XTCread("../mjh/traj.xtc")
x2 = XTCread("traj.xtc")
if (x1.nframes != x2.nframes) or (x1.nframes == 0) or (x2.nframes == 0):
    print("ERROR")
    print("Frame count mismatch (path test)")
    sys.exit(1)

if os.path.exists(".traj.xtc"):
    os.unlink(".traj.xtc")

x1 = XTCread("traj.xtc")
os.unlink(".traj.xtc")
nframes = len(x1.step)

for i in range(nframes):
    print(i)
    x2 = XTCread("traj.xtc", frames=i)

    #    print(x2.coords.shape)
    #    xx = (x1.coords[:,:,i].reshape(( x2.coords +  (1,))))
    #    print(xx.shape)
    if not np.array_equal(np.squeeze(x2.coords), np.squeeze(x1.coords[:, :, i])):
        print(x1.coords[:, :, i])
        print(x2.coords)
        print("ERROR " + str(i))
        sys.exit(1)

os.unlink(".traj.xtc")

# Test 3

os.chdir('/webdata/li198ha8nfoiw90y2/')
full = XTCread(
    filename='datasets/2//filtered/e52s7_e16s4f11/e52s7_e16s4f11-SDOERR_thrombinLig6x2-0-1-RND6631_9.filtered.xtc')
part = XTCread(
    filename='datasets/2//filtered/e52s7_e16s4f11/e52s7_e16s4f11-SDOERR_thrombinLig6x2-0-1-RND6631_9.filtered.xtc',
    frames=[168, 169])

if not np.array_equal(full.coords[:, :, 168], part.coords[:, :, 0]):
    print("ERROR frame 168")
    sys.exit(1)
if not np.array_equal(full.coords[:, :, 169], part.coords[:, :, 1]):
    print("ERROR frame 169")
    sys.exit(1)

sys.exit(0)
