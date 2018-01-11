# (c) 2015-2017 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#

# Gerard's test case

from htmd import *
import pickle
import sys
sys.setrecursionlimit(10000)

#Input and 1st visualization
mol = Molecule('4uai')
mol.filter("protein and chain A")
prep, prepdata = proteinPrepare(mol, pH=3, returnDetails=True)
d = prepdata.data

toChange = (d.resid == int(52)) & (d.chain == str("A"))

print(d.values[51])
d[toChange]
d.loc[toChange, "forced_protonation"] = "ASH"

#I do a reprepare before the proper reprepare
reprep, reprepdata = prepdata.reprepare()
u = reprepdata.data
print(u.values[51])

reprep, reprepdata = prepdata.reprepare()
u = reprepdata.data
print(u.values[51])

