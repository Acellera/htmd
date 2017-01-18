# Gerard's test case

from htmd import *
import pickle
import sys
import numpy as np
sys.setrecursionlimit(10000)

#Input and 1st visualization
mol = Molecule('3ptb')
mol.filter("protein and chain A")

pmol, pdata = proteinPrepare(mol, pH=3.4, returnDetails=True)
d = pdata.data

toChange = (d.resid == int(80)) & (d.chain == str("A"))

print("Before reprepare: ",d.loc[toChange])
print("Residue has this many atoms... ", np.sum(pmol.atomselect("resid 80 and chain A")))

d[toChange]
d.loc[toChange, "forced_protonation"] = "GLH"

pmol2, pdata2 = pdata.reprepare()
u = pdata2.data

print("After reprepare: ",u.loc[toChange])
print("Residue has this many atoms... ", np.sum(pmol2.atomselect("resid 80 and chain A")))

u[toChange]




