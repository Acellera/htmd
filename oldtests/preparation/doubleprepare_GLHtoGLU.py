# (c) 2015-2018 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#

# Gerard's test case

from htmd import *
import pickle
import sys
import numpy as np
sys.setrecursionlimit(10000)

#Input and 1st visualization
mol = Molecule('4uai')
mol.filter("protein and chain A")
prep, prepdata = proteinPrepare(mol, pH=3, returnDetails=True)
d = prepdata.data

toChange = (d.resid == int(60)) & (d.chain == str("A"))

print("Before reprepare: ",d.values[59])
print("GLU 50 has atoms... ",np.sum(prep.atomselect("resid 60 and chain A")))

d[toChange]
d.loc[toChange, "forced_protonation"] = "GLU"
# or, faster -- d.set_value(toChange,"forced_protonation","ASH")

# prepdata.reprepare()

reprep, reprepdata = prepdata.reprepare()
u = reprepdata.data
print("After single reprepare: ",u.values[59])
print(u.iloc[59])
print("GLU 50 has atoms... ",np.sum(reprep.atomselect("resid 60 and chain A")))
#print(u.loc[toChange,"forced_protonation"])




#Input and 1st visualization
mol = Molecule('4uai')
mol.filter("protein and chain A")
prep, prepdata = proteinPrepare(mol, pH=3, returnDetails=True)
d = prepdata.data

toChange = (d.resid == int(60)) & (d.chain == str("A"))

print("Before reprepare: ",d.values[59])
print("GLU 50 has atoms... ",np.sum(prep.atomselect("resid 60 and chain A")))

d[toChange]
d.loc[toChange, "forced_protonation"] = "GLU"
# or, faster -- d.set_value(toChange,"forced_protonation","ASH")

prepdata.reprepare()

reprep, reprepdata = prepdata.reprepare()
u = reprepdata.data
print("After double reprepare: ",u.values[59])
print(u.iloc[59])
print("GLU 50 has atoms... ",np.sum(reprep.atomselect("resid 60 and chain A")))
#print(u.loc[toChange,"forced_protonation"])


