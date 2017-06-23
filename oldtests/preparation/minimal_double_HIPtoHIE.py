# (c) 2015-2017 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
from htmd import *
import pickle
import sys
sys.setrecursionlimit(10000)

#Input and 1st visualization
mol = Molecule('4uai')
mol.filter("protein")
prep, prepdata = proteinPrepare(mol, pH=3, returnDetails=True)
d = prepdata.data
print(d.values[24])
d.loc[(d.resid == int(25)) & (d.chain == str("A")), "forced_protonation"] = "HIE"
#I do a reprepare before the proper reprepare
prepdata.reprepare()
reprep, reprepdata = prepdata.reprepare()
u = reprepdata.data
print(u.values[24])

