# (c) 2015-2018 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#

# https://github.com/Acellera/htmd/issues/200#issuecomment-265410515

from htmd import *

m = Molecule("4UAI")
m.filter("protein and chain A")

prot = m.copy()
m, pd = proteinPrepare(prot, returnDetails=True, pH=3.0)

m.write("4UAI-out.pdb")
pd.data.to_excel("4UAI-out.xlsx")

pd.data.forced_protonation[pd.data.resid == 63] = "GLU"
nm, npd = pd.reprepare()


# Outstanding bugs: pKa are not carried over
# npd.data.protonation (and nm.resname) columns do not reflect the changes




try:
    pd.data.forced_protonation[pd.data.resid == 63] = "GLN"
    nm, npd = pd.reprepare()
except:
    print("Got the expected failure")
