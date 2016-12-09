# https://github.com/Acellera/htmd/issues/200#issuecomment-265410515

from htmd import *

m = Molecule("4UAI")
m.filter("protein and chain A")

prot = m.copy()
m, pd = proteinPrepare(prot, returnDetails=True, pH=3.0)

m.write("4UAI-out.pdb")
pd.data.to_excel("4UAI-out.xlsx")

pd.data.forced_protonation[pd.data.resid == 62] = "GLN"

nm, npd = pd.reprepare()
