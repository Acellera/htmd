# (c) 2015-2017 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#

# https://github.com/Acellera/htmd/issues/200#issuecomment-265410515

from htmd import *

m = Molecule("4UAI")
m.filter("protein and chain A")

prot = m.copy()
x, y = proteinPrepare(prot, returnDetails=True)

x.write("4UAI-out.pdb")
y.data.to_excel("4UAI-out.xlsx")


