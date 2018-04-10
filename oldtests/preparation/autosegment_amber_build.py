# (c) 2015-2018 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
from htmd import *

p, t = htmd.util.opm("3am6")
pc = p.copy()

p.filter("protein")
p.write("/tmp/p.pdb")

ps = autoSegment(p)
ps.write("/tmp/ps.pdb")

psa = proteinPrepare(ps)
psa.write("/tmp/psa.pdb")

# This fails
# psaa = amber.build(psa, ionize=False, outdir="/tmp/builda")

ps2 = autoSegment(p, field="both")
ps2.write("/tmp/ps2.pdb")

ps2a = proteinPrepare(ps2)
ps2a.write("/tmp/ps2a.pdb")

ps2aa = amber.build(ps2a, ionize=False, outdir="/tmp/builda2")
