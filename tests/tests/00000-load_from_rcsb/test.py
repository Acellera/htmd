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
#!/usr/bin/env python
import os
from htmd import *
from htmd.molecule import *

DATADIR=os.environ['DATADIR']

a=Molecule('3eko')
a.writePDB( 'test.pdb', 'all' )


