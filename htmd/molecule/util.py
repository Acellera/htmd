# (c) 2015-2022 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
from moleculekit.util import *
from moleculekit.dihedral import *
import logging


logger = logging.getLogger(__name__)

logger.warning(
    "Please do not import from htmd.molecule.util. This will be deprecated. "
    "Use moleculekit.util instead or moleculekit.dihedral for the dihedral functionality."
)
