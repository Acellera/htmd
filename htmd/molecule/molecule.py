# (c) 2015-2022 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
from moleculekit.molecule import Molecule
import logging

logger = logging.getLogger(__name__)

logger.warning(
    "Please do not import htmd.molecule.molecule. It has been moved to the moleculekit package. This import will be deprecated."
)
