"""
HTMD provides a variety of protocols for setting up different types of simulations.
The most basic protocols include equilibration and production runs, both of which support flat-bottom potentials and
constraints for atoms.
"""

# (c) 2015-2022 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import logging

logger = logging.getLogger(__name__)

logger.warning(
    "HTMD protocols are deprecated and will be removed in a future release. "
    "Please refer to the MD tutorials of HTMD for the latest protocols using ACEMD 4."
)
