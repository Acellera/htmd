# (c) 2015-2018 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import os
import logging
from jobqueues.config import config as jqconfig

logger = logging.getLogger(__file__)

_config = {
    "viewer": "VMD",
    "njobs": 1,
    "acemdversion": 3,
    "configfile": os.getenv("HTMD_CONFIG") if os.getenv("HTMD_CONFIG") else None,
    "lsf": None,
    "slurm": None,
}


def config(
    viewer=_config["viewer"],
    ncpus=None,
    njobs=_config["njobs"],
    acemdversion=_config["acemdversion"],
    configfile=_config["configfile"],
    lsf=_config["lsf"],
    slurm=_config["slurm"],
):
    """
    Function to change HTMD configuration variables.

    Parameters
    ----------
    viewer : str
        Defines the backend viewer for molecular visualization
    njobs : int
        Defines the number of parallel jobs spawned for several HTMD operations.
        Negative numbers are used for spawning jobs as many as CPU threads. 
        -1: for all CPUs -2: for all except one etc.
    configfile : str
        Defines the HTMD configuration file that is called at the beginning of importing
    lsf : str
        Defines a YAML file that can contain default profile configurations for an LsfQueue
    slurm : str
        Defines a YAML file that can contain default profile configurations for an SlurmQueue
    """
    _config["viewer"] = viewer
    if ncpus is not None:
        logger.warning(
            "The ncpus config option has been renamed to njobs. Please use njobs instead."
        )
        _config["njobs"] = ncpus
    else:
        _config["njobs"] = njobs
    _config["acemdversion"] = acemdversion
    _config["configfile"] = configfile

    _config["lsf"] = lsf  # DEPRECATE THIS
    _config["slurm"] = slurm  # DEPRECATE THIS
    jqconfig(lsf=lsf, slurm=slurm)
