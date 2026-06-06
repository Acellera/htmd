# (c) 2015-2022 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import os
import logging

logger = logging.getLogger(__file__)

_config = {
    "njobs": 1,
    "configfile": os.getenv("HTMD_CONFIG") if os.getenv("HTMD_CONFIG") else None,
    "lsf": None,
    "slurm": None,
}


def config(
    ncpus: int | None = None,
    njobs: int = _config["njobs"],
    configfile: str | None = _config["configfile"],
    lsf: str | None = _config["lsf"],
    slurm: str | None = _config["slurm"],
) -> None:
    """Change HTMD configuration variables.

    Parameters
    ----------
    ncpus : int, optional
        Deprecated. Use ``njobs`` instead.
    njobs : int, optional
        Number of parallel jobs spawned for several HTMD operations.
        Negative numbers spawn as many jobs as CPU threads: -1 for all CPUs,
        -2 for all except one, etc.
    configfile : str, optional
        HTMD configuration file called at the beginning of importing.
    lsf : str, optional
        YAML file that can contain default profile configurations for an
        LsfQueue.
    slurm : str, optional
        YAML file that can contain default profile configurations for a
        SlurmQueue.
    """
    from jobqueues.config import config as jqconfig

    if ncpus is not None:
        logger.warning(
            "The ncpus config option has been renamed to njobs. Please use njobs instead."
        )
        _config["njobs"] = ncpus
    else:
        _config["njobs"] = njobs
    _config["configfile"] = configfile

    _config["lsf"] = lsf  # DEPRECATE THIS
    _config["slurm"] = slurm  # DEPRECATE THIS
    jqconfig(lsf=lsf, slurm=slurm)
