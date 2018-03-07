# (c) 2015-2017 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import os

_config = {'viewer': 'VMD',
           'ncpus': -2,
           'configfile': os.getenv('HTMD_CONFIG') if os.getenv('HTMD_CONFIG') else None,
           'lsf': None,
           'slurm': None
           }


def config(viewer=_config['viewer'],
           ncpus=_config['ncpus'],
           configfile=_config['configfile'],
           lsf=_config['lsf'],
           slurm=_config['slurm']):
    """
    Function to change HTMD configuration variables.

    Parameters
    ----------
    viewer : str
        Defines the backend viewer for molecular visualization
    ncpus : int
        Defines the number of cpus available for several HTMD operations
    configfile : str
        Defines the HTMD configuration file that is called at the beginning of importing
    lsf : str
        Defines a YAML file that can contain default profile configurations for an LsfQueue
    slurm : str
        Defines a YAML file that can contain default profile configurations for an SlurmQueue
    """
    _config['viewer'] = viewer
    _config['ncpus'] = ncpus
    _config['configfile'] = configfile
    _config['lsf'] = lsf
    _config['slurm'] = slurm
