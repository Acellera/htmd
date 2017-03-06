# (c) 2015-2017 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#

# Author: Juan Eiros Zamora (je714@ic.ac.uk) & Gil Ferreira Hoben (gf712@ic.ac.uk) Organisation: Imperial College London
# Modified by: Stefan Doerr

from htmd.apps.localqueue import LocalGPUQueue, _executeMDcommand
from shutil import which
import logging
import os
logger = logging.getLogger(__name__)


class PmemdLocal(LocalGPUQueue):
    """
    Parameters
    ----------
    pmemd: str
        Path to any executable pmemd binary. If None, it will try to detect
        the pmemd.cuda_SPFP binary
    ngpus : int
        Number of GPU devices that PmemdLocal will use. Each simulation will be
        run on a different GPU. PmemdLocal will use the first `ngpus` devices
        of the machine.
    devices : list
        A list of GPU device indexes on which PmemdLocal is allowed to run
        simulations. Mutually exclusive with `ngpus`
    datadir : str
        A folder to which completed simulations will be moved. If None they
        will be written in the input directory.
    """
    def __init__(self, pmemd=None, ngpus=None, devices=None, datadir=None):
        self.system_name = 'Test'

        if not pmemd:
            try:
                # Try to find pmemd.cuda_SPFP in the path
                pmemd = which("pmemd.cuda_SPFP", mode=os.X_OK)
                if pmemd:
                    logger.info("Found pmemd.cuda_SPFP at '" + pmemd + "'")
            except:
                raise NameError("""Cannot find 'pmemd.cuda_SPFP' in the PATH
                                Set its location with the 'pmemd=' named
                                argument""")

        if not os.access(pmemd, os.X_OK):
            raise NameError("pmemd file '" + pmemd + "' not executable")

        super().__init__(jobfunc, (pmemd, datadir), ngpus=ngpus, devices=devices)


def jobfunc(pmemd, datadir, path, gpuid):
    # path and gpuid arguments are provided by default by the localqueue class as last arguments

    # AMBER automatically runs the calculation on the GPU with the
    # most memory even if that GPU is already in use
    # See: http://ambermd.org/gpus/#Running
    # We assume the machine has the exhaustive mode setup
    # So if other GPUs are being used when we launch the adaptive
    # runs, jobs are only run in the free GPUs
    # You set Persistence and Compute Exclusive Modes by running the following as root:
    # $ nvidia-smi -pm 1
    # $ nvidia-smi -c 3

    # need to tell the shell script what engine we are using
    with open(os.path.join(path, 'MD.sh'), 'r') as bash:
        bash_file = bash.read()

    bash_file = bash_file.replace('ENGINE', pmemd)
    # logger.info('BASH script: {}'.format(bash))

    with open(os.path.join(path, 'MD.sh'), 'w') as equil:
        equil.write(bash_file)

    cmd = """cd {} && bash {} > log.txt 2>&1""".format(os.path.normpath(path), 'MD.sh')
    _executeMDcommand(cmd, path, datadir, 'pmemd.cuda', '*.nc')



