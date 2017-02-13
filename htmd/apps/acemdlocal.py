# (c) 2015-2017 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
from htmd.apps.localqueue import LocalGPUQueue, _executeMDcommand
from shutil import which, move
import os
from glob import glob as glob
import logging
logger = logging.getLogger(__name__)


class AcemdLocal(LocalGPUQueue):
    """ Class for running local ACEMD simulations on GPUs. Uses a queuing system.

    Parameters
    ----------
    acemd : str
        Path to ACEMD executable. If None, will try to detect it.
    ngpus : int
        Number of GPU devices that AcemdLocal will use. Each simulation will be run on a different GPU. AcemdLocal will
        use the first `ngpus` devices of the machine.
    devices : list
        A list of GPU device indexes on which AcemdLocal is allowed to run simulations. Mutually exclusive with `ngpus`
    datadir : str
        A folder to which completed simulations will be moved. If None they will be written in the input directory.
    inputfile : str
        The name of the input file for the simulation.
    timeout : float
        Timeout for simulations. Can be used to run simulations only for X amounts of minutes before exiting.

    .. currentmodule:: htmd.apps.acemdlocal.AcemdLocal
    .. rubric:: Methods
    .. autoautosummary:: htmd.apps.acemdlocal.AcemdLocal
        :methods:
        :inherited-members:
    .. rubric:: Attributes
    .. autoautosummary:: htmd.apps.acemdlocal.AcemdLocal
        :attributes:
    """
    def __init__(self, acemd=None, ngpus=None, devices=None, datadir=None, inputfile='input', timeout=None):
        if not acemd:
            try:
                # Try to find acemd in the path
                acemd = which("acemd", mode=os.X_OK)
                if acemd:
                    logger.info("Found ACEMD at '" + acemd + "'")
                else:
                    acemd = which("acemdhtmd", mode=os.X_OK)
                    if acemd:
                        logger.info("Found ACEMD at '" + acemd + "'")
            except:
                pass

        if not acemd:
            raise NameError("Cannot find 'acemd' in the PATH. Set its location with the 'acemd=' named argument")

        if not os.access(acemd, os.X_OK):
            raise NameError("ACEMD file '" + acemd + "' not executable")

        super().__init__(jobfunc, (acemd, datadir, inputfile, timeout), ngpus=ngpus, devices=devices)


def jobfunc(acemd, datadir, inputfile, timeout, path, gpuid):
    # path and gpuid arguments are provided by default by the localqueue class as last arguments
    import os
    timeoutstr = ''
    if timeout:
        timeoutstr = 'timeout {}'.format(timeout)
    cmd = 'cd {}; {} {} --device {} {} > log.txt 2>&1'.format(os.path.normpath(path), timeoutstr, acemd, gpuid, inputfile)
    logger.debug(cmd)
    _executeMDcommand(cmd, path, datadir, 'ACEMD', '*.xtc')


if __name__ == "__main__":
    from time import sleep
    # TODO: Fix this to work
    """
    a = AcemdLocal(acemd="/shared/acemd/bin/acemd", ngpus=2)
    a.submit("/tmp/job/1")
    a.submit("/tmp/job/2")
    a.submit("/tmp/job/3")
    a.submit("/tmp/job/4")
    a.wait()

    r = a.retrieve()
    if len(r):
        print("Retrieved: ")
        print(r)
    sleep(1)
    a.stop()
    sleep(10)
    print("Done")
    """
