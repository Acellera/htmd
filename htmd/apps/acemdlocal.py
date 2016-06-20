# (c) 2015-2016 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
from htmd.apps.localqueue import LocalGPUQueue
from shutil import which, move
import os
from glob import glob as glob
import logging
logger = logging.getLogger(__name__)


class AcemdLocal(LocalGPUQueue):
    """
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
    from subprocess import PIPE, Popen, TimeoutExpired, CalledProcessError
    import os

    timeoutstr = ''
    if timeout:
        timeoutstr = 'timeout {}'.format(timeout)
    cmd = 'cd {}; {} {} --device {} {} > log.txt 2>&1'.format(os.path.normpath(path), timeoutstr, acemd, gpuid, inputfile)
    try:
        #check_output(cmd, shell=True, timeout=obj.timeout)
        proc = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
        proc.communicate()
    except CalledProcessError:
        logger.error('Error in ACEMD for path: {}. Check the {} file.'.format(path, os.path.join(path, 'log.txt')))
        proc.kill()
        return
    except TimeoutExpired:
        proc.kill()
        return

    # If a datadir is provided, copy finished trajectories there. Only works for xtc files.
    if datadir is not None:
        if not os.path.isdir(datadir):
            os.mkdir(datadir)
        simname = os.path.basename(os.path.normpath(path))
        odir = os.path.join(datadir, simname)
        os.mkdir(odir)
        finishedtraj = glob(os.path.join(path, '*.xtc'))
        logger.info("Moving simulation {} to {}.".format(finishedtraj[0], odir))
        move(finishedtraj[0], odir)


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
