from htmd.apps.app import App
from htmd.apps.acemdlocal import AcemdLocal
from shutil import which, move
from subprocess import call, check_output, CalledProcessError
import threading
import logging
import os
logger = logging.getLogger(__name__)


class PmemdLocal(App):
    """
    Parameters
    ----------
    pmemd_cuda: str
        Path to any executable pmemd binary. If None, it will try to detect
        the pmemd.cuda_SPFP binary
    ngpus : int
        Number of GPU devices that AcemdLocal will use. Each simulation will be
        run on a different GPU. AcemdLocal will use the first `ngpus` devices
        of the machine.
    devices : list
        A list of GPU device indexes on which AcemdLocal is allowed to run
        simulations. Mutually exclusive with `ngpus`
    datadir : str
        A folder to which completed simulations will be moved. If None they
        will be written in the input directory.
    """
    def __init__(self, pmemd_cuda=None, ngpus=None, devices=None, datadir=None):
        self.states = dict()
        self.queue = queue.Queue()
        self.threads = []
        self.shutdown = False

        if ngpus is not None and devices is not None:
            raise ValueError('Parameters `ngpus` and `devices` are mutually exclusive.')

        if ngpus is None and devices is None:
            try:
                devices = range(int(check_output("nvidia-smi -L | wc -l",
                                shell=True).decode("ascii")))
            except:
                raise
        elif ngpus is not None:
            devices = range(ngpus)

        if devices is None:
            raise NameError("Could not determine which GPUs to use. "
                            "Specify the GPUs with the `ngpus=` or `devices=` parameters")
        else:
            logger.info("Using GPU devices {}".format(','.join(map(str, devices))))

        if not pmemd_cuda:
            try:
                # Try to find pmemd.cuda_SPFP in the path
                pmemd_cuda = which("pmemd.cuda_SPFP", mode=os.X_OK)
                if pmemd_cuda:
                    logger.info("Found pmemd.cuda_SPFP at '" + pmemd_cuda + "'")
            except:
                raise NameError("""Cannot find 'pmemd.cuda_SPFP' in the PATH
                                Set its location with the 'pmemd_cuda=' named
                                argument""")

        if not os.access(pmemd_cuda, os.X_OK):
            raise NameError("pmemd_cuda file '" + pmemd_cuda + "' not executable")

        self.threads = []
        for d in devices:
            t = threading.Thread(target=run_job, args=(self, d, pmemd_cuda, datadir))
            t.daemon = True
            t.start()
            self.threads.append(t)

    retrieve = AcemdLocal.retrieve
    submit = AcemdLocal.submit
    inprogress = AcemdLocal.inprogress
