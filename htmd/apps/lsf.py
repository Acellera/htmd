# (c) 2015-2017 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import shutil
from os.path import isdir
from subprocess import check_output
import os
import string
import random
from htmd.userinterface import UserInterface
from htmd.decorators import _Deprecated


@_Deprecated('1.5.15', 'htmd.queues.lsfqueue.LsfQueue')
class LSF(UserInterface):
    _commands = {
        'name': "",
        'queue': "gpu_priority",
        'ncpus': 1,
        'resources': 'select[ngpus>0] rusage[ngpus_excl_p=1]',
        'memory': "4000",  # MB
        'walltime': '23:59',  # hh:mm:ss
        'executable': 'acemd',  # The thing to run
        'app': "",
    }

    def __init__(self, **kw):
        for i in self._commands.keys():
            self.__dict__[i] = self._commands[i]

        for key, val in kw.items():
            if key in self._commands:
                self.__dict__[key] = val
            else:
                raise ValueError("Invalid configuration option [{}]".format(key))

        # Find executables
        self._bsub = self._find_binary("bsub")
        self._bjobs = self._find_binary("bjobs")

        self._dirs = []
        try:
            self._exe = self._find_binary(self.executable)
        except:
            self._exe = self.executable

        lst = [random.choice(string.ascii_letters + string.digits) for _ in range(10)]
        self.__dict__["name"] = "".join(lst)

    def _find_binary(self, binary):
        ret = shutil.which(binary, mode=os.X_OK)
        if not ret:
            raise ValueError("Could not find required executable [{}]".format(binary))
        ret = os.path.abspath(ret)
        return ret

    def retrieve(self):
        # Nothing to do as everything is in the same filesystem
        pass

    def wait(self):
        from time import sleep
        import sys
        while self.inprogress() != 0:
            sys.stdout.flush()
            sleep(1)

    def submit(self, mydirs, debug=False):
        """ SUBMIT - Submits all work units in a given directory list to
                     the engine with the options given in the constructor
                     opt.
             INPUT
                   - dirs: Directories which to submit.
        """

        if isinstance(mydirs, str):
            mydirs = [mydirs]
        self._dirs.extend(mydirs)
        for d in mydirs:
            if not isdir(d):
                raise NameError('Submit: directory ' + d + ' does not exist.')

        # if all folders exist, submit
        for d in mydirs:
            dirname = os.path.abspath(d)
            #            self.b = dirname
            if debug:
                print('Submitting ' + dirname)

            js = self._make_jobscript(dirname, self._exe)

            if not self.app:
                cmd = [
                    self._bsub,
                    "-J", self.name,
                    "-M", self.memory,
                    "-W", self.walltime,
                    "-q", self.queue,
                    "-R", self.resources,
                    "-n", str(self.ncpus),
                    js
                ]
            else:
                cmd = [
                    self._bsub,
                    "-J", self.name,
                    "-M", self.memory,
                    "-W", self.walltime,
                    "-q", self.queue,
                    "-R", self.resources,
                    "-n", str(self.ncpus),
                    "-app", self.app,
                    js
                ]
            if debug:
                print(cmd)
            try:
                ret = check_output(cmd)
                if debug:
                    print(ret)
            except:
                raise
                # raise ValueError( "Submission of [%s] failed" % (dirname) )

    def _make_jobscript(self, path, exe):
        fn = os.path.join(path, "run.sh")
        f = open(fn, "w")
        print("#!/bin/sh", file=f)
        print("cd \"" + path + "\"", file=f)
        print("module load acemd", file=f)
        print("module load acellera/test", file=f)
        print("module load gaussian", file=f)
        if "acemd" in exe:
            print("{} --device $CUDA_VISIBLE_DEVICES > log.txt 2>&1".format(exe), file=f)
        else:
            print("{}".format(exe), file=f)
        print("touch .done", file=f)

        f.close()
        os.chmod(fn, 0o700)
        return os.path.abspath(fn)

    def inprogress(self, debug=False):
        inprogress = 0
        for i in self._dirs:
            if not os.path.exists(os.path.join(i, ".done")):
                inprogress += 1
        return inprogress

    def __inprogress(self, debug=False):
        """ INPROGRESS - Returns the sum of the number of running and queued
                      workunits of the specific group in the engine.

             SYNOPSIS
                   grid.inprogress()
             EXAMPLE
                   grid.inprogress()
        """
        cmd = [self._bjobs, "-J", self.name, "-sum", "-noheader"]
        if debug:
            print(cmd)

        ret = check_output(cmd)
        if debug:
            print(ret.decode("ascii"))
        # TODO: check lines and handle errors
        l = ret.decode("ascii").split()
        return int(l[0]) + int(l[4])

    # if __name__ == "__main__":
    #    s=LSF( name="adaptivetest" )
    #    print(s.name)
    """
#    s.submit( "test/dhfr", debug=True )
#    s.submit( "test/dhfr", debug=True )
    ret= s.inprogress( debug=False )
    print(ret)
    pass
    """
