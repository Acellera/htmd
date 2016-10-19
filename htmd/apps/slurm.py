# (c) 2015-2016 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
# (c) 2015 Acellera Ltd
# All Rights Reserved
# Distributed under HTMD Software Academic License Agreement v1.0
# No redistribution in whole or part
#
import os
import pwd
import shutil
from os.path import isdir
from subprocess import check_output

from htmd import UserInterface


class Slurm(UserInterface):
    _commands = {
        'name': None,  # whatever identifier you want for the job
        'partition': "gpu",  # the 'queue' to run on
        'priority': 'gpu_priority',
        'resources': 'gpu:1',
        'memory': 4000,  # MB
        'walltime': '23:59:00',  # hh:mm:ss
        'executable': 'acemd',  # The thing to run
        'environment': "ACEMD_HOME,HTMD_LICENSE_FILE"  # Envvars to propagate to the job
    }

    def __init__(self, **kw):
        for i in self._commands.keys():
            self.__dict__[i] = self._commands[i]

        for key, val in kw.items():
            if key in self._commands:
                self.__dict__[key] = val
            else:
                raise ValueError("Invalid configuration option [{}]".format(key))

        if not self.name:
            raise ValueError("Name must be set")
        # Find executables
        self._sbatch = self._find_binary("sbatch")
        self._squeue = self._find_binary("squeue")

        self._exe = self._find_binary(self.executable)

    def _find_binary(self, binary):
        ret = shutil.which(binary, mode=os.X_OK)
        if not ret:
            raise ValueError("Could not find required executable [{}]".format(binary))
        ret = os.path.abspath(ret)
        return ret

    def retrieve(self):
        # Nothing to do as everything is in the same filesystem
        pass

    def submit(self, mydirs, debug=False):
        """SUBMIT - Submits all work units in a given directory list to
                     the engine with the options given in the constructor
                     opt.
             INPUT
                   - dirs: Directories which to submit.
        """

        if isinstance(mydirs, str):
            mydirs = [mydirs]

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

            # sbatch -J $jobname -p $partition /share/PI/rondror/software/Slurm/acemd.sbatch $basedir $jobname
            # $partition $TEMPBASE
            command = [self._sbatch,
                       "-J", str(self.name),
                       "-p", str(self.partition),
                       "--qos", str(self.partition),
                       "--gres", str(self.resources),
                       "--time", str(self.walltime),
                       "--mem", str(self.memory),
                       "--priority", str(self.priority),
                       "--export", str(self.environment),
                       "-D", str(dirname),
                       js
                       ]
            if debug:
                print(command)
            try:
                ret = check_output(command)
                if debug:
                    print(ret)
            except:
                raise
                # raise ValueError( "Submission of [%s] failed" % (dirname) )

    def _make_jobscript(self, path, exe):
        fn = os.path.join(path, "run.sh")
        f = open(fn, "w")
        print("#!/bin/sh", file=f)
        print("{} > log.txt 2>&1".format(exe), file=f)
        f.close()
        os.chmod(fn, 0o700)
        return os.path.abspath(fn)

    def inprogress(self, debug=False):
        """ INPROGRESS - Returns the sum of the number of running and queued
                      workunits of the specific group in the engine.

             SYNOPSIS
                   grid.inprogress()
             EXAMPLE
                   grid.inprogress()
        """

        user = pwd.getpwuid(os.getuid()).pw_name
        cmd = [self._squeue, "-n", self.name, "-u", user, "-p", self.partition]
        if debug:
            print(cmd)

        ret = check_output(cmd)
        if debug:
            print(ret.decode("ascii"))
        # TODO: check lines and handle errors
        l = ret.decode("ascii").split("\n")
        l = len(l) - 2
        if l < 0:
            l = 0  # something odd happened
        return l


if __name__ == "__main__":
    """
    s=Slurm( name="testy", partition="gpu")
    s.submit("test/dhfr1" )
    ret= s.inprogress( debug=False)
    print(ret)
    print(s)
    pass
    """
