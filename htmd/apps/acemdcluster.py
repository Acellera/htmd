# (c) 2015-2017 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import htmd
from htmd.apps.app import App
from os.path import isdir, abspath, basename, exists,join
from os import getcwd, getlogin, makedirs, chdir
from subprocess import call,check_output
from htmd import UserInterface


class AcemdCluster(UserInterface, App):
    _commands = {'J': None, 'b': None, 'q': 'gpu_priority', 'n': 1, 'W': '23:59', 'R': 'select[ngpus>0] rusage[ngpus_excl_p=1]'}
    _bsubcmd = '''
                #!/bin/csh
                hostname > machine.log
                /usr/bin/nvidia-smi >> machine.log
                acemd input > log
                #TODO: check if success then
                rm restart.*
                '''

    def retrieve(self):
        """ Retrieves completed work units from the engine

        Retrieves completed work units from the engine into a folder provided in the constructor options.

        Example
        -------
        >>> grid.retrieve()
        """
        #Nothing to do as everything is in the same filesystem
        pass

    def submit(self, mydirs):
        """ Submits all work units in a given directory list

        Submits all work units in a given directory list to the engine with the options given in the constructor opt.

        Parameters
        ----------
        mydirs : list of str
            A list or ndarray of directory paths

        Examples
        --------
        >>> grid.submit(glob('input/e2*'));
        """

        if isinstance(mydirs, str): mydirs = [mydirs]

        for d in mydirs:
            if not isdir(d):
                raise NameError('Submit: directory ' + d + ' does not exist.')

        # if all folders exist, submit
        for d in mydirs:
            dirname = abspath(d)
            self.b = dirname
            print('Submitting ' + dirname)
            command = ['bsub']
            command.append(self.cmdlinelist())
            command.append(self._bsubcmd)
            if self.dryrun:
                print(command)
            else:
                call(command)

    def inprogress(self):
        """ Get the number of simulations in progress

        Returns the sum of the number of running and queued workunits of the specific group in the engine.

        Example
        -------
        >>> grid.inprogress()
        """
        cmd = ['bjobs', '-J', self.J, '-noheader', '2']
        output = check_output(cmd)
        # TODO: check lines and handle errors
        return output #wrong

if __name__ == "__main__":
    pass