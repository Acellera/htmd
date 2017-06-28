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
from htmd.userinterface import UserInterface
from htmd.decorators import _Deprecated


@_Deprecated('1.5.15', 'htmd.queues.acecloudqueue.AceCloudQueue')
class AWS(UserInterface, App):

    def __init__(self, name=None):
        from acecloud import Cloud, Job
        self._cloud  = Cloud()
        self._group  = name

    def wait( self ):
        from time import sleep
        import sys
        while self.inprogress() != 0:
            sys.stdout.flush()
            sleep(1)

    def submit(self, dirs ):
        from acecloud import Cloud, Job
        if isinstance(mydirs, str): mydirs = [mydirs]

        for d in mydirs:
           name = os.path.basename( os.path.abspath(d))
           Job( ngpus=1, cloud = self._cloud, group=self._group, name=name, path=path )
 
    def retrieve(self):
        jj = self._cloud.getJobs( group=self._group, name="*" )
        for j in jj:
           j.retrieve()

    def inprogress(self):
        cnt=0
        jj = self._cloud.getJobs( group=self._group, name="*" )
        for j in jj:
           s=j.status()
           if( s == Status.RUNNING or s == Status.QUEUED ): cnt=cnt+1
        return cnt

