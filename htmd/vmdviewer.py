#!/usr/bin/env python
# (c) 2015-2017 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import subprocess
import inspect
import threading
import platform

try:
    # Python 2
    import Queue as queue
except:
    # Python 3
    import queue

    pass

import time
import os
import sys
from htmd.molecule.support import *
import numpy as np
import tempfile

_viewers = np.empty(0, dtype=object)


def _enqueue_output(obj, vmd, queue):
    while vmd.poll() is None:
        line = vmd.stdout.readline()
    # print("READIN: [" + line.decode("ascii") + "]", end="" );
    #			queue.put(line.decode("ascii"))
    #	print( "Process finished" )
    vmd.stdout.close()
    obj.done = 1


class VMD:
    """ Please do not directly call this class constructor. Use the `viewer` or `getCurrentViewer` function instead.
    """

    def __init__(self, vmd=None, host=None, dispdev='win'):
        self.done = 0
        vmd = getVMDpath(vmd=vmd)

        args = [vmd]
        if host:
            args.append('--host')
            args.append(host)
        args.append('--dispdev')
        args.append(dispdev)
        self.vmd = subprocess.Popen(args, stdin=subprocess.PIPE, stdout=subprocess.PIPE, bufsize=0, close_fds=True,
                                    shell=False)
        self.queue = queue.Queue()
        self.thread = threading.Thread(target=_enqueue_output, args=(self, self.vmd, self.queue))
        self.thread.daemon = True
        self.thread.start()
        #	self.vmd.stdin.write( b"chan configure stdout -buffering none\n" )
        time.sleep(2)
        self.send('display reposition 500 1000')
        self.send('display resize 800 800')
        self.send('menu main on')
        self.send('display depthcue off')
        self.send('axes location Off')
        # self.send('color Display Background white')
        self.send('display projection Orthographic')

    def send(self, command):
        """ Send a tcl command to VMD

        Parameters
        ----------
        command : str
            The tcl command which to send
        """
        if self.done:
            return
        # print( command )
        fn = string_to_tempfile(command, "tcl")
        #		print(fn)

        cc = "if { [ catch { source " + fn + "} ] } { unlink " + fn + " } else { unlink " + fn + " }\n"
        self.vmd.stdin.write(cc.encode("ascii"))
        self.vmd.stdin.flush()

        while os.path.isfile(fn) and not self.done:
            time.sleep(0.01)

    def loadMol(self, mol, name=None):
        """ Load a :class:`Molecule <htmd.molecule.molecule.Molecule>` object into VMD

        Parameters
        ----------
        mol : :class:`Molecule <htmd.molecule.molecule.Molecule>` object
            The Molecule to load into VMD
        name : str
            The name to give to the Molecule
        """
        # Getting name of variable passed to method
        '''
        caller = inspect.currentframe().f_back
        try:
            value = caller.f_locals[mol]
        except KeyError:
            value = caller.f_globals[mol]
        if name is None:
            name = value
        # -----------------------------------------
        '''
        mol.view(name=name)

    def rep(self, mode, sel='resname MOL', color=0):
        """ Modify representations for the top molecule in VMD

        Parameters
        ----------
        mode : ['protein','ligand']
            If set to 'protein' it will show a single conformation of the protein with sequence coloring
            If set to 'ligand' it will show all ligands in the specified `color` and the protein in teal
        sel : str
            Atomselection for the ligand
        color : int
            Color for the ligand. Use color numbers of VMD
        """
        if mode == 'ligand':
            # Protein representation
            self.send('mol modstyle 0 top NewCartoon')
            self.send('mol modselect 0 top "protein"')
            # self.send('mol modcolor 0 top Index')
            # Ligand representation
            self.send('set rep [molinfo top get numreps]')
            self.send('mol color ColorID ' + str(color))
            self.send('mol representation Lines 1.000000')
            self.send('mol selection ' + sel)
            self.send('mol material Opaque')
            self.send('mol addrep top')
            self.send('set end [molinfo top get numframes]')
            self.send('mol drawframes top $rep "0:$end"')
        elif mode == 'protein':
            # Protein representation
            self.send('mol modstyle 0 top NewCartoon')
            self.send('mol modselect 0 top "protein"')
            self.send('mol modcolor 0 top Index')
        else:
            raise NameError('Invalid mode. Choose between ''ligand'' and ''protein''')

    def completed(self):
        """ Check if the viewer has been closed

        Returns
        -------
        comp : bool
            Returns True if the VMD viewer has closed
        """
        return self.done

    def copy(self):
        return None

    def render(self, outfile, renderer='tachyon', resolution=(1920, 1080), aasamples=36, skylight=1, tachyon=None, convert=None, trim=False):
        """ Renders the current VMD scene into a file.

        Parameters
        ----------
        outfile : str
            File to which to render image
        renderer : ('tachyon', 'snapshot')
            Which renderer to use
        resolution : tuple
            X,Y resolution of the output image
        aasamples : int
            Number of anti-aliasing samples
        skylight : float
            Add a skylight
        tachyon : str
            Path to tachyon renderer executable
        convert : bool
            Attempts to convert the image to the datatype of the `outfile` extension
        trim : bool
            Trims the whitespace of the image
        """
        import shutil
        outfile = os.path.abspath(outfile)
        outname, ext = os.path.splitext(outfile)

        if renderer == 'tachyon':
            tmpext = '.psd'
            if tachyon is None:
                try:
                    tachyon = shutil.which('tachyon', mode=os.X_OK)
                except:
                    raise FileNotFoundError("Could not find `tachyon` executable, or no execute permissions are given. Try using renderer='snapshot' instead.")
            rendercommand = 'render Tachyon {out} "{tachyon}" -res {resx} {resy} -aasamples {aa} -add_skylight {sl} {out} -format PSD48 -o {out}.psd'.format(tachyon=tachyon, resx=resolution[0], resy=resolution[1], aa=aasamples, sl=skylight, out=outname)
        elif renderer == 'snapshot':
            tmpext = '.tga'
            rendercommand = 'render snapshot {}{}'.format(outname, tmpext)

        self.send(rendercommand)
        print(rendercommand)
        if not os.path.exists(outname+tmpext):
            raise RuntimeError('Rendering failed to produce image with following command: {}'.format(rendercommand))
        if os.path.exists(outname):
            os.remove(outname)

        if ext != tmpext:
            if convert is None:
                try:
                    convert = shutil.which('convert', mode=os.X_OK)
                except:
                    raise FileNotFoundError('Could not find `convert` executable, or no execute permissions are given. You can find the temporary render file in {}', outname+tmpext)
            os.system('{convert} {outname}{tmpext} {outname}{ext}'.format(convert=convert, tmpext=tmpext, outname=outname, ext=ext))
            if trim:
                os.system('{convert} {outname}{ext} -trim {outname}{ext}'.format(convert=convert, outname=outname, ext=ext))
            os.remove(outname+tmpext)


def getVMDpath(vmd=None):
    sys = platform.system()
    if not vmd:
        if sys == "Linux" or sys == "Darwin":
            vmd = os.path.join(os.path.dirname(inspect.getfile(VMD)), "vmd_wrapper")
        elif sys == "Windows":
            vmd = os.path.join(os.path.dirname(inspect.getfile(VMD)), "vmd_wrapper.bat")
        else:
            raise OSError("Don't know how to run VMD on platform [" + sys + "]")
    if (not vmd) or (not os.access(vmd, os.X_OK)):
        raise OSError("Cannot find VMD. Specify the location with 'vmd=' argument")
    return vmd


def getCurrentViewer(dispdev='win'):
    """ Get the handle to the current molecular viewer

    Returns
    -------
    viewer : :class:`VMD` object
        Returns a handle to the current viewer
    """
    global _viewers
    # Check for dead viewers and remove them
    todrop = np.empty(0, dtype=int)
    for i in range(len(_viewers)):
        if _viewers[i].completed():
            todrop = np.append(todrop, i)
    _viewers = np.delete(_viewers, todrop)
    # Creating a new viewer if none exist or otherwise returning the last viewer
    if len(_viewers) == 0:
        _viewers = np.append(_viewers, VMD(dispdev=dispdev))
    return _viewers[-1]


def viewer(dispdev='win'):
    """ Start a new molecular viewer

    Returns
    -------
    viewer : :class:`VMD` object
        Returns a handle to a new viewer
    """
    global _viewers
    # Check for dead viewers and remove them
    todrop = np.empty(0, dtype=int)
    for i in range(len(_viewers)):
        if _viewers[i].completed:
            todrop = np.append(todrop, i)
    _viewers = np.delete(_viewers, todrop)
    # Creating a new viewer
    _viewers = np.append(_viewers, VMD(dispdev=dispdev))
    print(_viewers)


def _tempfilename():
    return os.path.join(tempfile._get_default_tempdir(), next(tempfile._get_candidate_names()))


if __name__ == "__main__":
    vmd = VMD()
    vmd.send('menu main off')
    vmd.send('menu main on')
    vmd.send('exit')
