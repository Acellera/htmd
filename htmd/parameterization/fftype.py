# (c) 2015-2017 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#

import tempfile
import shutil
import subprocess
import os
from htmd.parameterization.ff import RTF, PRM, AmberRTF, AmberPRM
from enum import Enum


class FFTypeMethod(Enum):
    CHARMM = 1
    AMBER = 2
    CGenFF_2b6 = 1000
    GAFF = 1001
    GAFF2 = 1002


class FFType:
    def __init__(self, mol, method=FFTypeMethod.CGenFF_2b6, rtf=None, prm=None):
        self.frcmod = None
        self.prepi = None
        self.rtf = None
        self.prm = None
        if method == FFTypeMethod.CHARMM:
            self._rtf = RTF(rtf)
            self._prm = PRM(prm)
            # self._makeTopoFromCharmm()
        # elif (method == FFTypeMethod.AMBER):
        #    self._prepi = PREPI(prepi)
        #    self._frcmod = FRCMOD(frcmod)
        #    self._makeTopoFromAmber()
        elif method == FFTypeMethod.GAFF or method == FFTypeMethod.GAFF2:
            antechamber_binary = shutil.which("antechamber")
            if not antechamber_binary:
                raise RuntimeError("antechamber executable not found")
            parmchk2_binary = shutil.which("parmchk2")
            if not parmchk2_binary:
                raise RuntimeError("parmchk2 executable not found")

            cwd = os.getcwd()
            tmpdir = tempfile.mkdtemp()
            try:
                os.chdir(tmpdir)
                mol.write("mol.mol2")
                atomtype = "gaff"
                if method == FFTypeMethod.GAFF2:
                    atomtype = "gaff2"

                subprocess.call(
                    [antechamber_binary, "-at", atomtype, "-nc", str(mol.netcharge), "-fi", "mol2", "-i", "mol.mol2",
                     "-fo", "prepi", "-o", "mol.prepi"])
                subprocess.call([parmchk2_binary, "-f", "prepi", "-i", "mol.prepi", "-o", "mol.frcmod", "-a", "Y"])
                self._rtf = AmberRTF(mol, "mol.prepi", "mol.frcmod")
                self._prm = AmberPRM("mol.prepi", "mol.frcmod")
                os.chdir(cwd)
                shutil.rmtree(tmpdir)
            except:
                os.chdir(cwd)
                raise RuntimeError("FFTyping failed running Antechamber and Parmchk2")
            if not self._rtf or not self._prm:
                raise RuntimeError("FFTyping failed reading Antechamber/Parmchk2 output: see {}".format(tmpdir))

            pass

        elif method == FFTypeMethod.CGenFF_2b6:
            match_binary = shutil.which("match-typer")
            if not match_binary:
                raise RuntimeError("match executable not found")

            cwd = os.getcwd()
            tmpdir = tempfile.mkdtemp()
            try:
                os.chdir(tmpdir)
                mol.write("mol.pdb")
                subprocess.call(
                    [match_binary, "-charge", str(mol.netcharge), "-forcefield", "top_all36_cgenff_new", "mol.pdb"])
                self._rtf = RTF("mol.rtf")
                self._prm = PRM("mol.prm")
                os.chdir(cwd)
                shutil.rmtree(tmpdir)
            except:
                os.chdir(cwd)
                raise RuntimeError("FFTyping failed running Match")
            if not self._rtf or not self._prm:
                raise RuntimeError("FFTyping failed reading Match output: see {}".format(tmpdir))
        else:
            raise RuntimeError("Unknown method for FFType: {}".format(method))
