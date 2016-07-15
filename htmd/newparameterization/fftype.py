# (c) 2015-2016 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#

from htmd.molecule.molecule import Molecule
import tempfile
import shutil
import subprocess
import os
from htmd.newparameterization.ff import RTF,PRM
from enum import Enum

class FFTypeMethod(Enum):
  CHARMM     = 1
  AMBER      = 2
  CGenFF_2b6 = 1000
  GAFF       = 1001
  GAFF2      = 1002

class FFType:
  def __init__(self, mol, method=FFTypeMethod.CGenFF_2b6, rtf=None, prm=None, frcmod=None, prepi=None ):
    self.frcmod = None
    self.prepi  = None
    self.rtf = None
    self.prm = None
    if(method == FFTypeMethod.CHARMM ):
       self._rtf = RTF(rtf)
       self._prm = PRM(prm)
       self._makeTopoFromCharmm()
    elif(method == FFTypeMethod.AMBER ):
       self._prepi = PREPI(prepi)
       self._frcmod= FRCMOD( frcmod )
       self._makeTopoFromAmber()
    elif(method == FFTypeMethod.GAFF or method == FFTypeMethod.GAFF2):
      antechamber_binary = shutil.which("antechamber")
      if( not amberchamber_binary ): raise RuntimeError("antechamber executable not found")
      self._makeTopoFromAmber()
      pass
    elif(method == FFTypeMethod.CGenFF_2b6 ):
      match_binary = shutil.which("match")
      if( not match_binary ): raise RuntimeError("match executable not found")

      cwd = os.getcwd()
      tmpdir = tempfile.mkdtemp()
      self._makeTopoFromCharmm()
      try:
        os.chdir( tmpdir )
        mol.write( "mol.pdb" )
        subprocess.call( [match_binary, "-charge", str(mol.netcharge), "-forcefield", "top_all36_cgenff_new", "mol.pdb" ] ) 
        self._rtf = RTF( "mol.rtf" )
        self._prm = PRM( "mol.prm" )
        os.chdir(cwd)
        shutil.rmtree(tmpdir)     
      except:
        os.chdir(cwd)
        raise RuntimeError("FFTyping failed running Match")
      if not self._rtf or not self._prm: 
        raise RuntimeError("FFTyping failed reading Match output: see %s" % (tmpdir))

    else:
      raise RuntimeError("Unknown method for FFType: %s" %(method) )

  def _makeTopoFromCharmm( self ):
    pass     

  def _makeTopoFromAmber( self ):
    pass     
