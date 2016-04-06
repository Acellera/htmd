# -*- coding: utf-8 -*-
"""
Created on Tue Mar  8 14:42:58 2016

@author: sdoerr
"""

from htmd import *
from systempreparation import systemPreparation
import shutil

mol = Molecule('3ptb')
mol.filter('protein or water')
#mol.view()
molp = systemPreparation(mol, ffout='AMBER')
#molp.view()

molp.set('segid','P','protein')
molp.set('segid','W','water')

shutil.rmtree("build-charmm",ignore_errors=True)
mol_charmm = charmm.build(molp, outdir='./build-charmm/', ionize=False)

shutil.rmtree("build-amber",ignore_errors=True)
mol_amber  = amber.build(molp, outdir='./build-amber/', ionize=False)

