#!/usr/bin/env python3

import htmd
import shutil
from glob import glob
import os
from htmd.protocols.Amber_protocols.AmberProduction import Production

ProdTest = Production()
ProdTest.amber.nstlim = 12500000  #25 ns 
ProdTest.amber.ntpr = 100
ProdTest.amber.ntwx = 100
ProdTest.amber.gamma_ln = 2.0
ProdTest.amber.ntf = 2

ProdTest.amber.parmfile = 'Ala2.prmtop'
ProdTest.amber.coordinates = 'Ala2.rst'
ProdTest.write('./', './ready')


adapt = htmd.AdaptiveRun()
adapt.nmin = 1
adapt.nmax = 4
adapt.nepochs = 10
adapt.metricsel1 = 'name CA'
adapt.metrictype = 'distances'
adapt.ticadim = 5
adapt.updateperiod = 1800 # 30min update 
adapt.filtersel = 'not water'
adapt.app = htmd.apps.pmemdlocal.PmemdLocal(
    pmemd_cuda='/usr/local/amber/bin/pmemd.cuda_SPFP',
    datadir='./data')
adapt.generatorspath = './ready'
adapt.inputpath = './input'
adapt.datapath = './data'
adapt.filteredpath = './filtered'
adapt.run()

