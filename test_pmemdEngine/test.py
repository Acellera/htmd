import htmd
import sys
import shutil
from glob import glob

# delete previous test
for directory_i in glob('./test_pmemdEngine/*_test_pmemdEngine'):
    shutil.rmtree(directory_i)

# run adaptive sampling with MSM
adapt = htmd.AdaptiveRun()
adapt.app = htmd.apps.pmemdlocal.PmemdLocal(
    pmemd_cuda='/usr/local/amber/bin/pmemd.cuda_SPFP')
adapt.project = 'Test'
adapt.nmin = 4
adapt.nmax = 4
adapt.nepochs = 10
adapt.inputpath = 'test_pmemdEngine'
adapt.generatorspath = 'test_pmemdEngine'
adapt.metricsel1 = 'name CA'
adapt.datapath = 'test_pmemdEngine'
adapt.ticadim = 0
adapt.run()
