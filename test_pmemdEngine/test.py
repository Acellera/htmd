import htmd
import sys
import shutil
from glob import glob

# delete previous test
for directory_i in glob('./test_pmemdEngine/*_test_pmemdEngine'):
	shutil.rmtree(directory_i)

adapt = htmd.AdaptiveRun(inputpath='./test_pmemdEngine')
adapt.nmin = 2
adapt.nmax = 3
adapt.nepochs = 2
adapt.ticadim = 0
adapt.metricsel1 = 'name CA'
adapt.generatorspath = './test_pmemdEngine'
adapt.app = htmd.apps.pmemdlocal.PmemdLocal(pmemd_cuda='/usr/local/amber/bin/pmemd.cuda_SPFP')
adapt.run()
