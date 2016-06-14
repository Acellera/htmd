import htmd
import sys
import shutil

# delete previous test
if os.path.exists('./test_pmemdEngine'):
	shutil.rmtree('./test_pmemdEngine')

adapt = htmd.AdaptiveRun(inputpath='./test_pmemdEngine')
adapt.nmin = 2
adapt.nmax = 3
adapt.nepochs = 2
adapt.ticadim = 0
adapt.metricsel1 = 'name CA'
adapt.generatorspath = './test_pmemdEngine'
adapt.app = htmd.apps.pmemdlocal.PmemdLocal(pmemd_cuda='/usr/local/amber/bin/pmemd.cuda_SPFP')
adapt.run()
