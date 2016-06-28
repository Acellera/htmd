import htmd
import shutil
from glob import glob
import os

if os.path.isdir(os.path.join(htmd.home(), 'TrajData')):
    for directory_i in glob(os.path.join(htmd.home(), 'TrajData')):
	     shutil.rmtree(directory_i)
else:
    os.makedirs(os.path.join(htmd.home(), 'TrajData'))

# Create Equilibration Protocol
from htmd.protocols.Amber_protocols.AmberProduction import Production
ProdTest = Production()
ProdTest.amber.nstlim = int(12500000 / 5)
# for now let us leave it as is and create a folder with all that we need
# to run this test we need Test_Protocol with a prmtop and rst files
# and a (empty) Test_Protocol_out directory, both in htmd.home()
ProdTest.write(os.path.join(htmd.home(), 'Test_Protocol'),
                 os.path.join(htmd.home(), 'Test_Protocol_out'))

# delete previous test
for directory_i in glob(adapt.inputpath):
    shutil.rmtree(directory_i)

# run adaptive sampling with MSM
adapt = htmd.AdaptiveRun()
#adapt.datapath= os.path.join(htmd.home(), 'TrajData')
#adapt.project = 'Test'
adapt.nmin = 2
adapt.nmax = 4
adapt.nepochs = 10
#adapt.inputpath = os.path.join(htmd.home(), 'test_pmemdEngine')
adapt.generatorspath = os.path.join(htmd.home(), 'Test_Protocol_out')
adapt.metricsel1 = 'name CA'
adapt.ticadim = 2
adapt.updateperiod = 3600
#adapt.filteredpath = os.path.join(htmd.home(), 'Filtered')
adapt.filtersel='name CA'
adapt.app = htmd.apps.pmemdlocal.PmemdLocal(
    pmemd_cuda='/usr/local/amber/bin/pmemd.cuda_SPFP',
    datadir='data')
adapt.run()