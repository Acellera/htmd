import htmd
import shutil
from glob import glob
import os

# delete previous test
 for directory_i in glob(os.path.join(htmd.home(), 'test_pmemdEngine/*_Test_Protocol_out')):
     shutil.rmtree(directory_i)

# Create Equilibration Protocol
from htmd.protocols.Amber_protocols.AmberProduction import Production
ProdTest = Production()
ProdTest.amber.nstlim = 100000
# for now let us leave it as is and create a folder with all that we need
# to run this test we need Test_Protocol with a prmtop and rst files
# and a (empty) Test_Protocol_out directory, both in htmd.home()
ProdTest.write(os.path.join(htmd.home(), 'Test_Protocol'),
                 os.path.join(htmd.home(), 'Test_Protocol_out'))

# run adaptive sampling with MSM
adapt = htmd.AdaptiveRun()
adapt.app = htmd.apps.pmemdlocal.PmemdLocal(
    pmemd_cuda='/usr/local/amber/bin/pmemd.cuda_SPFP',
    datadir=os.path.join(htmd.home(), 'TrajData'))
adapt.datapath= os.path.join(htmd.home(), 'TrajData')
#adapt.project = 'Test'
adapt.nmin = 2
adapt.nmax = 4
adapt.nepochs = 2
adapt.inputpath = os.path.join(htmd.home(), 'test_pmemdEngine')
adapt.generatorspath = os.path.join(htmd.home(), 'Test_Protocol_out')
adapt.metricsel1 = 'name CA'
#adapt.datapath = 'test_pmemdEngine'
#adapt.ticadim = 0
adapt.updateperiod = 60
adapt.run()