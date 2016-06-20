import htmd
import shutil
from glob import glob
import os

# delete previous test
for directory_i in glob(os.path.join(htmd.home(),'test_pmemdEngine/*_Test_Protocol_out')):
    shutil.rmtree(directory_i)

# Create Equilibration Protocol
from htmd.protocols.Amber_protocols.AmberEquilibration import Equilibration
test_equil = Equilibration()
# for now let us leave it as is and create a folder with all that we need
# to run this test we need Test_Protocol with a prmtop and rst files
# and a (empty) Test_Protocol_out directory, both in htmd.home()
test_equil.write(os.path.join(htmd.home(), 'Test_Protocol'),
				 os.path.join(htmd.home(),'Test_Protocol_out'))

# run adaptive sampling with MSM
adapt = htmd.AdaptiveRun()
adapt.app = htmd.apps.pmemdlocal.PmemdLocal(
    pmemd_cuda='/usr/local/amber/bin/pmemd.cuda_SPFP')
#adapt.project = 'Test'
adapt.nmin = 1
adapt.nmax = 1
adapt.nepochs = 2
adapt.inputpath = 'test_pmemdEngine'
adapt.generatorspath = os.path.join(htmd.home(),'Test_Protocol_out')
adapt.metricsel1 = 'name CA'
#adapt.datapath = 'test_pmemdEngine'
#adapt.ticadim = 0
#adapt.updateperiod = 1
adapt.run()
