# (c) 2015-2017 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import os
import re
import unittest
from tempfile import TemporaryDirectory
from subprocess import call


class TestParametrize(unittest.TestCase):

    def test_doc(self):

        reflines = \
            ['\n',
             'HTMD License accepted automatically. Check license here: https://raw.githubusercontent.com/Acellera/htmd/master/htmd/LICENCE.txt\n',
             '\n',
             'For advanced features (e.g. parameterize) and to remove this message, we recommend registering. Run htmd_register in your terminal.\n',
             '\n',
             'Please cite HTMD: Doerr et al.(2016)JCTC,12,1845. https://dx.doi.org/10.1021/acs.jctc.6b00049\n',
             '\n',
             'HTMD Documentation at: https://www.htmd.org/docs/latest/\n',
             '\n',
             'You are on the latest HTMD version (unpackaged : /home/raimis/opt/miniconda3/envs/htmd/lib/python3.6/site-packages/htmd-1.9.4-py3.6.egg/htmd).\n',
             '\n',
             'HTMD: Logging setup failed\n',
             'Deprecation warning: "from htmd import *" will be deprecated. \n',
             'To import all HTMD shortcuts please from now on use "from htmd.ui import *"\n',
             'usage: parameterize [-h] -m <input.mol2> [-l] [-c CHARGE] [--rtf RTF]\n',
             '                    [--prm PRM] [-o OUTDIR] [-t A1-A2-A3-A4] [-n NCPUS]\n',
             '                    [-f {GAFF,GAFF2,CGENFF,all}] [-b {6-31g-star,cc-pVDZ}]\n',
             '                    [--theory {RHF,B3LYP}] [--vacuum] [--no-min] [--no-esp]\n',
             '                    [--no-torsions] [-e {inline,LSF,Slurm}]\n',
             '                    [--qmcode {Gaussian,PSI4,TeraChem}] [--freeze-charge A1]\n',
             '                    [--no-geomopt]\n',
             '\n',
             'Acellera Small Molecule Parameterization Tool\n',
             '\n',
             'optional arguments:\n',
             '  -h, --help            show this help message and exit\n',
             '  -m <input.mol2>, --mol2 <input.mol2>\n',
             '                        Molecule to parameterise, in mol2 format\n',
             '  -l, --list, --list-torsions\n',
             '                        List parameterisable torsions\n',
             '  -c CHARGE, --charge CHARGE\n',
             '                        Net charge on molecule (default: sum of the partial\n',
             '                        charges on the .mol2 file)\n',
             '  --rtf RTF             Inital RTF parameters (req --prm)\n',
             '  --prm PRM             Inital PRM parameters (req --rtf)\n',
             '  -o OUTDIR, --outdir OUTDIR\n',
             '                        Output directory (default: ./)\n',
             '  -t A1-A2-A3-A4, --torsion A1-A2-A3-A4\n',
             '                        Torsion to parameterise (default: all)\n',
             '  -n NCPUS, --ncpus NCPUS\n',
             '                        Number of CPUs to use (default: 8)\n',
             '  -f {GAFF,GAFF2,CGENFF,all}, --forcefield {GAFF,GAFF2,CGENFF,all}\n',
             '                        Inital FF guess to use (default: all)\n',
             '  -b {6-31g-star,cc-pVDZ}, --basis {6-31g-star,cc-pVDZ}\n',
             '                        QM Basis Set (default: cc-pVDZ)\n',
             '  --theory {RHF,B3LYP}  QM Theory (default: B3LYP)\n',
             '  --vacuum              Perform QM calculations in vacuum (default: True)\n',
             '  --no-min              Do not perform QM minimisation (default: False)\n',
             '  --no-esp              Do not perform QM charge fitting (default: False)\n',
             '  --no-torsions         Do not perform torsion fitting (default: False)\n',
             '  -e {inline,LSF,Slurm}, --exec {inline,LSF,Slurm}\n',
             '                        Mode of execution for the QM calculations (default:\n',
             '                        inline)\n',
             '  --qmcode {Gaussian,PSI4,TeraChem}\n',
             '                        QM code (default: PSI4)\n',
             '  --freeze-charge A1    Freeze the charge of the named atom (default: None)\n',
             '  --no-geomopt          Do not perform QM geometry optimisation when fitting\n',
             '                        torsions (default: True)\n']

        with TemporaryDirectory() as tmpDir:
            with open(os.path.join(tmpDir, 'stdout'), mode='w+') as stdout,\
                 open(os.path.join(tmpDir, 'stderr'), mode='w+') as stderr:

                returncode = call(('parameterize', '-h'), stdout=stdout, stderr=stderr, cwd=tmpDir)
                self.assertEqual(returncode, 0)

                stdout.seek(0)
                for ref, res in zip(reflines, stdout.readlines()):
                    if re.search('HTMD version', res):
                        continue
                    self.assertEqual(ref, res)

                stderr.seek(0)
                self.assertEqual(stderr.readlines(), [])

    def test_torsions(self):

        reflines = \
            ['Got charge 0\n',
             'Set charge 0\n',
             'O1     O.1    OG312  O.1-  \n',
             'O2     O.1    OG312  O.1-  \n',
             'H1     H.1    HGA4   H.1() \n',
             'H2     H.1    HGA4   H.1() \n',
             'mol Success!\n',
             '\n',
             'HTMD License accepted automatically. Check license here: https://raw.githubusercontent.com/Acellera/htmd/master/htmd/LICENCE.txt\n',
             '\n',
             'For advanced features (e.g. parameterize) and to remove this message, we recommend registering. Run htmd_register in your terminal.\n',
             '\n',
             'Please cite HTMD: Doerr et al.(2016)JCTC,12,1845. https://dx.doi.org/10.1021/acs.jctc.6b00049\n',
             '\n',
             'HTMD Documentation at: https://www.htmd.org/docs/latest/\n',
             '\n',
             'You are on the latest HTMD version (unpackaged : /home/raimis/opt/miniconda3/envs/htmd/lib/python3.6/site-packages/htmd-1.9.4-py3.6.egg/htmd).\n',
             '\n',
             'HTMD: Logging setup failed\n',
             'Deprecation warning: "from htmd import *" will be deprecated. \n',
             'To import all HTMD shortcuts please from now on use "from htmd.ui import *"\n',
             ' === Listing soft torsions of /home/raimis/prj/htmd.git/htmd/data/test-param/H2O2.mol2 ===\n',
             '\n',
             'Dihedral 0: 2-0-1-3\n',
             'Net Charge: 0\n',
             'Equivalent atom groups:\n',
             ' O1 O2\n',
             ' H1 H2\n',
             'Soft torsions:\n',
             ' H1 O1 O2 H2\n',
             'Detected soft torsions:\n',
             '\tH1-O1-O2-H2\n']

        molFile = os.path.abspath(os.path.join('..', 'data', 'test-param', 'H2O2.mol2'))

        with TemporaryDirectory() as tmpDir:
            with open(os.path.join(tmpDir, 'stdout'), mode='w+') as stdout, \
                    open(os.path.join(tmpDir, 'stderr'), mode='w+') as stderr:

                returncode = call(('parameterize', '-m', molFile, '-l'), stdout=stdout, stderr=stderr, cwd=tmpDir)
                self.assertEqual(returncode, 0)

                stdout.seek(0)
                for ref, res in zip(reflines, stdout.readlines()):
                    if re.search('HTMD version', res):
                        continue
                    self.assertEqual(ref, res)

                stderr.seek(0)
                self.assertEqual(stderr.readlines(), [])


if __name__ == '__main__':
    unittest.main()
