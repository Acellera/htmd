# (c) 2015-2017 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import fnmatch
import os
from subprocess import call, check_output
import sys
import time

excludedfolders = ('./oldtests', './htmd/tests', './doc', './htmdlib', './package', './continuous-integration', './tutorials')
excludedfiles = ('__init__.py', 'license_headers.py', 'setup.py', 'sync_acellera_conda_channel_deps.py', 'makerelease.py')


def excluded(name, exclusionlist):
    for e in exclusionlist:
        if name.startswith(e):
            return True
    return False


def containsTests(fname):
    with open(fname, 'r') as f:
        lines = f.readlines()
        for l in lines:
            if l.startswith('if __name__'):
                return True
        return False


# Collecting py files to run
filestotest = []
for root, dirnames, filenames in os.walk('.'):
    if excluded(root, excludedfolders):
        continue

    for filename in fnmatch.filter(filenames, '*.py'):
        if not excluded(filename, excludedfiles):
            filestotest.append(os.path.join(root, filename))

# Avoid network communication at each import
os.environ['HTMD_NONINTERACTIVE'] = '1'

# Running py files
failed = []
times = []
for f in filestotest:
    if not containsTests(f):
        continue
    t = time.time()
    print(' ************************  Running "{}"  ************************'.format(f))
    if f.endswith('amber.py') or f.endswith('charmm.py') or f.endswith('preparation.py'):
        out = call('export PYTHONHASHSEED=1; python {}'.format(f), shell=True)
    else:
        out = call('python {}'.format(f), shell=True)
    finishtime = time.time() - t
    print(' ************************  Result : {} Time : {} ************************'.format(out, finishtime))
    times.append([f, finishtime])
    if out != 0:
        failed.append(f)

for p in sorted(times, key=lambda x: x[1]):
    print(p)

# Reporting errors
if len(failed) != 0:
    errormsg = 'Tests failed for the following files:\n'
    for f in failed:
        errormsg += '{}\n'.format(f)
    raise AssertionError(errormsg)

print(' All tests passed! ')



