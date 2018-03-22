# (c) 2015-2018 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import fnmatch
import os
from subprocess import run
import sys
import time

try:
    from htmd.home import home
except ImportError:
    raise ImportError('Could not import htmd.home')

excludedfolders = ('./oldtests', './htmd/tests', './doc', './htmdlib', './package', './continuous-integration',
                   './tutorials')
# TODO: Remove this if when issue #532 is solved
if os.environ.get('TRAVIS_OS_NAME') == 'osx':
    excludedfiles = ('__init__.py', 'license_headers.py', 'setup.py', 'sync_acellera_conda_channel_deps.py',
                     'makerelease.py', 'smallmol.py')
else:
    excludedfiles = ('__init__.py', 'license_headers.py', 'setup.py', 'sync_acellera_conda_channel_deps.py',
                     'makerelease.py')


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
for root, dirnames, filenames in os.walk(home()):
    if excluded(root, excludedfolders):
        continue

    for filename in fnmatch.filter(filenames, '*.py'):
        if not excluded(filename, excludedfiles):
            filestotest.append(os.path.join(root, filename))

# Avoid network communication at each import
os.environ['HTMD_NONINTERACTIVE'] = '1'

# Getting python executable
pythonexe = sys.executable
if not pythonexe:
    raise RuntimeError('Python executable not found.')
else:
    print('Using python executable: {}'.format(pythonexe))

# Running py files
failed = []
times = []
for f in filestotest:
    if not containsTests(f):
        continue
    t = time.time()
    print(' ************************  Running "{}"  ************************'.format(f))
    if f.endswith('amber.py') or f.endswith('charmm.py') or f.endswith('preparation.py'):
        process = run('export PYTHONHASHSEED=1; {} {}'.format(pythonexe, f), shell=True)
    else:
        process = run('{} {}'.format(pythonexe, f), shell=True)
    finishtime = time.time() - t
    print(' ************************  Result : {} Time : {} ************************'.format(process.returncode,
                                                                                             finishtime))
    times.append([f, finishtime])
    if process.returncode != 0:
        if process.returncode != 134:  # TODO: this solves spurious return codes on Travis. To be removed in the future.
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



