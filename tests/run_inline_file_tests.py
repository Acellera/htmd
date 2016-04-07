import fnmatch
import os
from subprocess import call, check_output
import sys

excludedfolders = ('./tests',  './doc')
excludedfiles = ('__init__.py')


def excluded(name, exlusionlist):
    for e in exlusionlist:
        if name.startswith(e):
            return True
    return False

# Collecting py files to run
filestotest = []
for root, dirnames, filenames in os.walk('.'):
    if excluded(root, excludedfolders):
        continue

    for filename in fnmatch.filter(filenames, '*.py'):
        if excluded(filename, excludedfiles):
            continue
        filestotest.append(os.path.join(root, filename))

# Running py files
failed = []
for f in filestotest:
    print(' ************************  Running "{}"  ************************'.format(f))
    out = call('python {}'.format(f), shell=True)
    print(' ************************  Result : {}  ************************'.format(out))
    if out != 0:
        failed.append(f)

# Reporting errors
if len(failed) != 0:
    errormsg = 'Tests failed for the following files:\n'
    for f in failed:
        errormsg += '{}\n'.format(f)
    raise AssertionError(errormsg)

print(' All tests passed! ')



