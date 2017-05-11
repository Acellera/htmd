# (c) 2015-2017 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import os
import shutil
from subprocess import call, check_output
from tutorials.test_tutorials import test_tutorials


def run_tests(testfolder):
    excode = test_tutorials(testfolder)
    if excode == 0:
        print('All tests were collected and passed successfully.')
    else:
        print('Error occured in running tests.')

    # Error codes taken from pytest documentation
    if excode == 1:
        raise RuntimeError('Tests were collected and run but some of the tests failed.')
    elif excode == 2:
        raise RuntimeError('Test execution was interrupted by the user.')
    elif excode == 3:
        raise RuntimeError('Internal error happened while executing tests.')
    elif excode == 4:
        raise RuntimeError('pytest command line usage error.')
    elif excode == 5:
        raise RuntimeError('No tests were collected.')


def makeMajorRelease(version, message, testfolder):
    run_tests(testfolder)
    relname = 'rel-{}'.format(version)
    call(['git', 'checkout', 'master'])
    call(['git', 'fetch'])
    call(['git', 'pull'])
    call(['git', 'checkout', '-b', relname])
    call(['git', 'tag', version, '-m', message])
    call(['git', 'push', '--tags', 'origin', relname])
    call(['git', 'checkout', 'master'])
    call(['git', 'tag', version, '-m', message])
    call(['git', 'push', '--tags'])


def makeMinorDevRelease(version, message, testfolder):
    run_tests(testfolder)
    call(['git', 'checkout', 'master'])
    call(['git', 'fetch'])
    call(['git', 'pull'])
    call(['git', 'tag', version, '-m', message])
    call(['git', 'push', '--tags', 'origin', 'master'])


def makeStableBugFix(versiontofix, newversion, message, testfolder):
    run_tests(testfolder)
    relname = 'rel-{}'.format(versiontofix)
    call(['git', 'checkout', 'master'])
    call(['git', 'fetch'])
    call(['git', 'pull'])
    call(['git', 'checkout', versiontofix])
    call(['git', 'fetch'])
    call(['git', 'pull'])
    call(['git', 'tag', newversion, '-m', message])
    call(['git', 'push', '--tags', 'origin', relname])
    call(['git', 'checkout', 'master'])
