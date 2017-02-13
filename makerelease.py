# (c) 2015-2017 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import os
import shutil
from subprocess import call, check_output


def makeMajorRelease(version, message):
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


def makeMinorDevRelease(version, message):
    call(['git', 'checkout', 'master'])
    call(['git', 'fetch'])
    call(['git', 'pull'])
    call(['git', 'tag', version, '-m', message])
    call(['git', 'push', '--tags', 'origin', 'master'])


def makeStableBugFix(versiontofix, newversion, message):
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
