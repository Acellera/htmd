# (c) 2015-2016 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#


def compareVersions():
    from htmd.version import version
    from natsort import natsorted
    import os
    import time
    #t = time.time()
    from os.path import expanduser
    __home = expanduser("~")
    __htmdconf = os.path.join(__home, '.htmd')
    if not os.path.exists(__htmdconf):
        try:
            os.makedirs(__htmdconf)
        except:
            print('Unable to create {} folder. Will not check for new HTMD versions.'.format(__htmdconf))
            return
    __file = os.path.join(__htmdconf, '.latestversion')

    # Check if one day has passed since last version check. If yes, get new version and write to file
    if not os.path.isfile(__file) or time.time() > os.path.getmtime(__file) + 86400:
        _writeLatestVersionFile(__file)

    try:
        f = open(__file, 'r')
    except:
        print('Unable to open {} file for reading. Will not check for new HTMD versions.'.format(__file))
        return
    latestversions = f.readlines()
    f.close()

    currver = version()
    if currver != 'unpackaged' and len(latestversions) == 2:  # Supporting users which still haven't passed the one day limit. Can remove in next
        if _is_stable(currver):
            latest = latestversions[0].strip()
            verstring = 'stable'
        else:
            latest = latestversions[1].strip()
            verstring = 'devel'
    elif len(latestversions) == 1:  # Supporting users which still haven't passed the one day limit. Can remove in next version
        latest = latestversions[0].strip()
        verstring = ''
    else:
        return
    if currver != 'unpackaged' and natsorted((latest, currver))[1] != currver:
        print('New {} HTMD version ({}) is available. You are currently on ({}). Use \'conda update htmd\' to '
              'update to the new version.'.format(verstring, latest, currver))
    else:
        print('You are on the latest HTMD version ({}).'.format(currver))


def _writeLatestVersionFile(fname):
    import os
    try:
        f = open(fname, 'w')
    except:
        print('Unable to open {} file for writing. Will not check for new HTMD versions.'.format(fname))
        return
    
    try:
        stable, dev = _release_version('acellera', 'htmd')
    except Exception as err:
        print("Failed at checking latest conda version. {}".format(err))
        f.close()
        return
    
    f.write('{}\n{}'.format(stable, dev))
    os.utime(fname, None)
    f.close()


def _release_version(user, package):
    import requests
    from binstar_client.utils import get_server_api

    api = get_server_api()
    package = api.package(user, package)

    versionlist = package['versions']
    laststable = None
    lastdev = None
    for ver in versionlist[::-1]:  # Iterate in reverse due to sorting of conda versions
        if _is_stable(ver):
            laststable = ver
        else:
            lastdev = ver
        if laststable and lastdev:
            break

    return laststable, lastdev


def _is_stable(ver):
    import numpy as np
    if np.mod(int(ver.split('.')[1]), 2) == 0:  # Even versions are stable (i.e. 1.2.x 3.4.x etc)
        return True
    else:
        return False

