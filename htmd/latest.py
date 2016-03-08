# (c) 2015-2016 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
def compareVersions():
    from htmd.home import home
    from htmd.version import version
    from natsort import natsorted
    import os
    import time
    #t = time.time()
    __file = os.path.join(home(), '.latestversion')

    if not os.path.isfile(__file) or time.time() > os.path.getmtime(__file) + 86400: #86400:  # Check if one day has passed since last version check
        _writeLatestVersionFile(__file)

    f = open(__file, 'r')
    latestver = f.read()
    f.close()
    currver = version()
    if currver != 'unpackaged' and natsorted((latestver, currver))[1] != currver:
        print('New HTMD version ({}) is available. You are currently on ({}). Use \'conda update htmd\' to update to the new version.'.format(latestver, currver))
    else:
        print('You are on the latest HTMD version ({}).'.format(currver))
    #elapse = time.time() - t
    #print(elapse)


def _writeLatestVersionFile(fname):
    import os
    try:
        ver = _release_version('acellera', 'htmd')
    except Exception as err:
        print("{}".format(err))
        return

    f = open(fname, 'w')
    f.write(ver)
    os.utime(fname, None)
    f.close()


def _release_version(user, package):
    import requests
    from binstar_client.utils import get_binstar
    try:
        requests.get('http://conda.anaconda.org', timeout=1.)
        # Do conda version check
    except requests.Timeout:
        raise NameError('Failed connecting to http://conda.anaconda.org. Cannot check for new HTMD versions.')

    binstar = get_binstar()
    package = binstar.package(user, package)
    return package['versions'][-1]

