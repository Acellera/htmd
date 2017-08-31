#!/usr/bin/env python3
# (c) 2015-2017 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#

import os
import sys
from subprocess import call, check_output
import json

from binstar_client.utils import get_server_api

if 'TRAVIS' not in os.environ:
    print('This script should be run on Travis, with a fresh conda installation')
    sys.exit(1)

if os.getenv('TRAVIS_PYTHON_VERSION') == '3.6':
    python_version = 'py36'
elif os.getenv('TRAVIS_PYTHON_VERSION') == '3.5':
    python_version = 'py35'
else:
    print('Python version {} not supported'.format(os.getenv('TRAVIS_PYTHON_VERSION')))
    sys.exit(1)

osname = os.getenv('TRAVIS_OS_NAME')

api = get_server_api()

os.chdir('/tmp')

# Add packages to sync in this list here

# Get all non-defaults installed packages from conda
packages = check_output(['conda', 'list', '--json']).decode('utf8')
packages = json.loads(packages)

nd_packages = list()
for package in packages:
    if package['channel'] != 'defaults':
        nd_packages.append(package)

for package in nd_packages:
    name = package['name']
    channel = package['channel']
    version = package['version']
    build_number = package['build_number']

    remote = api.package(channel, name)
    try:
        acellera = api.package('acellera', name)
    except:
        acellera = None

    # For each file on the remote channel of the given os, version, build, and python version
    for rf in remote['files']:
        rf_version = rf['version']
        rf_build_number = rf['attrs']['build_number']
        rf_build = rf['attrs']['build']
        rf_platform = rf['attrs']['platform']
        rf_basename = rf['basename']

        if rf_version != version or rf_build_number != build_number or python_version not in rf_build or \
                        rf_platform != osname:
            continue

        # Searching for the file in acellera channel
        found = False
        if acellera is not None:
            for af in acellera['files']:
                if af['version'] == rf_version and af['basename'] == rf['basename']:
                    found = True
                    break

        if not found:
            print('\nSyncing {}/{} file "{}" to acellera channel'.format(channel, package, rf['basename']))

            url = 'https:' + rf['download_url']
            print('Downloading {}'.format(url))

            odir = os.path.dirname(os.path.abspath(rf['basename']))
            if not os.path.exists(odir):
                os.makedirs(odir)

            call(['curl', '-L', '-s', url, '-o', rf['basename']])
            print('Uploading...')
            try:
                os.getenv('ANACONDA_TOKEN_BASIC')
                call(['anaconda', 'upload', '--force', '-t', os.getenv('ANACONDA_TOKEN_BASIC'), '-u', 'acellera', rf['basename']])
            except:
                try:
                    call(['anaconda', 'upload', '--force', '-u', 'acellera', rf['basename']])
                except:
                        print('Failed to sync')
            os.remove(rf['basename'])

