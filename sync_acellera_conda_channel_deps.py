#!/usr/bin/env python3
# (c) 2015-2017 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#

import os
from subprocess import call

from binstar_client.utils import get_server_api

api = get_server_api()

os.chdir('/tmp')

# Add packages to sync in this list here

packages = [
    ['omnia', 'fftw3f'],
    ['omnia', 'openmm'],
    ['omnia', 'ambermini'],
    ['omnia', 'bhmm'],
    ['omnia', 'funcsigs'],
    ['omnia', 'mdtraj'],
    ['omnia', 'msmtools'],
    ['openbabel', 'openbabel'],
    ['omnia', 'pint'],
    ['omnia', 'progress_reporter'],
    ['omnia', 'pyemma'],
    ['omnia', 'thermotools'],
    ['psi4', 'dftd3'],
    ['bioconda', 'nglview'],
    ['omnia', 'nglview'],  # Adding both omnia and bioconda nglview because they provide different builds
]

for channel, package in packages:
    remote = api.package(channel, package)
    try:
        acellera = api.package('acellera', package)
    except:
        acellera = None

    # For each file on the remove channel of the given version
    for rf in remote['files']:
        if (rf['version'] != remote['latest_version']) or ('py27' in rf['basename']):
            continue

        # Searching for the file in acellera channel
        found = False
        if acellera is not None:
            for af in acellera['files']:
                if (af['version'] == rf['version']) and (af['basename'] == rf['basename']):
                    found = True
                    break

        if not found:
            print('\nSyncing {}/{} file "{}" to acellera'.format(channel, package, rf['basename']))

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

