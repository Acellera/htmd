# (c) 2015-2022 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import os


def compareVersions():
    from htmd import __version__ as currver
    from natsort import natsorted
    import os
    import time
    from os.path import expanduser

    __home = expanduser("~")
    __htmdconf = os.path.join(__home, ".htmd")
    if not os.path.exists(__htmdconf):
        try:
            os.makedirs(__htmdconf)
        except Exception:
            print(
                f"Unable to create {__htmdconf} folder. Will not check for new HTMD versions."
            )
            return
    __file = os.path.join(__htmdconf, ".latestversion")

    # Check if one day has passed since last version check. If yes, get new version and write to file
    if (
        not os.path.isfile(__file)
        or time.time() > os.path.getmtime(__file) + 86400
        or os.stat(__file).st_size == 0
    ):
        _writeLatestVersionFile(__file)

    try:
        with open(__file, "r") as f:
            latest = f.read().strip()
    except Exception:
        print(
            f"Unable to open {__file} file for reading. Will not check for new HTMD versions."
        )
        return

    if natsorted((latest, currver))[1] != currver:
        print(
            f"New HTMD version ({latest}) is available. You are currently on ({currver}).\n"
            "We recommend you create a new conda environment with the latest HTMD version.\n"
            f"Run: `conda create -n htmd{latest} htmd={latest} -c acellera -c conda-forge`\n"
        )
    else:
        print(f"You are on the latest HTMD version ({currver}).")

    print("")


def _writeLatestVersionFile(fname):
    import urllib.request
    import json

    version = None
    try:
        url = "https://api.anaconda.org/package/acellera/htmd"
        with urllib.request.urlopen(url) as r:
            version = json.loads(r.read())["latest_version"]
    except Exception:
        return

    try:
        with open(fname, "w") as f:
            f.write(version)
    except Exception:
        print(
            f"Unable to open {fname} file for writing. Will not check for new HTMD versions."
        )
        return

    os.utime(fname, None)
