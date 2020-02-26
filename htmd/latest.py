# (c) 2015-2018 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#


def compareVersions():
    from htmd.version import version
    from natsort import natsorted
    import os
    import time
    from os.path import expanduser

    __home = expanduser("~")
    __htmdconf = os.path.join(__home, ".htmd")
    if not os.path.exists(__htmdconf):
        try:
            os.makedirs(__htmdconf)
        except:
            print(
                "Unable to create {} folder. Will not check for new HTMD versions.".format(
                    __htmdconf
                )
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
        f = open(__file, "r")
    except:
        print(
            "Unable to open {} file for reading. Will not check for new HTMD versions.".format(
                __file
            )
        )
        return
    latestversions = f.readlines()
    f.close()

    if len(latestversions) != 2:
        print(
            "There is something wrong with your {} file. Will not check for new HTMD versions.".format(
                __file
            )
        )
        return

    currver = version()
    if currver != "unpackaged":
        if _is_stable(currver):
            pieces = latestversions[0].split()
            latest = pieces[0].strip()
            verstring = "stable"
        else:
            pieces = latestversions[1].split()
            latest = pieces[0].strip()
            verstring = "devel"
        pydeps = ""
        if len(pieces) > 1:
            pydeps = " python[{}]".format(pieces[1])
    elif currver == "unpackaged":
        pass
    else:
        return

    if currver != "unpackaged" and natsorted((latest, currver))[1] != currver:
        print(
            "New {verstring} HTMD version ({latest}{pydeps}) is available. You are currently on ({currver})."
            "There are several methods to update:"
            "    - Create a new conda env. using `conda create -n htmd{latest} htmd={latest} -c acellera -c psi4 -c conda-forge`"
            "    - Create a brand new conda installation and run `conda install htmd -c acellera -c psi4 -c conda-forge`"
            "    - Run: `conda update htmd -c acellera -c psi4 -c conda-forge` (NOT RECOMMENDED)"
            "".format(
                verstring=verstring, latest=latest, pydeps=pydeps, currver=currver
            )
        )
    else:
        if currver != "unpackaged":
            print("You are on the latest HTMD version ({}).".format(currver))
        else:
            from htmd.home import home

            print(
                "You are on the latest HTMD version ({} : {}).".format(currver, home())
            )

    print("")


def _writeLatestVersionFile(fname):
    import os

    try:
        f = open(fname, "w")
    except:
        print(
            "Unable to open {} file for writing. Will not check for new HTMD versions.".format(
                fname
            )
        )
        return

    try:
        from binstar_client.utils import get_server_api

        api = get_server_api(log_level=0)
        package = api.package("acellera", "htmd")
    except Exception as err:
        print(
            "Failed at checking latest conda version. ({})".format(type(err).__name__)
        )
        return

    stable, dev, stabledeps, devdeps = _release_version(package)

    f.write("{} {}\n{} {}".format(stable, ",".join(stabledeps), dev, ",".join(devdeps)))
    os.utime(fname, None)
    f.close()


def _release_version(package):

    versionlist = package["versions"]
    laststable = None
    lastdev = None
    for ver in versionlist[::-1]:  # Iterate in reverse due to sorting of conda versions
        if laststable is None and _is_stable(ver):
            laststable = ver
        elif lastdev is None:
            lastdev = ver
        if laststable and lastdev:
            break

    stabledeps = _release_python_dep(package, laststable)
    devdeps = _release_python_dep(package, lastdev)

    return laststable, lastdev, stabledeps, devdeps


def _is_stable(ver):
    import numpy as np

    if (
        np.mod(int(ver.split(".")[1]), 2) == 0
    ):  # Even versions are stable (i.e. 1.2.x 3.4.x etc)
        return True
    else:
        return False


def _release_python_dep(package, version, opersys=None):
    import platform

    if opersys is None:
        opersys = platform.system().lower()
        if opersys == "windows":
            opersys = "win"
    try:
        versions = []
        for f in package["files"]:
            if f["version"] == version and f["attrs"][
                "operatingsystem"
            ].lower().startswith(opersys.lower()):
                for d in f["dependencies"]["depends"]:
                    if d["name"].lower() == "python":
                        versions.append(d["specs"][0][1])
        if len(versions):
            return versions
        else:
            return " does not exist for your platform. Please create an issue on HTMD git issue tracker."
    except:
        return

