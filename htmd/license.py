# (c) 2015-2017 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
from htmd.home import home
from ctypes import *
import os
import platform


def licenseEntitlements():
    libdir = home(libDir=True)

    # No liblicense.so for OS X yet. Need to buy it
    if platform.system() == "Darwin":
        return dict()

    if platform.system() == "Windows":
        cdll.LoadLibrary(os.path.join(libdir, "libgcc_s_seh-1.dll"))
        if os.path.exists(os.path.join(libdir, "psprolib.dll")):
            cdll.LoadLibrary(os.path.join(libdir, "psprolib.dll"))

    lib = cdll.LoadLibrary(os.path.join(libdir, "liblicense.so"))
    arg1 = create_string_buffer(2048)
    arg2 = c_int(2048)
    lib.license_check(arg1, arg2)
    val = (arg1.value.decode("utf-8")).split(":")

    toks = dict()
    for v in val:
        x = v.split("-")
        if len(x) > 1:
            if x[0] not in toks:
                toks[x[0]] = dict()
            if (x[1] not in toks[x[0]]):
                toks[x[0]][x[1]] = 0
            try:
                toks[x[0]][x[1]] = toks[x[0]][x[1]] + int(x[2])
            except:
                pass
    return toks


if __name__ == "__main__":
    print("License entitlements:");
    print(licenseEntitlements())
