import platform
import os
import inspect
from htmd.home import home


def path():
    libdir = os.path.join(home(), "lib")
    if os.path.exists(os.path.join(libdir, "basic")):
        libdir = os.path.join(libdir, "basic", platform.system())
    elif os.path.exists(os.path.join(libdir, "pro")):
        libdir = os.path.join(libdir, "pro", platform.system())
    else:
        raise FileNotFoundError('Could not find libs.')
    return libdir