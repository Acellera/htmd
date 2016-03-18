import platform
import os
import inspect


def path():
    # TODO: Replace the inspect with something better
    libdir = os.path.join(os.path.dirname(inspect.getfile(path)))
    if os.path.exists(os.path.join(libdir, "basic")):
        libdir = os.path.join(libdir, "basic", platform.system())
    elif os.path.exists(os.path.join(libdir, "pro")):
        libdir = os.path.join(libdir, "pro", platform.system())
    else:
        raise FileNotFoundError('Could not find libs.')
    return libdir