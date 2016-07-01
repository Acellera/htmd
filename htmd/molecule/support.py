# (c) 2015-2016 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
from tempfile import NamedTemporaryFile
import htmd
import numpy
import os
from ctypes import cdll


def xtc_lib():
    lib = {}
    libdir = htmd.home(libDir=True)

    import platform
    if platform.system() == "Windows":
        lib['libc'] = cdll.msvcrt
        cdll.LoadLibrary(os.path.join(libdir, "libgcc_s_seh-1.dll"))
        if os.path.exists(os.path.join(libdir, "psprolib.dll")):
            cdll.LoadLibrary(os.path.join(libdir, "psprolib.dll"))
    else:
        # lib['libc'] = cdll.LoadLibrary("libc.so.6")
        lib['libc'] = cdll.LoadLibrary("libc.{}".format("so.6" if platform.uname()[0] != "Darwin" else "dylib"))

    lib['libxtc'] = cdll.LoadLibrary(os.path.join(libdir, "libxtc.so"))
    return lib


def string_to_tempfile(content, ext):
    f = NamedTemporaryFile(delete=False, suffix="." + ext)
    f.write(content.encode("ascii"))
    f.close()
    return f.name


def pack_string_buffer(data):
    '''data1 = data.astype(dtype=numpy.string_)

#    buf = [create_string_buffer(data[i].encode('ascii')) for i in range(len(data))]
    buf = [create_string_buffer(data1[i]) for i in range(len(data1))]
    ptr = (c_char_p * len(data1))()

    for i in range(0, len(data1)):
        ptr[i] = addressof(buf[i])'''

    # Stefan alternative
    ptr = (c_char_p * len(data))()
    ptr[:] = data.astype(dtype=numpy.string_).tolist()
    #ptr[:] = [x.encode() for x in data]

    return ptr


def pack_int_buffer(data):
    ptr = (c_int * len(data))()
    ptr[:] = data
    '''i = 0
    for d in data:
        ptr[i] = d
        i += 1'''
    return ptr


def pack_ulong_buffer(data):
    ptr = (c_ulong * len(data))()
    ptr[:] = data
    '''i = 0
    for d in data:
        ptr[i] = d
        i += 1'''
    return ptr


def pack_double_buffer(data):
    ptr = (c_double * len(data))()
    ptr[:] = data
    '''i = 0
    for d in data:
        ptr[i] = d
        i += 1'''
    return ptr
