# (c) 2015-2016 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
from tempfile import NamedTemporaryFile
from ctypes import *
import numpy


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
