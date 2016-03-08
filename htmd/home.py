# (c) 2015-2016 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import htmd
import os
import inspect


def home():
    return os.path.dirname(inspect.getfile(htmd))

#Don't know how to do this
# def modulehome(modname):
#    return os.path.dirname(os.path.dirname(inspect.getfile(modname)))

if __name__ == "__main__":
    h = htmd.home()
    print(h)