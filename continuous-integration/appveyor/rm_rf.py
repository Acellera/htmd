# (c) 2015-2018 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
from __future__ import print_function
import os
import sys
import stat
import shutil

def remove_readonly(func, path, excinfo):
    os.chmod(path, stat.S_IWRITE)
    func(path)

def main():
    print(sys.executable)
    try:
        shutil.rmtree(sys.argv[1], onerror=remove_readonly)
    except Exception as e:
        print("Error")
        print(e)

if __name__ == '__main__':
    main()
