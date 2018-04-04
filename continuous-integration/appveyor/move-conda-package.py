# (c) 2015-2018 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import os
import glob
import shutil
from conda_build.config import Config


config = Config()
binary_package_glob = os.path.join(config.bldpkgs_dir, 'htmd*.tar.bz2')
binary_package = glob.glob(binary_package_glob)[0]

shutil.move(binary_package, '.')
