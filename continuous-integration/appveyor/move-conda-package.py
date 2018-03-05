import os
import glob
import shutil
from conda_build.config import Config


config = Config()
binary_package_glob = os.path.join(config.bldpkgs_dir, 'htmd*.tar.bz2')
binary_package = glob.glob(binary_package_glob)[0]

shutil.move(binary_package, '.')
