# (c) 2015-2018 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
from setuptools import setup, find_packages
import subprocess

version = subprocess.Popen(["git", "describe", "--tags"], stdout=subprocess.PIPE).stdout.read().decode("utf8")
version = version.split("-")
version = version[0]

f = open("package/htmd-deps/DEPENDENCIES", "r")
deps = None

print("Version [%s]" % (version))
print("Dependencies:")
print(deps)

setup(name='htmd',
      version=version,
      description='HTMD',
      packages=find_packages(exclude=['data']),
      install_requires=deps,
      zip_safe=False,
      url="https://www.htmd.org",
      maintainer="Acellera Ltd",
      maintainer_email="noreply@acellera.com",
      entry_points={
          "console_scripts": [
              'htmdnb    = htmdx.cli:main_htmd_notebook',
              'htmd      = htmdx.cli:main_htmd',
              'htmd_register   = htmdx.cli:main_do_nothing',
              'parameterize = htmd.parameterization.cli:main_parameterize',
              'activate_license  = htmdx.cli:main_activate'
          ]
      },
      package_data={
        'htmd': ['share/*/*/*/*'],
      },
      include_package_data=True
      )
