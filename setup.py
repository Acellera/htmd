# (c) 2015-2018 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
from setuptools import setup, find_packages
from subprocess import check_output

version = check_output(['git', 'describe', '--tags']).decode('utf8').split('-')[0]

with open('package/htmd-deps/DEPENDENCIES', 'r') as f:
    requirements = f.read().splitlines()
requirements = None  # Some of them are on conda and cannot work with pip install

with open('README.md', 'r') as fh:
    long_description = fh.read()

print('Version [%s]' % (version))
print('Dependencies:')
print(requirements)

setup(name='htmd',
      version=version,
      description='HTMD',
      long_description=long_description,
      zip_safe=False,
      url='https://www.htmd.org',
      maintainer='Acellera Ltd',
      maintainer_email='noreply@acellera.com',
      entry_points={
          'console_scripts': [
              'htmd_register   = htmdx.cli:htmd_do_register',
              'parameterize = htmd.parameterization.cli:main_parameterize',
              'activate_license  = htmdx.cli:main_activate'
          ]
      },

      install_requires=requirements,
      packages=find_packages(include=['htmd*'], exclude=['data']),
      package_data={
        'htmd': ['share/*/*/*/*'],
      }
)
