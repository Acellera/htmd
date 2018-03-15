import json
from subprocess import check_output
import sys
import os

workdir = sys.argv[1]

text = """
package:
  name: htmd-deps
  version: {{ environ.get('BUILD_VERSION') }}

source:
   path: .

requirements:
  build:
    - python
    - requests

  run:
"""

# Read in all dependencies
deps = []
with open(os.path.join(workdir, 'DEPENDENCIES'), 'r') as f:
    for line in f:
        # Remove comments and split
        packagematch = line.split('#')[0].split()
        # Append dict of the conda package match specification
        if packagematch:
            deps.append(dict(zip(['name', 'version', 'build_string'], packagematch)))

# Get all installed packages from conda
packages = check_output(['conda', 'list', '--json']).decode('utf8')
packages = json.loads(packages)

# Find the version of each dependency and add them to the meta.yaml file
found = {dep['name']: False for dep in deps}
for p in packages:
    name = p['name'].lower()
    if name in [dep['name'] for dep in deps]:
        dep = next(dep for dep in deps if dep['name'] == name)
        if 'version' not in dep:
            text += '    - {} >={}\n'.format(name, p['version'])
        else:
            text += '    - {} {}\n'.format(name, dep['version'])
        found[name] = True

# Check if all dependencies were found. If not print which and exit with error
if not all(found):
    for key in found:
        if not found[key]:
            print('Could not find dependency "{}". Exiting with error.'.format(key))
    sys.exit(1)

# Write htmd-deps meta.yaml file
with open(os.path.join(workdir, 'meta.yaml'), 'w') as f:
    f.write(text)
