import subprocess
import sys

if sys.argv[1] == 'tag':
    output = subprocess.check_output("git describe --tags", shell=True).decode('utf8')
    tag = output.split('-')[0]
    print(tag)

if sys.argv[1] == 'branch':
    output = subprocess.check_output("git rev-parse --abbrev-ref HEAD", shell=True).decode('utf8')
    if output.startswith('master'):
        print('devel')
    elif output.startswith('rel-'):
        print('stable')
