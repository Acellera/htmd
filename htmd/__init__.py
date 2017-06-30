from htmd.ui import *  # TODO: This will be deprecated in a few versions
from htmd.config import config
from htmd.version import version as _version
import os

config()
__version__ = _version()

if os.getenv('HTMD_CONFIG'):
    configfile = os.getenv('HTMD_CONFIG')
    if not os.path.isfile(configfile):
        raise FileNotFoundError('HTMD Config file {} does not exist (Check HTMD_CONFIG).'.format(configfile))
    elif not configfile.endswith('.py'):
        raise Warning('HTMD Config file {} may not be a python file (Check HTMD_CONFIG).'.format(configfile))
    else:
        try:
            with open(configfile) as f:
                code = compile(f.read(), configfile, 'exec')
                exec(code)
        except:
            raise RuntimeError('Failed to execute the HTMD Config file {}.'.format(configfile))
        else:
            print('\nHTMD Config file {} executed.'.format(configfile))

print('Deprecation warning: "from htmd import *" will be deprecated. '
      '\nTo import all HTMD shortcuts please from now on use "from htmd.ui import *"')
