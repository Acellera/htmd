from htmd.ui import *  # TODO: This will be deprecated in a few versions
from htmd.version import version as _version
import os
from htmd.config import _config

__version__ = _version()


if _config['configfile']:
    if not os.path.isfile(_config['configfile']):
        raise FileNotFoundError('HTMD Config file {} does not exist.'.format(_config['configfile']))
    elif not _config['configfile'].endswith('.py'):
        raise Warning('HTMD Config file {} may not be a python file.'.format(_config['configfile']))
    else:
        try:
            with open(_config['configfile']) as f:
                exec(compile(f.read(), _config['configfile'], 'exec'))
        except:
            raise RuntimeError('Failed to execute the HTMD Config file {}.'.format(_config['configfile']))
        else:
            print('\nHTMD Config file {} executed.'.format(_config['configfile']))

print('Deprecation warning: "from htmd import *" will be deprecated. '
      '\nTo import all HTMD shortcuts please from now on use "from htmd.ui import *"')
