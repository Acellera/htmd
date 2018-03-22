# (c) 2015-2018 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
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
