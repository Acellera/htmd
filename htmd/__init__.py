# (c) 2015-2018 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
from htmd.version import version as _version
from htmd.versionwarnings import _issueWarnings, _disableWarnings
import __main__ as main
import os.path
from htmd.config import _config
import htmd.home
import logging.config

__version__ = _version()

try:
    logging.config.fileConfig(os.path.join(htmd.home.home(), 'logging.ini'), disable_existing_loggers=False)
except:
    print("HTMD: Logging setup failed")

reg_file = os.path.join(os.path.expanduser('~'), '.htmd', '.registered-htmd', 'registration')
if (not hasattr(main, '__file__')) and \
        ((not os.path.isfile(reg_file)) or os.getenv("LICENCE_ACCEPTED") == "YES" or os.getenv("TRAVIS_REPO_SLUG")):
    print('\nCopyright by Acellera Ltd. By executing you are accepting the License. In order to register, '
          'run htmd_register on your terminal.')
    print('The registration information must be valid so that it might be verified.')

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

_issueWarnings()