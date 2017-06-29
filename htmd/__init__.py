from htmd.version import version as _version
from htmdx.cli import check_registration, show_news
from htmd.latest import compareVersions
import os
from htmd.config import config

config()
__version__ = _version()

if not (os.getenv("HTMD_NONINTERACTIVE")):
    check_registration(product='htmd')
    show_news()
    compareVersions()

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

print('Deprecation warning: To import all HTMD shortcuts please use "from htmd.ui import *"')