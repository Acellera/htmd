import os
import warnings
from pathlib import Path

try:
    _htmdconfigfolder = os.path.join(os.path.expanduser("~"), '.htmd')
    os.makedirs(_htmdconfigfolder, exist_ok=True)
    _warningsfile = os.path.join(_htmdconfigfolder, '.disablewarnings')
    if not os.path.exists(_warningsfile):
        Path(_warningsfile).touch()
except:
    pass

def _issueWarnings():
    if not os.path.exists(_warningsfile):
        return

    with open(_warningsfile, 'r') as f:
        disabledversions = f.read().splitlines()

    if '1.16' not in disabledversions:
        warnings.warn('As of HTMD 1.16 the default ACEMD version for all protocols has changed to version 3. ' \
                    'If you want to use version 2 protocols change the _version argument in the protocols ' \
                    'or add `config(acemdversion=2)` to the beginning of your scripts. ' \
                    'To disable this warning run once: `from htmd import _disableWarnings; _disableWarnings(\'1.16\');`'
        , UserWarning)

    if '1.17' not in disabledversions:
        warnings.warn('As of HTMD 1.17 the default number of threads HTMD spawns for calculations is set to 1. ' \
                    'You can enable parallelism at your own risk using `config(njobs=-2)` in the beginning of your scripts. ' \
                    'To disable this warning run once: `from htmd import _disableWarnings; _disableWarnings(\'1.17\');`'
        , UserWarning)


def _disableWarnings(version):
    if not os.path.exists(_warningsfile):
        return

    with open(_warningsfile, 'r') as f:
        versions = f.read().splitlines()

    if version not in versions:
        versions.append(version)

    with open(_warningsfile, 'w') as f:
        f.write('\n'.join(versions))