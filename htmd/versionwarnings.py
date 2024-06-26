# (c) 2015-2022 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import os
from pathlib import Path

_htmdconfigfolder = os.path.join(os.path.expanduser("~"), ".htmd")
_warningsfile = os.path.join(_htmdconfigfolder, ".disablewarnings")
try:
    os.makedirs(_htmdconfigfolder, exist_ok=True)
    if not os.path.exists(_warningsfile):
        Path(_warningsfile).touch()
except Exception:
    pass


def _issueWarnings():
    if not os.path.exists(_warningsfile):
        return  # Don't issue warnings if we failed to create the disable file

    if ("CI" in os.environ) and os.environ["CI"]:
        return  # Don't issue warnings if running in CI


def _disableWarnings(version):
    if not os.path.exists(_warningsfile):
        return

    with open(_warningsfile, "r") as f:
        versions = f.read().splitlines()

    if version not in versions:
        versions.append(version)

    with open(_warningsfile, "w") as f:
        f.write("\n".join(versions))
