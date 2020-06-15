import os
import warnings
from pathlib import Path

_htmdconfigfolder = os.path.join(os.path.expanduser("~"), ".htmd")
_warningsfile = os.path.join(_htmdconfigfolder, ".disablewarnings")
try:
    os.makedirs(_htmdconfigfolder, exist_ok=True)
    if not os.path.exists(_warningsfile):
        Path(_warningsfile).touch()
except:
    pass


def _issueWarnings():
    if not os.path.exists(_warningsfile):
        return  # Don't issue warnings if we failed to create the disable file

    if ("CI" in os.environ) and os.environ["CI"]:
        return  # Don't issue warnings if running in CI

    with open(_warningsfile, "r") as f:
        disabledversions = f.read().splitlines()

    if "1.16" not in disabledversions:
        warnings.warn(
            "As of HTMD 1.16 the default number of threads HTMD spawns for calculations is set to 1. "
            "You can enable parallelism at your own risk using `config(njobs=-2)` in the beginning of your scripts. "
            "To disable this warning run once: `from htmd import _disableWarnings; _disableWarnings('1.16');`",
            UserWarning,
        )

    if "1.21" not in disabledversions:
        warnings.warn(
            "As of HTMD 1.21 support for ACEMD v2 has stopped. Please use ACEMD3 instead"
            " as well as the corresponding equilibration and production protocols. "
            "To disable this warning run once: `from htmd import _disableWarnings; _disableWarnings('1.21');`",
            UserWarning,
        )


def _disableWarnings(version):
    if not os.path.exists(_warningsfile):
        return

    with open(_warningsfile, "r") as f:
        versions = f.read().splitlines()

    if version not in versions:
        versions.append(version)

    with open(_warningsfile, "w") as f:
        f.write("\n".join(versions))

