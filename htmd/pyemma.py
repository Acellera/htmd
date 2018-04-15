# This file is necessary to remove the version warning of pyEMMA which should not appear in HTMD
import warnings

with warnings.catch_warnings():
    warnings.simplefilter('ignore', category=UserWarning)
    from pyemma import plots, msm, coordinates

    import sys
    sys.modules['htmd.pyemma.msm'] = msm
    sys.modules['htmd.pyemma.plots'] = plots
    sys.modules['htmd.pyemma.coordinates'] = coordinates