# (c) 2015-2018 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
from pint import UnitRegistry
import numpy as np


def convert(fromunit, tounit, value, fstep=1, timestep=4):
    """ Converts values between units

    Parameters
    ----------
    fromunit : str
        Unit to convert from
    tounit : str
        Unit to convert to
    value : scalar
        The value to convert
    fstep : float
        The sampling frame step of the simulation in nanoseconds
    timestep : int
        The timestep of the simulation in femtoseconds

    Returns
    -------
    conv : scalra
        The converted value
    """

    ureg = UnitRegistry()
    ureg.define('frame = {} * ns'.format(fstep))
    ureg.define('step = ({} / 1000000) * ns = timestep'.format(timestep))

    q = ureg.Quantity(value, fromunit)
    convval = q.to(tounit)
    if convval.units == 'frame' or convval.units == 'step':
        vals = np.round(convval.magnitude).astype(int)
        if vals.size == 1:  # Fix for PyEMMA tica. remove in future
            return int(vals)
        return vals
    else:
        return convval.magnitude
