# (c) 2015-2022 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
from pint import UnitRegistry
import numpy as np


def convert(
    fromunit: str,
    tounit: str,
    value: float | list | range | np.ndarray,
    fstep: float = 1,
    timestep: float = 4,
):
    """Convert values between units.

    Parameters
    ----------
    fromunit : str
        Unit to convert from.
    tounit : str
        Unit to convert to.
    value : scalar or list or range or numpy.ndarray
        The value or array of values to convert.
    fstep : float, optional
        Sampling frame step of the simulation in nanoseconds.
    timestep : float, optional
        Timestep of the simulation in femtoseconds.

    Returns
    -------
    conv : scalar or numpy.ndarray
        The converted value. When converting to ``frame`` or ``step`` units,
        the result is rounded to the nearest integer.

    Examples
    --------
    >>> convert("ns", "frame", 100, fstep=0.1)
    1000
    """

    ureg = UnitRegistry()
    ureg.define("frame = {} * ns".format(fstep))
    ureg.define("step = ({} / 1000000) * ns = timestep".format(timestep))

    q = ureg.Quantity(value, fromunit)
    convval = q.to(tounit)
    if convval.units == "frame" or convval.units == "step":
        vals = np.round(convval.magnitude).astype(int)
        if vals.size == 1:  # Fix for tica. remove in future
            return int(vals)
        return vals
    else:
        return convval.magnitude
