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
    fstep : int
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
        return np.round(convval.magnitude).astype(int)
    else:
        return convval.magnitude
