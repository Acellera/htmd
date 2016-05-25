from pint import UnitRegistry
import numpy as np


def convert(fromunit, tounit, value, fstep=1):
    ureg = UnitRegistry()
    ureg.define('frame = {} * ns = step'.format(fstep))

    q = ureg.Quantity(value, fromunit)
    convval = q.to(tounit)
    if convval.units == 'frame':
        return np.round(convval.magnitude).astype(int)
    else:
        return convval.magnitude
