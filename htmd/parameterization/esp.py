# (c) 2015-2017 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#

from math import pi as PI
from math import sqrt, sin, cos, acos
import numpy as np
from numpy.random import uniform as rand
from htmd.molecule.vdw import radiusByElement

class ESP:
    # A temporary class to collect ESP stuff

    @staticmethod
    def _dist(a, b):
        c = a - b
        return sqrt(c.dot(c))

    @staticmethod
    def _dist2(a, b):
        c = a - b
        return c.dot(c)

    @staticmethod
    def _rand_sphere_sample(centre, r, density):
        # Produce a set of points on the sphere of radius r centred on centre
        # with ~density points / unit^2

        surface_area = 4. / 3. * PI * r * r
        n_points = int(density * surface_area)
        area_per_point = 1. / density  # surface_area / n_points
        mindist = sqrt(area_per_point / PI)

        points = np.zeros((n_points, 3))

        i = 0
        mindist2 = mindist * mindist
        pos = np.zeros(3)
        while i < n_points:
            z = 2. * rand() - 1.
            lon = 2. * PI * rand()
            lat = acos(z)
            x = cos(lon) * sin(lat)
            y = sin(lon) * sin(lat)

            pos[0] = x * r
            pos[1] = y * r
            pos[2] = z * r

            # Crudely test to see if it is in range of other points
            too_close = False
            for j in range(i):
                if ESP._dist2(points[j, :], pos) < mindist2:
                    too_close = True
                    break
            if not too_close:
                points[i, :] = pos
                i += 1
        points = points
        points = points + centre
        return points

    @staticmethod
    def _vdw_radii(elements):
        radii = np.zeros(elements.shape[0], dtype=np.float32)
        i = 0
        for e in elements:
            radii[i] = radiusByElement(e)
            i += 1
        return radii

    @staticmethod
    def _points(coords, radii, multipliers, density):
        points = []
        np.random.seed(0)
        # Make a set of points in a vdw shell around each atom
        for m in multipliers:
            for i in range(coords.shape[0]):
                p = ESP._rand_sphere_sample(coords[i, :], radii[i] * m, density)
                # remove any points that are within radii[i]*m of i-th atom
                for pp in p:
                    too_close = False
                    for j in range(coords.shape[0]):
                        if ESP._dist(coords[j, :], pp) < radii[j] * m:
                            too_close = True
                            break
                    if not too_close:
                        points.append(pp)

        return np.asarray(points, dtype=np.float32)

    @staticmethod
    def _generate_points(molecule, vdw_radii=[1.4, 1.6, 1.8, 2.0, 2.2], density=10):

        points = []
        for frame in range(molecule.coords.shape[2]):
            pp = ESP._points(molecule.coords[:, :, frame],
                             ESP._vdw_radii(molecule.element),
                             vdw_radii, density)
            points.append(pp)

        return points


if __name__ == '__main__':
    pass
