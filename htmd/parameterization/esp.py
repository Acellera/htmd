# (c) 2015-2017 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import os
from math import pi as PI
from math import sqrt, sin, cos, acos
import numpy as np
from numpy.random import uniform as rand
from scipy import constants as const
from scipy.spatial.distance import cdist
import nlopt
import unittest

from htmd.molecule.vdw import radiusByElement

class ESP:
    """
    Electrostatic potential (ESP) charge fitting

    Capabilities
    ------------
    - Generate points
    - Consider equivalent atoms
    - Impose total molecule charge
    - Impose boundaries for charge values
    - Charge values can be frozen

    The charges are fitting to reproduce ESP at the reference point computed by QM. The fitting is performed with
    COBYLA algorithm considering the equivalent atoms and imposing the total charge of the molecule.

    The charge values are confined to [-1.25; 1.25] ([0.001; 1.25] for hydrogen) range to prevent non-physical results.
    Also, the specific charges can be fixed to the orginal values (as defined in the molecule object).

    Attributes
    ----------
    molecule : FFMolecule
        Molecule object
    qm_results : List of QMResult
        Reference QM results
    fixed : list of ints
        List of fixed atom indices

    Examples
    --------

    Load water molecule
    >>> import os
    >>> from htmd.home import home
    >>> from htmd.parameterization.ffmolecule import FFMolecule, FFTypeMethod
    >>> molFile = os.path.join(home('test-qm'), 'H2O.mol2')
    >>> mol = FFMolecule(molFile, method=FFTypeMethod.GAFF2)

    Set up and run a QM (B3LYP/6-31G*) calculation of ESP
    >>> from htmd.qm import Psi4
    >>> from tempfile import mkdtemp
    >>> qm = Psi4()
    >>> qm.molecule = mol
    >>> qm.esp_points = ESP.generate_points(mol)[0]
    >>> qm.directory = mkdtemp()
    >>> qm_results = qm.run()
    >>> qm_results[0].errored
    False

    Create an ESP charge fitting object
    >>> from htmd.parameterization.esp import ESP
    >>> esp = ESP()
    >>> esp # doctest: +ELLIPSIS
    <htmd.parameterization.esp.ESP object at 0x...>

    Set up and run charge fitting
    >>> esp.molecule = mol
    >>> esp.qm_results = qm_results
    >>> esp_results = esp.run()

    ESP charges for water molecule
    >>> esp_results['charges'] # doctest: +ELLIPSIS
    array([-0.3940...,  0.1970...,  0.1970...])
    """

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
    def generate_points(molecule, vdw_radii=(1.4, 1.6, 1.8, 2.0, 2.2), density=10):
        """
        Generate points for ESP fitting around a molecule.

        The points are distributed on concentric spheres with centres at the atom positions. The radii of the spheres are
        the computed as van der Waals radii multiplied by a factor. The point are close to other atoms that the the radius
        are excluded.

        Parameters
        ----------
        molecule : FFMolecule
            Molecule object
        vdw_radii : list of floats
            List of van der Waals factors
        density: float
            Density of points

        Return
        ------
        points : list of list
            Set of points
        """
        points = []
        for frame in range(molecule.coords.shape[2]):
            pp = ESP._points(molecule.coords[:, :, frame],
                             ESP._vdw_radii(molecule.element),
                             vdw_radii, density)
            points.append(pp)

        return points

    def __init__(self):

        self.molecule = None
        self.qm_results = None
        self.fixed = []

        self._reciprocal_distances = None

    @property
    def ngroups(self):
        """Number of charge groups"""

        return len(self.molecule._equivalent_atom_groups)

    def _map_groups_to_atoms(self, group_charges):

        charges = np.zeros(self.molecule.numAtoms)
        for atom_group, group_charge in zip(self.molecule._equivalent_atom_groups, group_charges):
            charges[atom_group] = group_charge

        return charges

    def _compute_constraint(self, group_charges, _):

        charges = self._map_groups_to_atoms(group_charges)
        constraint = np.sum(charges) - self.molecule.netcharge

        return constraint

    def _get_bounds(self):

        # Set bound arrays
        lower_bounds = np.ones(self.ngroups) * -1.25
        upper_bounds = np.ones(self.ngroups) * +1.25

        # If the restraint relates to an H, set the lower bound to 0
        for i in range(self.ngroups):
            if 'H' == self.molecule.element[self.molecule._equivalent_atom_groups[i][0]]:
                lower_bounds[i] = 0.001

        # Fix the charges of the specified atoms to those already set in the
        # charge array. Note this also fixes the charges of the atoms in the
        # same equivalency group.
        for atom in self.fixed:
            group = self.molecule._equivalent_group_by_atom[atom]
            lower_bounds[group] = self.molecule.charge[atom]
            upper_bounds[group] = self.molecule.charge[atom]

        return lower_bounds, upper_bounds

    def _compute_objective(self, group_charges, _):

        qm_result = self.qm_results[0]

        # Compute the reciprocal distances between ESP points and atomic charges
        if self._reciprocal_distances is None:
            distances = cdist(qm_result.esp_points, self.molecule.coords[:, :, 0])
            distances *= const.physical_constants['Bohr radius'][0]/const.angstrom  # Angstrom --> Bohr
            self._reciprocal_distances = np.reciprocal(distances)

        charges = self._map_groups_to_atoms(group_charges)
        esp_values = np.dot(self._reciprocal_distances, charges)
        rms = np.sqrt(np.mean((esp_values - qm_result.esp_values)**2))

        return rms

    def run(self):
        """
        Run ESP charge fitting

        Return
        ------
        results : dict
            Dictionary with the fitted charges and fitting loss value
        """

        # Get charge bounds
        lower_bounds, upper_bounds = self._get_bounds()

        # Set up NLopt
        opt = nlopt.opt(nlopt.LN_COBYLA, self.ngroups)
        opt.set_min_objective(self._compute_objective)
        opt.set_lower_bounds(lower_bounds)
        opt.set_upper_bounds(upper_bounds)
        opt.add_equality_constraint(self._compute_constraint)
        opt.set_xtol_rel(1.e-6)
        opt.set_maxeval(1000*self.ngroups)
        opt.set_initial_step(0.001)

        # Optimize the charges
        group_charges = opt.optimize(np.zeros(self.ngroups) + 0.001) # TODO: a more elegant way to set initial charges
        # TODO: check optimizer status
        charges = self._map_groups_to_atoms(group_charges)
        loss = self._compute_objective(group_charges, None)

        return {'charges': charges, 'loss': loss}


class TestESP(unittest.TestCase):

    def setUp(self):

        from htmd.home import home
        from htmd.parameterization.ffmolecule import FFMolecule, FFTypeMethod

        molFile = os.path.join(home('test-param'), 'H2O2.mol2')
        self.mol = FFMolecule(molFile, method=FFTypeMethod.GAFF2)
        self.esp = ESP()
        self.esp.molecule = self.mol

    def test_ngroups(self):

        self.assertEqual(self.esp.ngroups, 2)

    def test_mapping(self):

        charges = self.esp._map_groups_to_atoms([1, 2])
        self.assertListEqual(list(charges), [1, 1, 2, 2])

    def test_constraint_function(self):

        self.assertEqual(self.esp._compute_constraint([1, 2], None), 6)

    def test_get_bounds(self):

        lower_bounds, upper_bounds = self.esp._get_bounds()
        self.assertEqual(list(lower_bounds), [-1.25, 0.001])
        self.assertEqual(list(upper_bounds), [1.25, 1.25])

        self.esp.fixed = [0]
        lower_bounds, upper_bounds = self.esp._get_bounds()
        self.assertEqual(list(lower_bounds), [-0.25279998779296875, 0.001])
        self.assertEqual(list(upper_bounds), [-0.25279998779296875, 1.25])

        self.esp.fixed = [3]
        lower_bounds, upper_bounds = self.esp._get_bounds()
        self.assertEqual(list(lower_bounds), [-1.25, 0.25279998779296875])
        self.assertEqual(list(upper_bounds), [1.25, 0.25279998779296875])

    def test_compute_objective(self):

        from htmd.qm import QMResult

        np.random.seed(20170901) # Make the test deterministic

        result = QMResult
        result.esp_points = np.random.normal(size=(100, 3))
        result.esp_values = np.random.normal(size=100)

        self.esp.qm_results = [result]

        rms = self.esp._compute_objective([1, 2], None)
        self.assertAlmostEqual(rms, 6.0746907953331517)

def load_tests(loader, tests, ignore):
    """Load DocTests into as a Unittest suite"""

    import doctest

    if os.environ.get('TRAVIS_OS_NAME') != 'osx':  # Psi4 does not work on Mac
        tests.addTests(doctest.DocTestSuite(__name__))

    return tests

if __name__ == '__main__':

    unittest.main(verbosity=2)
