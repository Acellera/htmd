# (c) 2015-2018 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import logging
import os
import unittest

import numpy as np
from scipy import constants as const
from scipy.spatial.distance import cdist
import nlopt

from moleculekit.molecule import Molecule
from moleculekit.periodictable import periodictable
from htmd.parameterization.detect import detectEquivalentAtoms

logger = logging.getLogger(__name__)


def randomPointsOnSphere(num_points):
    """
    Generate random points on a sphere

    >>> np.random.seed(20181113)
    >>> randomPointsOnSphere(2)
    array([[-0.8483968 , -0.52844579, -0.03111143],
           [ 0.37728796, -0.925473  ,  0.03396356]])
    """
    points = np.random.normal(size=(num_points, 3))
    points = np.where(points == [0, 0, 0], [1, 0, 0], points)

    points /= np.linalg.norm(points, axis=1, keepdims=True)

    assert points.shape == (num_points, 3)
    assert np.allclose(np.linalg.norm(points, axis=1), 1)

    return points

class MoleculeGrid:
    """
    Molecular grid for RESP charge fitting

    The grid points are distributed on concentric spheres (shells) with centres at each atom position.
    The radii of the spheres are the computed as van der Waals radii multiplied by a shell factor.
    The point that are closer to other atoms than the the smallest shell are excluded.

    Parameters
    ----------
    molecule : FFMolecule
        Molecule object
    shell_factors : list of floats
        List of van der Waals factors for each shell
    density: float
        Density of points (in points/Ang**2)

    Examples
    --------

    Load water molecule
    >>> import os
    >>> from htmd.home import home
    >>> from moleculekit.molecule import Molecule
    >>> molFile = os.path.join(home('test-charge'), 'H2O.mol2')
    >>> mol = Molecule(molFile, guessNE='bonds', guess=('angles', 'dihedrals'))
    >>> mol.write('H2O.xyz') # doctest: +SKIP

    >>> np.random.seed(20181113)
    >>> grid = MoleculeGrid(mol)
    >>> len(grid.getPoints())
    4157
    >>> grid.getPoints() # doctest: +NORMALIZE_WHITESPACE
    array([[-0.88818835, -1.04753265, -0.02590512],
           [-0.03612161, -0.95868347,  1.63620777],
           [ 0.62043862,  2.06506096, -0.65822007],
           ...,
           [ 1.79185447, -2.69472971,  0.04730653],
           [ 0.15053   , -2.95065264,  0.56205264],
           [ 1.34573805, -1.44031174, -1.75836006]])
    >>> grid.writeXYZ('H2O_default.xyz') # doctest: +SKIP

    >>> np.random.seed(20181113)
    >>> grid = MoleculeGrid(mol, shell_factors=(1, 2), density=50)
    >>> len(grid.getPoints())
    4275
    >>> grid.getPoints() # doctest: +NORMALIZE_WHITESPACE
    array([[ 0.70522759,  1.49704354, -0.45864291],
           [ 1.63974354,  1.22579616, -0.64424237],
           [ 1.74152223,  0.47726292, -1.17241699],
           ...,
           [-0.70728694, -2.44544594,  0.88320856],
           [ 0.59560643, -2.38238637,  1.79258897],
           [ 2.11277479, -1.56282107, -1.2258255 ]])
    >>> grid.writeXYZ('H2O_1_2__50.xyz') # doctest: +SKIP
    """
    def __init__(self, molecule, shell_factors=(1.4, 1.6, 1.8, 2.0), density=25):

        self._molecule = molecule
        if not isinstance(self._molecule, Molecule):
            raise TypeError('"molecule" has to be instance of {}'.format(Molecule))
        if self._molecule.numFrames != 1:
            raise ValueError('"molecule" can have just one frame, but it has {}'.format(self._molecule.numFrames))

        self._shell_factors = list(map(float, shell_factors))
        for factor in self._shell_factors:
            if factor <= 0:
                raise ValueError('The elements of "shell_factors" have to be positive, but get {}'.format(factor))

        self._density = density
        if self._density <= 0:
            raise ValueError('"density" has to be positive, but get {}'.format(self._density))

        self._points = self._generatePoints()
        self._points = self._filterPoints(self._points)

    def _generatePoints(self):

        all_points = []

        for element, coord in zip(self._molecule.element, self._molecule.coords[:, :, 0]):
            vdw_radius = periodictable[element].vdw_radius
            for factor in self._shell_factors:

                # Compute the number of point for each shell
                radius = factor * vdw_radius
                area = 4/3 * np.pi * radius**2
                num_points = int(self._density * area)

                # Generate points
                points = radius * randomPointsOnSphere(num_points) + coord
                all_points.append(points)

        return np.concatenate(all_points)

    def _filterPoints(self, points):

        # Compute distance threshold for each atom
        thresholds = np.array([periodictable[element].vdw_radius for element in self._molecule.element])
        thresholds *= min(self._shell_factors) - 0.001

        # Detect the points further away for each atom than its threshold
        distances = cdist(points, self._molecule.coords[:, :, 0])
        is_valid = np.all(distances > thresholds, axis=1)

        return points[is_valid]

    def getPoints(self):
        """
        Return grid points
        """
        return self._points

    def writeXYZ(self, file):
        """
        Write the molecular grid in XYZ format
        """
        if isinstance(file, str):
            with open(file, 'w') as stream:
                self.writeXYZ(stream)

        else:
            points = self.getPoints()
            file.write('{}\n\n'.format(len(points)))
            for point in points:
                file.write('X {:10.6f} {:10.6f} {:10.6f}\n'.format(*tuple(point)))


class ESP:
    """
    Electrostatic potential (ESP) charge fitting

    Capabilities
    ------------
    - Consider equivalent atoms
    - Impose total molecule charge
    - Impose boundaries for charge values
    - Charge values can be frozen

    The charges are fitting to reproduce ESP at the reference point computed by QM. The fitting is performed with
    COBYLA algorithm considering the equivalent atoms and imposing the total charge of the molecule.

    The charge values are confined to [-1.25; 1.25] ([0.0; 1.25] for hydrogen) range to prevent non-physical results.
    Also, the specific charges can be fixed to the orginal values (as defined in the molecule object).

    Attributes
    ----------
    molecule : FFMolecule
        Molecule object
    qm_results : List of QMResult
        Reference QM results
    apply_bounds: boolean
        Apply bounds to atomic charges
    restraint_factor: float
        Restraint factor for heavy elements
    fixed : list of ints
        List of fixed atom indices

    Examples
    --------

    Load water molecule
    >>> import os
    >>> from htmd.home import home
    >>> from moleculekit.molecule import Molecule
    >>> molFile = os.path.join(home('test-charge'), 'H2O.mol2')
    >>> mol = Molecule(molFile, guessNE='bonds', guess=('angles', 'dihedrals'))

    Generate points
    >>> np.random.seed(20181113)
    >>> grid = MoleculeGrid(mol)
    >>> len(grid.getPoints())
    4157

    Set up and run a QM (B3LYP/6-31G*) calculation of ESP
    >>> from htmd.qm import Psi4
    >>> from tempfile import mkdtemp
    >>> qm = Psi4()
    >>> qm.molecule = mol
    >>> qm.esp_points = grid.getPoints()
    >>> qm.directory = mkdtemp()
    >>> qm_results = qm.run()
    >>> qm_results[0].errored
    False

    Fit ESP charges
    >>> from htmd.charge.esp import ESP
    >>> esp = ESP()
    >>> esp # doctest: +ELLIPSIS
    <htmd.charge.esp.ESP object at 0x...>
    >>> esp.molecule = mol
    >>> esp.qm_results = qm_results
    >>> esp_results = esp.run()
    >>> esp_results['charges'] # doctest: +ELLIPSIS
    array([-0.3936...,  0.1968...,  0.1968...])
    >>> esp_results['loss'] # doctest: +ELLIPSIS
    1.835...e-05
    >>> esp_results['RMSD'] # doctest: +ELLIPSIS
    0.004284...

    >>> esp = ESP()
    >>> esp.molecule = mol
    >>> esp.qm_results = qm_results
    >>> esp.restraint_factor = 0.001
    >>> esp_results = esp.run()
    >>> esp_results['charges'] # doctest: +ELLIPSIS
    array([-0.3811...,  0.1905...,  0.1905...])
    >>> esp_results['loss'] # doctest: +ELLIPSIS
    6.836...e-05
    >>> esp_results['RMSD'] # doctest: +ELLIPSIS
    0.004465...
    """
    def __init__(self):

        self.molecule = None
        self.qm_results = None
        self.apply_bounds = True
        self.restraint_factor = 0
        self.fixed = []

        self._molecular_charge = 0

        self._reciprocal_distances = None

        self._equivalent_atom_groups = None
        self._equivalent_group_by_atom = None

        self._restraint_factors = None

    @property
    def ngroups(self):
        """Number of charge groups"""
        return len(self._equivalent_atom_groups)

    def _map_groups_to_atoms(self, group_charges):

        charges = np.zeros(self.molecule.numAtoms)
        for atom_group, group_charge in zip(self._equivalent_atom_groups, group_charges):
            charges[list(atom_group)] = group_charge

        return charges

    def _constraint(self, group_charges, _):

        charges = self._map_groups_to_atoms(group_charges)
        constraint = np.sum(charges) - self._molecular_charge

        return constraint

    def _get_bounds(self):

        # Set very loose bounds, i.e. effectively no bounds
        lower_bounds = -10 * np.ones(self.ngroups)
        upper_bounds = +10 * np.ones(self.ngroups)

        if self.apply_bounds:

            # Set reasonable bounds
            lower_bounds = np.ones(self.ngroups) * -1.25
            upper_bounds = np.ones(self.ngroups) * +1.25

            # Bond hydrogen charges to be positive
            for i in range(self.ngroups):
                element = self.molecule.element[self._equivalent_atom_groups[i][0]]
                if element == 'H':
                    lower_bounds[i] = 0

        # Fix atom charges considering equivalent groups
        for atom in self.fixed:
            group = self._equivalent_group_by_atom[atom]
            lower_bounds[group] = self.molecule.charge[atom]
            upper_bounds[group] = self.molecule.charge[atom]

        return lower_bounds, upper_bounds

    def _objective(self, group_charges, _):

        points = self.qm_results[0].esp_points
        target_values = self.qm_results[0].esp_values

        # Compute the reciprocal distances between ESP points and atomic charges
        if self._reciprocal_distances is None:
            distances = cdist(points, self.molecule.coords[:, :, 0])
            distances *= const.physical_constants['Bohr radius'][0]/const.angstrom  # Angstrom --> Bohr
            self._reciprocal_distances = np.reciprocal(distances)

        charges = self._map_groups_to_atoms(group_charges)
        actual_values = np.dot(self._reciprocal_distances, charges)
        loss = np.mean((actual_values - target_values)**2) + np.mean(self._restraint_factors * charges**2)

        return loss

    def run(self):
        """
        Run ESP charge fitting

        Return
        ------
        results : dict
            Dictionary with the fitted charges and fitting loss value
        """
        logger.info('Start RESP charge fitting')

        self._molecular_charge = self.qm_results[0].charge

        # Detect equivalent atoms
        equivalents = detectEquivalentAtoms(self.molecule)
        self._equivalent_atom_groups = equivalents[0]
        self._equivalent_group_by_atom = equivalents[2]

        # Set up heavy atom restrains
        self._restraint_factors = np.zeros(self.molecule.numAtoms)
        if self.restraint_factor > 0:
            for i, element in enumerate(self.molecule.element):
                if element != 'H':
                    self._restraint_factors[i] = self.restraint_factor

        # Get charge bounds
        lower_bounds, upper_bounds = self._get_bounds()

        logger.info('Atom charge boundaries and restraint factor:')
        for i, name in enumerate(self.molecule.name):
            lower = lower_bounds[self._equivalent_group_by_atom[i]]
            upper = upper_bounds[self._equivalent_group_by_atom[i]]
            factor = self._restraint_factors[i]
            logger.info('  {:4s}: {:7.3f} {:7.3f} {:10.6f}'.format(name, lower, upper, factor))

        # Set up optimizer
        opt = nlopt.opt(nlopt.LN_COBYLA, self.ngroups)
        logger.info('Optimizer: {}'.format(opt.get_algorithm_name()))
        opt.set_min_objective(self._objective)
        opt.add_equality_constraint(self._constraint)
        logger.info('Molecular charges constraint: {:.3f}'.format(self._molecular_charge))
        opt.set_lower_bounds(lower_bounds)
        opt.set_upper_bounds(upper_bounds)
        opt.set_xtol_rel(1.e-6)
        opt.set_maxeval(1000*self.ngroups)
        opt.set_initial_step(0.001)

        # Optimize the charges
        group_charges = opt.optimize(np.zeros(self.ngroups))
        status = opt.last_optimize_result()
        loss = self._objective(group_charges, None)
        charges = self._map_groups_to_atoms(group_charges)
        logger.info('Optimizer status: {}'.format(status))
        logger.info('Final loss: {:.6f}'.format(loss))

        # Compute RMSD
        self._restraint_factors[:] = 0
        msd = self._objective(group_charges, None)
        logger.info('Final RMSD: {:.6f} au'.format(np.sqrt(msd)))

        logger.info('Finish RESP charge fitting')

        return {'charges': charges, 'status': status, 'loss': loss, 'RMSD': np.sqrt(msd)}


class TestESP(unittest.TestCase):

    def setUp(self):
        from htmd.home import home
        from moleculekit.molecule import Molecule

        molFile = os.path.join(home('test-param'), 'H2O2.mol2')
        mol = Molecule(molFile, guessNE='bonds', guess=('angles', 'dihedrals'))
        self.mol = mol
        self.esp = ESP()
        self.esp.molecule = self.mol

        # Precalculating here for the tests
        equivalents = detectEquivalentAtoms(self.esp.molecule)
        self.esp._equivalent_atom_groups = equivalents[0]
        self.esp._equivalent_group_by_atom = equivalents[2]

    def test_ngroups(self):
        self.assertEqual(self.esp.ngroups, 2)

    def test_mapping(self):

        charges = self.esp._map_groups_to_atoms([1, 2])
        self.assertListEqual(list(charges), [1, 1, 2, 2])

    def test_constraint(self):

        self.assertEqual(self.esp._constraint([1, 2], None), 6)

    def test_get_bounds(self):

        lower_bounds, upper_bounds = self.esp._get_bounds()
        self.assertEqual(list(lower_bounds), [-1.25, 0.0])
        self.assertEqual(list(upper_bounds), [1.25, 1.25])

        self.esp.apply_bounds = False
        self.esp.fixed = []
        lower_bounds, upper_bounds = self.esp._get_bounds()
        self.assertEqual(list(lower_bounds), [-10, -10])
        self.assertEqual(list(upper_bounds), [10, 10])

        self.esp.apply_bounds = True
        self.esp.fixed = [0]
        lower_bounds, upper_bounds = self.esp._get_bounds()
        self.assertEqual(list(lower_bounds), [-0.25279998779296875, 0.0])
        self.assertEqual(list(upper_bounds), [-0.25279998779296875, 1.25])

        self.esp.apply_bounds = True
        self.esp.fixed = [3]
        lower_bounds, upper_bounds = self.esp._get_bounds()
        self.assertEqual(list(lower_bounds), [-1.25, 0.25279998779296875])
        self.assertEqual(list(upper_bounds), [1.25, 0.25279998779296875])

    def test_objective(self):

        from htmd.qm import QMResult

        np.random.seed(20170901)  # Make the test deterministic

        result = QMResult
        result.esp_points = np.random.normal(size=(100, 3))
        result.esp_values = np.random.normal(size=100)

        self.esp.qm_results = [result]

        self.esp._restraint_factors = [0, 0, 0, 0]
        loss = self.esp._objective([1, 2], None)
        self.assertAlmostEqual(loss, 36.90186825890532)

        self.esp._restraint_factors = [1, 1, 0, 0]
        loss = self.esp._objective([1, 2], None)
        self.assertAlmostEqual(loss, 37.40186825890532)


def load_tests(loader, tests, ignore):
    """Load DocTests into as a Unittest suite"""
    import doctest

    tests.addTests(doctest.DocTestSuite(__name__))

    return tests


if __name__ == '__main__':

    unittest.main(verbosity=2)
