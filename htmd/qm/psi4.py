# (c) 2015-2018 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import os
import numpy as np
from scipy import constants as const

import periodictable as pt

from htmd.qm.base import QMBase, QMResult


class Psi4(QMBase):
    """
    Class to set up and run QM calculations with Psi4

    Capabilities
    ------------
    - Single-point energy calculations with HF and DFT.
    - Electronic properties: dipole and quadrupole monents, Mulliken charges, ESP at given points
    - Geometry optimization with/without dihedral restraints
    - Can use various queuing systems

    Attributes
    ----------
    molecule : :class:`FFMolecule`
        Molecule
    multiplity : int, default = 1
        Multiplicity of the molecule (1 = singlet, 2 = doublet, 3 = triplet, etc.)
    theory : str, default = 'BLYP3'
        Level of theory. A full list is at Psi4.THEORIES
    basis : str, default = '6-31G*'
        Basis set. A full list is at Psi4.BASIS_SETS
    correction : str, default = 'none'
        Empirical dispersion correction. A full list is at Psi4.CORRECTIONS
    solvent : str, default = 'vacuum'
         Implicit solvent model. A full isit is at Psi4.SOLVENTS
    esp_points : nunpy.ndarray, default = None
        Point coordinates to compute ESP values. The array shape has to be (number_of_points, 3)
    optimize : bool, default = False
        Run geometery optimization
    restrained_dihedrals : nunpy.ndarray, default = None
        List of restrained dihedrals (0-based atom indices). The array shape has to be (number_of_restrained_dihedrals, 4).
    queue : :class:`SimQueue`, defualt = :class:`LocalCPUQueue`
        Queue object used to run simulations
    directory : str, default = '.'
        Working directoy

    Examples
    --------

    Create an object of H2 molecule
    >>> import os
    >>> from htmd.home import home
    >>> from moleculekit.molecule import Molecule
    >>> molFile = os.path.join(home('test-qm'), 'H2-0.74.mol2')
    >>> mol = Molecule(molFile)

    Create a Psi4 object
    >>> from htmd.qm import Psi4
    >>> qm = Psi4()
    >>> qm # doctest: +ELLIPSIS
    <htmd.qm.psi4.Psi4 object at 0x...>

    Run single-point QM calculation of H2 with BLYP and cc-pVDZ
    >>> from tempfile import TemporaryDirectory
    >>> with TemporaryDirectory() as tmp:
    ...     qm.molecule = mol
    ...     qm.theory = 'BLYP'
    ...     qm.basis = 'cc-pVDZ'
    ...     qm.directory = tmp
    ...     result = qm.run()

    The QM results are returned as a list of htmd.qm.QMResult objects. See htmd.qm.QMResult documentation for details.
    >>> result # doctest: +ELLIPSIS
    [<htmd.qm.base.QMResult object at 0x...>]
    >>> result[0].errored
    False
    >>> result[0].energy # doctest: +ELLIPSIS
    -728.97083177...
    >>> result[0].mulliken # doctest: +ELLIPSIS
    [...0.0, ...0.0]

    Run the geometry optimization of H2 with BLYP, but change basis to 3-21G.
    NOTE: the `directory` attribut needs to be set to empty or non-existing directory, overwise the previous
    calculation results are read.
    >>> with TemporaryDirectory() as tmp:
    ...     qm.basis = '3-21G'
    ...     qm.optimize = True
    ...     qm.directory = tmp
    ...     result = qm.run()

    The initial coordinates and optimize coordinates of H2 can be compared
    >>> mol.coords
    array([[[ 0.  ],
            [ 0.  ],
            [-0.37]],
    <BLANKLINE>
           [[ 0.  ],
            [ 0.  ],
            [ 0.37]]], dtype=float32)
    >>> result[0].coords
    array([[[ 0.        ],
            [ 0.        ],
            [-0.37527978]],
    <BLANKLINE>
           [[ 0.        ],
            [ 0.        ],
            [ 0.37527978]]])

    The QM calculations run using LocalCPUQueue by default, but this can be changed to the others.
    >>> from htmd.queues.slurmqueue import SlurmQueue
    >>> qm.queue = SlurmQueue() # doctest: +SKIP
    """

    @property
    def _command(self):
        return 'psi4 -i psi4.in -o psi4.out &> psi4.log'

    def _completed(self, directory):
        # Abuse "timer.dat" to detect if a Psi4 job has completed.
        # "timer.dat" is written up on completion (successful or failed),
        # but it will be missing a job has benn killed or crashed due to a machine problems.
        return os.path.exists(os.path.join(directory, 'timer.dat'))

    def _writeInput(self, directory, iframe):

        with open(os.path.join(directory, 'psi4.in'), 'w') as f:

            f.write('import psi4\n\n')

            f.write('set_num_threads(%d)\n' % self.queue.ncpu)
            f.write('set_memory(\'{} MB\')\n\n'.format(self.queue.memory))

            # Use more conservative memory estimate, otherwise Psi4 exceeds the memory limit of a queuing system
            f.write('set { scf_mem_safety_factor 0.7 }\n\n')

            reference = 'r' if self.multiplicity == 1 else 'u'
            reference += 'hf' if self.theory == 'HF' else 'ks'
            f.write('set { reference %s }\n\n' % reference)

            # Set basis sets
            atomic_number = lambda element: pt.elements.symbol(element).number
            elements = sorted(np.unique(self._molecule.element), key=atomic_number)
            element_basis = [self.substituteBasisSet(element, self.basis) for element in elements]
            f.write('basis = \'\'\n')
            for element, basis in zip(elements, element_basis):
                f.write(f'basis += \'assign {element} {basis}\'\n')
            f.write('psi4.basis_helper(basis)\n\n')

            if self.solvent == 'vacuum':
                pass
            elif self.solvent == 'PCM':
                # TODO check all agruments
                f.write('set { PCM true\n      PCM_scf_type total }\n')
                f.write('PCM = {\n')
                f.write('  Units = Angstrom\n')
                f.write('  Medium {\n    SolverType = IEFPCM\n    Solvent = Water\n  }\n')
                f.write('  Cavity {\n    RadiiSet = UFF\n    Type = GePol\n')
                f.write('    Scaling = False\n    Area = 0.3\n    Mode = Implicit\n  }\n')
                f.write('}\n\n')
            else:
                raise NotImplementedError

            # Write the molecule
            f.write('molecule MOL {\n' )
            f.write('    %d %d\n' % (self.charge, self.multiplicity))
            f.write('    noreorient\n')
            f.write('    nocom\n')
            f.write('    symmetry c1\n')
            elements = self._molecule.element
            coords = self._molecule.coords[:, :, iframe]
            for element, coord in zip(elements, coords):
                f.write('    %-2s %10f %10f %10f\n' % (element, coord[0], coord[1], coord[2]))
            f.write('}\n\n')

            if self._restrained_dihedrals is not None:
                dihedrals = ['%d %d %d %d' % tuple(dihedral) for dihedral in self._restrained_dihedrals]
                dihedrals = ', '.join(dihedrals)
                f.write('set optking { frozen_dihedral = (" %s ") }\n\n' % dihedrals)

            # Enable a dynamic optimization algorithm selection to converge problematic cases:
            # http://www.psicode.org/psi4manual/master/optking.html#dealing-with-problematic-optimizations
            if self.optimize:
                f.write('set optking { dynamic_level = 1 \n geom_maxiter = 250\n print_trajectory_xyz_file = True }\n\n')

            function = 'optimize' if self.optimize else 'energy'
            theory = 'SCF' if self.theory == 'HF' else self.theory
            theory += '' if self.correction == 'none' else '-%s' % self.correction
            f.write('energy, wfn = %s(\'%s\', return_wfn=True)\n\n' % (function, theory))

            f.write('oeprop(wfn, \'DIPOLE\', \'QUADRUPOLE\', \'MULLIKEN_CHARGES\')\n')
            if self.esp_points is not None:
                f.write('oeprop(wfn, \'GRID_ESP\')\n')
            f.write('\n')

            f.write('with open(\'psi4out.xyz\', \'w\') as f:\n')
            f.write('    f.write(\'%d \' )\n' % len(coords))
            f.write('    f.write(\'%.12f\\n\' % energy)\n')
            f.write('    f.write(MOL.save_string_xyz())\n')

    def _readOutput(self, directory):

        result = QMResult()
        result.completed = True

        xyzFile = os.path.join(directory, 'psi4out.xyz')
        if os.path.exists(xyzFile):
            with open(xyzFile) as f:
                result.energy = float(f.readline().split()[1])  # Read the 2nd number on the 1st line
            result.energy *= const.physical_constants['Hartree energy'][0]/(const.kilo*const.calorie/const.Avogadro) # Hartree to kcal/mol
            result.coords = np.loadtxt(xyzFile, skiprows=2, usecols=(1, 2, 3))
            result.coords = np.atleast_3d(result.coords)  # TODO get rid of this
        else:
            result.errored = True

        if self.esp_points is not None:
            espFile = os.path.join(directory, 'grid_esp.dat')
            if os.path.exists(espFile):
                result.esp_points = self.esp_points
                result.esp_values = np.loadtxt(espFile)
                result.esp_values *= const.angstrom/const.physical_constants['Bohr radius'][0] # 1/Bohr to 1/Angstrom
            else:
                result.errored = True

        outFile = os.path.join(directory, 'psi4.out')
        if os.path.exists(outFile):
            with open(outFile) as f:
                lines = f.readlines()

            for i in range(len(lines)):
                if lines[i].strip().startswith('Mulliken Charges:'):
                    result.mulliken = []
                    for j in range(self._natoms):
                        s = lines[i + 2 + j].split()
                        result.mulliken.append(float(s[5]))

                if lines[i].strip().startswith('Dipole Moment:'):
                    s = lines[i + 1].split()
                    result.dipole = [float(s[1]), float(s[3]), float(s[5]), float(s[7])]

                if lines[i].strip().startswith('Traceless Quadrupole Moment:'):
                    s1, s2 = lines[i + 1].split(), lines[i + 2].split()
                    result.quadrupole = [float(s1[1]), float(s1[3]), float(s1[5]),
                                         float(s2[1]), float(s2[3]), float(s2[5])]

        else:
            result.errored = True

        return result


if __name__ == '__main__':

    import sys
    import doctest

    if os.environ.get('TRAVIS_OS_NAME') != 'osx':  # Psi4 does not work in Mac
        if doctest.testmod().failed:
            sys.exit(1)
