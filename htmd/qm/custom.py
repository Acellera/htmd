# (c) 2015-2017 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import os
import time
import pickle
import logging
import abc

import numpy as np
import nlopt

from htmd.numbautil import dihedralAngle
from htmd.qm.base import QMBase, QMResult
from protocolinterface import val


logger = logging.getLogger(__name__)


class Minimizer(abc.ABC):
    def __init__(self):
        pass

    @abc.abstractmethod
    def minimize(self, coords, restrained_dihedrals):
        pass


class OMMMinimizer(Minimizer):
    def __init__(self, mol, prm, platform='CPU', device=0, buildff='AMBER', guessAnglesDihedrals=True):
        """ A minimizer based on OpenMM

        Parameters
        ----------
        mol : Molecule
            The Molecule object containing the topology of the molecule
        prm : parmed.ParameterSet
            A parmed ParameterSet object containing the parameters of the molecule
        platform : str
            The platform on which to run the minimization ('CPU', 'CUDA')
        device : int
            If platform is 'CUDA' this defines which GPU device to use
        buildff : str
            The forcefield for which to build the Molecule to then minimize it with OpenMM
        guessAnglesDihedrals : bool
            If the class should guess angles and dihedrals of the Molecule.

        Examples
        --------
        >>> from htmd.parameterization.fftype import fftype
        >>> from htmd.molecule.molecule import Molecule

        >>> molFile = os.path.join(home('test-qm'), 'H2O2-90.mol2')
        >>> mol = Molecule(molFile)
        >>> prm, mol = fftype(mol, method='GAFF2')
        >>> mini = OMMMinimizer(mol, prm)
        >>> minimcoor = mini.minimize(mol.coords, restrained_dihedrals=[0, 1, 6, 12])
        """
        super().__init__()

        import parmed
        from htmd.util import tempname
        import simtk.openmm as mm

        if buildff == 'AMBER':
            self.structure = self._get_prmtop(mol, prm)

        self.system = self.structure.createSystem()
        self.platform = mm.Platform.getPlatformByName(platform)
        self.platprop = {'CudaPrecision': 'mixed', 'CudaDeviceIndex': device} if platform == 'CUDA' else None


    def _get_prmtop(self, mol, prm):
        from htmd.parameterization.writers import writeFRCMOD, getAtomTypeMapping
        from tempfile import TemporaryDirectory
        from subprocess import call
        from simtk.openmm import app

        with TemporaryDirectory() as tmpDir:
            frcFile = os.path.join(tmpDir, 'mol.frcmod')
            mapping = getAtomTypeMapping(prm)
            writeFRCMOD(mol, prm, frcFile, typemap=mapping)
            mol2 = mol.copy()
            mol2.atomtype[:] = np.vectorize(mapping.get)(mol2.atomtype)
            molFile = os.path.join(tmpDir, 'mol.mol2')
            mol2.write(molFile)

            with open(os.path.join(tmpDir, 'tleap.inp'), 'w') as file:
                file.writelines(('loadAmberParams %s\n' % frcFile,
                                 'MOL = loadMol2 %s\n' % molFile,
                                 'saveAmberParm MOL mol.prmtop mol.inpcrd\n',
                                 'quit'))

            with open(os.path.join(tmpDir, 'tleap.out'), 'w') as out:
                call(('tleap', '-f', 'tleap.inp'), cwd=tmpDir, stdout=out)

            prmtop = app.AmberPrmtopFile(os.path.join(tmpDir, 'mol.prmtop'))

        return prmtop

    def minimize(self, coords, restrained_dihedrals=None):
        from simtk import unit
        from simtk.openmm import app
        import simtk.openmm as mm
        import nlopt

        if coords.ndim == 3:
            coords = coords[:, :, 0]

        integrator = mm.LangevinIntegrator(0, 0, 0)
        sim = app.Simulation(self.structure.topology, self.system, integrator, self.platform, self.platprop)

        otheridx = np.arange(coords.shape[0])
        if restrained_dihedrals is not None:
            dihidx = np.concatenate(restrained_dihedrals).flatten()
            otheridx = np.setdiff1d(np.arange(coords.shape[0]), dihidx)

        if len(otheridx) == 0:  # If the molecule only consists of the dihedral atoms return
            return coords

        def goalFunc(x, grad):
            currcoords = coords.copy()
            currcoords[otheridx] = x.reshape((len(otheridx), 3))

            sim.context.setPositions(currcoords * unit.angstrom)
            state = sim.context.getState(getEnergy=True, getForces=True)
            energy = state.getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole)
            forces = state.getForces(asNumpy=True).value_in_unit(unit.kilocalories_per_mole / unit.angstrom)
            grad[:] = -forces[otheridx].reshape(-1)
            return energy

        opt = nlopt.opt(nlopt.LD_LBFGS, len(otheridx)*3)
        opt.set_min_objective(goalFunc)
        opt.set_ftol_abs(1E-4)
        opt.set_xtol_abs(1E-6)
        x = opt.optimize(coords[otheridx].reshape(-1))
        endcoords = coords.copy()
        endcoords[otheridx] = x.reshape((len(otheridx), 3))

        return endcoords


class CustomEnergyBasedMinimizer(Minimizer):
    def __init__(self, mol, calculator):
        super().__init__()
        self.opt = nlopt.opt(nlopt.LN_COBYLA, mol.coords.size)

        def objective(x, _):
            return float(calculator.calculate(x.reshape((-1, 3, 1)), mol.element, units='kcalmol')[0])

        self.opt.set_min_objective(objective)

    def minimize(self, coords, restrained_dihedrals):
        if restrained_dihedrals is not None:
            for dihedral in restrained_dihedrals:
                indices = dihedral.copy()
                ref_angle = dihedralAngle(coords[indices, :, 0])

                def constraint(x, _):
                    coords = x.reshape((-1, 3))
                    angle = dihedralAngle(coords[indices])
                    return np.sin(.5 * (angle - ref_angle))

                self.opt.add_equality_constraint(constraint)

        self.opt.set_xtol_abs(1e-3)  # Similar to Psi4 default
        self.opt.set_maxeval(1000 * self.opt.get_dimension())
        self.opt.set_initial_step(1e-3)
        return self.opt.optimize(coords.ravel()).reshape((-1, 3, 1))


class CustomQM(QMBase):
    """
    Imitation of QM calculations with custom class

    >>> import os
    >>> import numpy as np
    >>> from tempfile import TemporaryDirectory
    >>> from htmd.home import home
    >>> from htmd.numbautil import dihedralAngle
    >>> from htmd.molecule.molecule import Molecule
    >>> from htmd.qm.custom import CustomQM
    >>> from acemdai.calculator import AAICalculator

    Create a molecule
    >>> molFile = os.path.join(home('test-qm'), 'H2O2-90.mol2')
    >>> mol = Molecule(molFile)

    Create the AcemdAI calculator
    >>> networkfile = './mynet.pkl'
    >>> aceai = AAICalculator(networkfile=networkfile, maxatoms=26, maxneighs=50)

    Run a single-point energy and ESP calculation
    >>> with TemporaryDirectory() as tmpDir:
    ...     qm = CustomQM()
    ...     qm.calculator = aceai
    ...     qm.molecule = mol
    ...     qm.esp_points = np.array([[1., 1., 1.]])
    ...     qm.directory = tmpDir
    ...     result = qm.run()[0]
    4 elements | 10 element pairs | 384 features
    CUDA: Allocating features array (1, 4, 384)
    CUDA: Allocating gradient array (1, 4, 4, 384, 3)
    >>> qm # doctest: +ELLIPSIS
    <htmd.qm.custom.CustomCalculator object at ...>
    >>> result # doctest: +ELLIPSIS
    <htmd.qm.base.QMResult object at ...
    >>> result.errored
    False
    >>> result.energy # doctest: +ELLIPSIS
    -94970.499...
    >>> np.rad2deg(dihedralAngle(result.coords[[2, 0, 1, 3], :, 0])) # doctest: +ELLIPSIS
    89.99...

    Run a minimization
    >>> with TemporaryDirectory() as tmpDir:
    ...     qm = CustomQM()
    ...     qm.calculator = aceai
    ...     qm.molecule = mol
    ...     qm.optimize = True
    ...     qm.directory = tmpDir
    ...     result = qm.run()[0]
    4 elements | 10 element pairs | 384 features
    CUDA: Allocating features array (1, 4, 384)
    CUDA: Allocating gradient array (1, 4, 4, 384, 3)
    >>> result.energy # doctest: +ELLIPSIS
    -95173.433...
    >>> np.rad2deg(dihedralAngle(result.coords[[2, 0, 1, 3], :, 0])) # doctest: +ELLIPSIS
    125.993...

    Run a constrained minimization
    >>> with TemporaryDirectory() as tmpDir:
    ...     qm = CustomQM()
    ...     qm.calculator = aceai
    ...     qm.molecule = mol
    ...     qm.optimize = True
    ...     qm.restrained_dihedrals = np.array([[2, 0, 1, 3]])
    ...     qm.directory = tmpDir
    ...     result = qm.run()[0]
    4 elements | 10 element pairs | 384 features
    CUDA: Allocating features array (1, 4, 384)
    CUDA: Allocating gradient array (1, 4, 4, 384, 3)
    >>> result.energy # doctest: +ELLIPSIS
    -95170.800...
    >>> np.rad2deg(dihedralAngle(result.coords[[2, 0, 1, 3], :, 0])) # doctest: +ELLIPSIS
    89.99...
    """

    def __init__(self, verbose=True):
        super().__init__()
        self._verbose = verbose
        self._arg('calculator', ':class: `Calculator`', 'Calculator object', default=None, validator=None, required=True)
        self._arg('minimizer', ':class: `Minimizer`', 'Minimizer object', default=None, validator=None)

    # Fake implementations of the abstract methods
    def _command(self): pass
    def _writeInput(self, directory, iframe): pass
    def _readOutput(self, directory): pass
    def setup(self): pass
    def submit(self): pass

    def _completed(self, directory):
        return os.path.exists(os.path.join(directory, 'data.pkl'))

    def retrieve(self):

        results = []
        for iframe in range(self.molecule.numFrames):
            self.molecule.frame = iframe

            directory = os.path.join(self.directory, '%05d' % iframe)
            os.makedirs(directory, exist_ok=True)
            pickleFile = os.path.join(directory, 'data.pkl')
            molFile = os.path.join(directory, 'mol.mol2')

            if self._completed(directory):
                with open(pickleFile, 'rb') as fd:
                    result = pickle.load(fd)
                logger.info('Loading data from %s' % pickleFile)

            else:
                start = time.clock()

                result = QMResult()
                result.errored = False
                result.coords = self.molecule.coords[:, :, iframe:iframe + 1].copy()

                if self.optimize:
                    if self.minimizer is None:
                        self.minimizer = CustomEnergyBasedMinimizer(self.molecule, self.calculator)
                    result.coords = self.minimizer.minimize(result.coords, self.restrained_dihedrals).reshape((-1, 3, 1))
                    mol = self.molecule.copy()
                    mol.frame = 0
                    mol.coords = result.coords
                    mol.write(molFile)

                result.energy = float(self.calculator.calculate(result.coords, self.molecule.element, units='kcalmol')[0])
                result.dipole = [0, 0, 0]

                #if self.optimize:
                #    assert opt.last_optimum_value() == result.energy # A self-consistency test

                finish = time.clock()
                result.calculator_time = finish - start
                if self._verbose:
                    logger.info('Custom calculator calculation time: %f s' % result.calculator_time)

                with open(pickleFile, 'wb') as fd:
                    pickle.dump(result, fd)

            results.append(result)

        return results


if __name__ == '__main__':

    import sys
    # TODO: Currently doctest is not working correctly, and qmml module is not made available either
    # import doctest
    #
    # if doctest.testmod().failed:
    #     sys.exit(1)
