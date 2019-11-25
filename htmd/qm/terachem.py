# (c) 2015-2019 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import os
import numpy as np
from scipy import constants as const

from moleculekit.molecule import Molecule
from htmd.qm.base import QMBase, QMResult


class TeraChem(QMBase):

    @property
    def _command(self):
        if 'TeraChem' not in os.environ:
            raise RuntimeError('Environment variable TeraChem is undefined')
        if not os.path.exists(os.path.join(os.environ['TeraChem'], 'SetTCVars.sh')):
            raise RuntimeError('$TeraChem/SetTCVars.sh does not exist')
        return 'source $TeraChem/SetTCVars.sh && ' \
               'terachem terachem.in &> terachem.out'

    def _completed(self, directory):

        outFile = os.path.join(directory, 'terachem.out')

        if not os.path.exists(outFile):
            return False

        with open(outFile) as fd:
            for line in fd.readlines():
                if line.startswith(' Job terminated:'):
                    return True

        return False

    def _writeInput(self, directory, frame):

        xyzFile = os.path.join(directory, 'terachem.xyz')
        self.molecule.frame = frame
        self.molecule.write(xyzFile)

        with open(os.path.join(directory, 'terachem.in'), 'w') as f:

            f.write(f'method {("R" if self.multiplicity == 1 else "U") + self.theory}\n')
            f.write(f'basis {self.basis}\n')
            if self.correction != 'none':
                f.write(f'dispersion {"D2" if self.correction == "D" else self.correction}\n')
            if self.solvent == 'PCM':
                f.write('pcm cosmo\n')

            f.write(f'charge {self.charge}\n')
            f.write(f'spinmult {self.multiplicity}\n')
            f.write('coordinates terachem.xyz\n')

            if self.esp_points is not None:
                raise NotImplemented('ESP is not available')

            if self.optimize:
                f.write('new_minimizer yes\n')
                if self._restrained_dihedrals is not None:
                    f.write('$constraint_freeze\n')
                    for dihedral in self._restrained_dihedrals:
                        f.write(f'  dihedral {"_".join(map(str, dihedral))}\n')
                    f.write('$end\n')
                f.write('run minimize\n')
            else:
                f.write('run energy\n')

            f.write('precision double\n')

            f.write(f'end\n')

    def _readOutput(self, directory):

        result = QMResult()
        result.completed = True

        with open(os.path.join(directory, 'terachem.out')) as fd:
            for line in fd.readlines():
                if line.startswith('FINAL ENERGY:'):
                    result.energy = float(line.split()[2])
                    result.energy *= const.physical_constants['Hartree energy'][0] /\
                                     (const.kilo * const.calorie / const.Avogadro)  # Hartree --> kcal/mol
                if line.startswith('DIPOLE MOMENT:'):
                    tokens = line.split()
                    result.dipole = [float(tokens[2][1:-1]),
                                     float(tokens[3][ :-1]),
                                     float(tokens[4][ :-1]),
                                     float(tokens[7][ :-1])]

        result.mulliken = np.loadtxt(os.path.join(directory, 'scr', 'charge_mull.xls'), usecols=(2,))

        if result.energy is None or result.dipole is None:
            result.errored = True

        geomFile = os.path.join(directory, 'scr', 'optim.xyz' if self.optimize else 'xyz.xyz')
        if os.path.exists(geomFile):
            result.coords = Molecule(geomFile).coords[:,:,-1][:,:,None]
        else:
            result.errored = True

        if self.esp_points is not None:
            raise NotImplemented('ESP is not available')

        return result


if __name__ == '__main__':

    import sys
    import doctest

    sys.exit(doctest.testmod().failed)
