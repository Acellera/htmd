# (c) 2015-2018 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import os
import numpy as np
from scipy import constants as const

from htmd.qm.base import QMBase, QMResult


class Gaussian(QMBase):
    """
    Class to set up and run QM calculations with Gaussian

   See htmd.qm.Psi4 documentation for more details.

    Examples
    --------

    Create a Gaussian object
    >>> from htmd.qm import Gaussian
    >>> qm = Gaussian()
    >>> qm # doctest: +ELLIPSIS
    <htmd.qm.gaussian.Gaussian object at 0x...>
    """

    @property
    def _command(self):
        return 'g09 < g09.in > g09.out 2>&1'

    def _completed(self, directory):
        return os.path.exists(os.path.join(directory, 'g09.out'))

    def _writeInput(self, directory, iframe):

        with open(os.path.join(directory, 'g09.in'), 'w') as f:

            f.write('%nprocshared = {}\n'.format(self.queue.ncpu))
            memory = int(np.floor(self.queue.memory/1000))
            f.write('%mem = {} GB\n\n'.format(memory))

            if self.theory == 'HF':
                dispersion = ''
            elif self.theory == 'B3LYP':
                dispersion = ' EmpiricalDispersion=GD3'
            else:
                raise NotImplementedError

            line = '# %s/%s NoSymm SCF=tight' % (self.theory, self.basis)
            line += dispersion

            if self.solvent == 'vacuum':
                pass
            elif self.solvent == 'PCM':
                line += ' SCRF=PCM'
            else:
                raise NotImplementedError

            if self.optimize:
                line += ' Opt=ModRedundant'

            if self.esp_points is not None:
                line +=  ' Prop=(read,field)'

            line += '\n\n'
            f.write(line)

            f.write('Mol\n%d %d\n' % (self._charge, self.multiplicity))
            elements = self._molecule.element
            coords = self._molecule.coords[:, :, iframe]
            for element, coord in zip(elements, coords):
                f.write('    %-2s %10f %10f %10f\n' % (element, coord[0], coord[1], coord[2]))
            f.write('\n')

            if self._restrained_dihedrals is not None:
                for dihedral in self._restrained_dihedrals:
                    f.write('%d %d %d %d F\n' % tuple(dihedral))
                f.write('\n')

            if self.esp_points is not None:
                f.write('@grid.dat /N\n')

    def _readOutput(self, directory):

        result = QMResult()

        with open(os.path.join(directory, 'g09.out')) as f:
            fl = f.readlines()

        try:
            data = {}

            for i in range(len(fl)):
                if "Dipole moment (field-independent basis, Debye):" in fl[i]:
                    s = fl[i + 1].split()
                    data['dipole'] = [float(s[1]), float(s[3]), float(s[5]), float(s[7])]
                if "Traceless Quadrupole moment (field-independent basis, Debye-Ang):" in fl[i]:
                    s1 = fl[i + 1].split()
                    s2 = fl[i + 2].split()
                    data['quadrupole'] = [float(s1[1]), float(s1[3]), float(s1[5]), float(s2[1]), float(s2[3]),
                                          float(s2[5])]
                if "Mulliken atomic charges:" in fl[i] or "Mulliken charges:" in fl[i]:
                    data['mulliken'] = []
                    for j in range(self._natoms):
                        data['mulliken'].append(float(fl[i + 2 + j].split()[2]))

            for l in fl:
                if "SCF Done:  E(RHF) = " in l:
                    ff = l.split()
                    data['energy'] = float(ff[4])
                if "SCF Done:  E(RB3LYP) = " in l:
                    ff = l.split()
                    data['energy'] = float(ff[4])
            i = 0
            while i < len(fl):
                if "Number     Number       Type             X           Y           Z" in fl[i]:
                    i += 2
                    data['coords'] = np.zeros((self._natoms, 3))
                    for j in range(self._natoms):
                        ff = fl[j + i].split()
                        data['coords'][j, 0] = float(ff[3])
                        data['coords'][j, 1] = float(ff[4])
                        data['coords'][j, 2] = float(ff[5])
                else:
                    i += 1

            i = 0
            while i < len(fl):
                if "               Potential          X             Y             Z" in fl[i]:
                    i = i + 2 + self._natoms
                    data['gridesp'] = []
                    while not ("----" in fl[i]):
                        ff = fl[i].split()
                        data['gridesp'].append(float(ff[1]))
                        i += 1
                    data['gridesp'] = np.asarray(data['gridesp'], dtype=np.float)

                else:
                    i += 1

        except:
            data = None

        if data is None:
            result.errored = True
        elif not ('energy' in data) or data['energy'] == 0.:
            result.errored = True
        else:
            result.energy = data['energy']
            result.energy *= const.physical_constants['Hartree energy'][0]/(const.kilo*const.calorie/const.Avogadro) # Hartree to kcal/mol
            result.coords = np.atleast_3d(data['coords'])
            if 'gridesp' in data:
                result.esp_points = np.squeeze(np.asarray(self.esp_points))
                result.esp_values = data['gridesp']
                result.esp_values *= const.angstrom/const.physical_constants['Bohr radius'][0] # 1/Bohr to 1/Angstrom
            if 'dipole' in data:
                result.dipole = data['dipole']
            if 'quadrupole' in data:
                result.quadrupole = data['quadrupole']
            if 'mulliken' in data:
                result.mulliken = data['mulliken']

        return result


if __name__ == '__main__':

    import sys
    import doctest

    if doctest.testmod().failed:
        sys.exit(1)