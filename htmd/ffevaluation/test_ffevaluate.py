# (c) 2015-2018 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
from htmdmol.molecule import Molecule
from htmd.ffevaluation.ffevaluate import FFEvaluate
import parmed
from glob import glob
import numpy as np
import os


def disableDispersionCorrection(system):
    # According to openMM:
    # The long range dispersion correction is primarily useful when running simulations at constant pressure, since it
    # produces a more accurate variation in system energy with respect to volume.
    # So I will disable it to avoid implementing it for now in ffevaluate
    from simtk.openmm import NonbondedForce
    for f in system.getForces():
        if isinstance(f, NonbondedForce):
            f.setUseDispersionCorrection(False)

def openmm_energy(prm, structure, coords, box=None, cutoff=None):
    import parmed
    from simtk import unit
    from simtk import openmm
    from simtk.openmm import app
    from parmed.amber import AmberParm

    if box is not None and not np.all(box == 0):
        if cutoff is None:
            raise RuntimeError('You need to provide a cutoff when passing a box')
        a = unit.Quantity((box[0] * unit.angstrom, 0 * unit.angstrom, 0 * unit.angstrom))
        b = unit.Quantity((0 * unit.angstrom, box[1] * unit.angstrom, 0 * unit.angstrom))
        c = unit.Quantity((0 * unit.angstrom, 0 * unit.angstrom, box[2] * unit.angstrom))
        structure.box_vectors = (a, b, c)
        if isinstance(structure, AmberParm):
            system = structure.createSystem(nonbondedMethod=app.CutoffPeriodic, nonbondedCutoff=cutoff*unit.angstrom)
        else:
            system = structure.createSystem(prm, nonbondedMethod=app.CutoffPeriodic, nonbondedCutoff=cutoff*unit.angstrom)
        system.setDefaultPeriodicBoxVectors(a, b, c)
    else:
        if isinstance(structure, AmberParm):
            system = structure.createSystem()
        else:
            system = structure.createSystem(prm)

    disableDispersionCorrection(system)
    integrator = openmm.LangevinIntegrator(300 * unit.kelvin, 1 / unit.picoseconds, 2 * unit.femtoseconds)
    platform = openmm.Platform.getPlatformByName('CPU')
    context = openmm.Context(system, integrator, platform)

    # Run OpenMM with given coordinates
    context.setPositions(coords * unit.angstrom)
    energies = parmed.openmm.energy_decomposition(structure, context)
    state = context.getState(getForces=True)
    forces = state.getForces(asNumpy=True).value_in_unit(unit.kilocalories_per_mole/unit.angstrom)

    return energies, forces


def keepForces(prm, psf, mol, forces=('angle', 'bond', 'dihedral', 'lennardjones', 'electrostatic'), verbose=False):
    from collections import OrderedDict
    if 'angle' not in forces:
        if verbose: print('Disabling angle forces')
        for type in prm.angle_types:
            prm.angle_types[type].k = 0
    if 'bond' not in forces:
        if verbose: print('Disabling bonded forces')
        for type in prm.bond_types:
            prm.bond_types[type].k = 0
    if 'dihedral' not in forces:
        if verbose: print('Disabling dihedral forces')
        for type in prm.dihedral_types:
            for dih in prm.dihedral_types[type]:
                dih.phi_k = 0
        if verbose: print('Disabling improper forces')
        for type in prm.improper_types:
            prm.improper_types[type].psi_k = 0
        for type in prm.improper_periodic_types:
            prm.improper_periodic_types[type].phi_k = 0
    if 'lennardjones' not in forces:
        if verbose: print('Disabling LJ forces')
        for type in prm.atom_types:
            prm.atom_types[type].epsilon = prm.atom_types[type].epsilon_14 = 0
            prm.atom_types[type].sigma = prm.atom_types[type].sigma_14 = 0
            prm.atom_types[type].nbfix = {}
        prm.nbfix_types = OrderedDict()
    if 'electrostatic' not in forces:
        if verbose: print('Disabling electrostatic forces')
        for res in prm.residues:
            for atom in prm.residues[res]:
                atom.charge = 0
        for a in psf.atoms:
            a.charge = 0
        mol.charge[:] = 0


def keepForcesAmber(struct, mol, forces=('angle', 'bond', 'dihedral', 'lennardjones', 'electrostatic'), verbose=False):
    if 'angle' not in forces:
        if verbose: print('Disabling angle forces')
        for i in range(len(struct.angle_types)):
            struct.angle_types[i].k = 0
    if 'bond' not in forces:
        if verbose: print('Disabling bonded forces')
        for i in range(len(struct.bond_types)):
            struct.bond_types[i].k = 0
    if 'dihedral' not in forces:
        if verbose: print('Disabling dihedral forces')
        for i in range(len(struct.dihedral_types)):
            struct.dihedral_types[i].phi_k = 0
        if verbose: print('Disabling improper forces')
        for i in range(len(struct.improper_types)):
            struct.improper_types[i].psi_k = 0
    if 'lennardjones' not in forces:
        if verbose: print('Disabling LJ forces')
        for i in range(len(struct.atoms)):
            struct.atoms[i].epsilon = struct.atoms[i].epsilon_14 = 0
            # struct.atoms[i].nbfix = {}
        # prm.nbfix_types = OrderedDict()
    if 'electrostatic' not in forces:
        if verbose: print('Disabling electrostatic forces')
        for res in range(len(struct.residues)):
            for atom in struct.residues[res]:
                atom.charge = 0
        for a in struct.atoms:
            a.charge = 0
        mol.charge[:] = 0


def compareEnergies(myenergies, omm_energies, verbose=False, abstol=1e-4):
    if 'angle' in omm_energies:
        d = myenergies['angle'] - omm_energies['angle']
        if abs(d) > abstol:
            raise RuntimeError('Too high difference in angle energies:', d)
        if verbose: print('Angle diff:', d)
    if 'bond' in omm_energies:
        d = myenergies['bond'] - omm_energies['bond']
        if abs(d) > abstol:
            raise RuntimeError('Too high difference in bond energies:', d)
        if verbose: print('Bond diff:', d)
    if 'dihedral' in omm_energies:
        d = (myenergies['dihedral'] + myenergies['improper']) - omm_energies['dihedral']
        if abs(d) > abstol:
            raise RuntimeError('Too high difference in dihedral energies:', d)
        if verbose: print('Dihedral diff:', d)
    if 'nonbonded' in omm_energies:
        d = (myenergies['vdw'] + myenergies['elec']) - omm_energies['nonbonded']
        if abs(d) > abstol:
            raise RuntimeError('Too high difference in nonbonded energies:', d)
        if verbose: print('Nonbonded diff:', d)
    if 'improper' in omm_energies:
        d = myenergies['improper'] - omm_energies['improper']
        if abs(d) > abstol:
            raise RuntimeError('Too high difference in improper energies:', d)
        if verbose: print('Improper diff:', d)
    d = myenergies['total'] - omm_energies['total']
    if abs(d) > abstol:
        raise RuntimeError('Too high difference in total energy:', d)
    if verbose: print('Total diff:', d)
    return d


def compareForces(forces1, forces2):
    return np.max(np.abs(forces1 - forces2).flatten())


def fixParameters(parameterfile):
    from htmd.util import tempname
    tmpfile = tempname(suffix='.prm')
    with open(parameterfile, 'r') as f:
        lines = f.readlines()

    with open(tmpfile, 'w') as f:
        for l in lines:
            l = l.replace('!MASS', 'MASS')
            l = l.replace('!ATOMS', 'ATOMS')
            f.write(l)
    return tmpfile


if __name__ == '__main__':
    from natsort import natsorted
    from htmd.home import home
    from htmdmol.molecule import Molecule
    from glob import glob
    import parmed
    import os
    import logging

    logging.getLogger('parmed.structure').setLevel('ERROR')

    for d in glob(os.path.join(home(dataDir='test-ffevaluate'), '*', '')):
        print('\nRunning test:', d)
        if os.path.basename(os.path.abspath(d)) == 'thrombin-ligand-amber':
            abstol = 1e-1
        elif os.path.basename(os.path.abspath(d)) == 'waterbox':
            abstol = 1e-3
        else:
            abstol = 1e-4

        prmtopFile = glob(os.path.join(d, '*.prmtop'))
        psfFile = glob(os.path.join(d, '*.psf'))
        pdbFile = glob(os.path.join(d, '*.pdb'))
        xtcFile = glob(os.path.join(d, '*.xtc'))
        if len(glob(os.path.join(d, '*.prm'))):
            prmFiles = [fixParameters(glob(os.path.join(d, '*.prm'))[0]), ]
        rtfFile = glob(os.path.join(d, '*.rtf'))
        if len(rtfFile):
            prmFiles.append(rtfFile[0])
        else:
            rtfFile = None

        if len(psfFile):
            mol = Molecule(psfFile[0])
        elif len(prmtopFile):
            mol = Molecule(prmtopFile[0])
        if len(xtcFile):
            mol.read(natsorted(xtcFile))
        elif len(pdbFile):
            mol.read(pdbFile[0])
        else:
            raise RuntimeError('No PDB or XTC')
        coords = mol.coords
        coords = coords[:, :, 0].squeeze()
        rfa = False
        cutoff = 0
        if not np.all(mol.box == 0):
            cutoff = np.min(mol.box) / 2 - 0.01
            rfa = True

        chargebackup = mol.charge.copy()
        for force in ('angle', 'bond', 'dihedral', 'lennardjones', 'electrostatic'):
            mol.charge = chargebackup.copy()
            if len(psfFile):
                struct = parmed.charmm.CharmmPsfFile(psfFile[0])
                prm = parmed.charmm.CharmmParameterSet(*prmFiles)
                keepForces(prm, struct, mol, forces=force)
            elif len(prmtopFile):
                struct = parmed.load_file(prmtopFile[0])
                prm = parmed.amber.AmberParameterSet().from_structure(struct)
                keepForces(prm, struct, mol, forces=force)
                keepForcesAmber(struct, mol, forces=force)

            energies, forces, atmnrg = FFEvaluate(mol, prm, cutoff=cutoff, rfa=rfa).calculate(mol.coords, mol.box)
            energies = FFEvaluate.formatEnergies(energies[:, 0])
            forces = forces[:, :, 0].squeeze()
            omm_energies, omm_forces = openmm_energy(prm, struct, coords, box=mol.box, cutoff=cutoff)
            ediff = compareEnergies(energies, omm_energies, abstol=abstol)
            print('  ', force, 'Energy diff:', ediff, 'Force diff:', compareForces(forces, omm_forces))

        if len(psfFile):
            struct = parmed.charmm.CharmmPsfFile(psfFile[0])
            prm = parmed.charmm.CharmmParameterSet(*prmFiles)
            keepForces(prm, struct, mol)
        elif len(prmtopFile):
            struct = parmed.load_file(prmtopFile[0])
            prm = parmed.amber.AmberParameterSet().from_structure(struct)
            keepForces(prm, struct, mol)
            keepForcesAmber(struct, mol)
        energies, forces, atmnrg = FFEvaluate(mol, prm, cutoff=cutoff, rfa=rfa).calculate(mol.coords, mol.box)
        energies = FFEvaluate.formatEnergies(energies[:, 0])
        forces = forces[:, :, 0].squeeze()
        omm_energies, omm_forces = openmm_energy(prm, struct, coords, box=mol.box, cutoff=cutoff)
        ediff = compareEnergies(energies, omm_energies, abstol=abstol)
        print('All forces. Total energy:', energies['total'], 'Energy diff:', ediff, 'Force diff:', compareForces(forces, omm_forces))




