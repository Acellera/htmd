from htmd.molecule.molecule import Molecule
from htmd.ffevaluation.ffevaluate import ffevaluate
import parmed
from glob import glob
import numpy as np
import os


def openmm_energy(prm, structure, coords, box=None):
    import parmed
    from simtk import unit
    from simtk import openmm
    from simtk.openmm import app

    if box is not None:
        a = unit.Quantity((box[0] * unit.angstrom, 0 * unit.angstrom, 0 * unit.angstrom))
        b = unit.Quantity((0 * unit.angstrom, box[1] * unit.angstrom, 0 * unit.angstrom))
        c = unit.Quantity((0 * unit.angstrom, 0 * unit.angstrom, box[2] * unit.angstrom))
        structure.box_vectors = (a, b, c)
        system = structure.createSystem(prm, nonbondedMethod=app.CutoffPeriodic)
    else:
        system = structure.createSystem(prm)
    integrator = openmm.LangevinIntegrator(300 * unit.kelvin, 1 / unit.picoseconds, 2 * unit.femtoseconds)
    platform = openmm.Platform.getPlatformByName('CPU')
    context = openmm.Context(system, integrator, platform)

    # Run OpenMM with given coordinates
    context.setPositions(coords * unit.angstrom)
    energies = parmed.openmm.energy_decomposition(psf, context)
    state = context.getState(getForces=True)
    forces = state.getForces(asNumpy=True).value_in_unit(unit.kilocalories_per_mole/unit.angstrom)

    return energies, forces


def drawForce(start, vec):
    assert start.ndim == 1 and vec.ndim == 1
    from htmd.vmdviewer import getCurrentViewer
    vmd = getCurrentViewer()
    vmd.send("""
    proc vmd_draw_arrow {start end} {
        # an arrow is made of a cylinder and a cone
        draw color green
        set middle [vecadd $start [vecscale 0.9 [vecsub $end $start]]]
        graphics top cylinder $start $middle radius 0.15
        graphics top cone $middle $end radius 0.25
    }
    """)
    vmd.send('vmd_draw_arrow {{ {} }} {{ {} }}'.format(' '.join(map(str, start)), ' '.join(map(str, start + vec))))


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
    if 'lennardjones' not in forces:
        if verbose: print('Disabling LJ forces')
        for type in prm.atom_types:
            prm.atom_types[type].epsilon = prm.atom_types[type].epsilon_14 = 0
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
        d = myenergies['dihedral'] - omm_energies['dihedral']
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


def viewForces(mol, forces, omm_forces):
    mol.view()
    for cc, ff in zip(mol.coords[:, :, 0], forces):
        drawForce(cc, ff)

    mol.view()
    for cc, ff in zip(mol.coords[:, :, 0], omm_forces):
        drawForce(cc, ff)


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
    from htmd.ffevaluation.ffevaluate import _formatEnergies
    from htmd.home import home
    from htmd.molecule.molecule import Molecule
    from glob import glob
    import os

    for d in glob(os.path.join(home(dataDir='test-ffevaluate'), '*', '')):
        print('\nRunning test:', d)
        if os.path.basename(os.path.abspath(d)) == 'waterbox':
            abstol = 1e-3
        else:
            abstol = 1e-4

        psfFile = glob(os.path.join(d, '*.psf'))[0]
        pdbFile = glob(os.path.join(d, '*.pdb'))[0]
        xtcFile = glob(os.path.join(d, '*.xtc'))
        prmFiles = [fixParameters(glob(os.path.join(d, '*.prm'))[0]), ]
        rtfFile = glob(os.path.join(d, '*.rtf'))
        if len(rtfFile):
            prmFiles.append(rtfFile[0])
        else:
            rtfFile = None

        mol = Molecule(psfFile)
        if len(xtcFile):
            mol.read(natsorted(xtcFile))
        else:
            mol.read(pdbFile)
        coords = mol.coords
        coords = coords[:, :, 0].squeeze()
        mol.box[:] = 0

        chargebackup = mol.charge.copy()
        for force in ('angle', 'bond', 'dihedral', 'lennardjones', 'electrostatic'):
            mol.charge = chargebackup.copy()
            psf = parmed.charmm.CharmmPsfFile(psfFile)
            prm = parmed.charmm.CharmmParameterSet(*prmFiles)
            keepForces(prm, psf, mol, forces=force)
            energies, forces, atmnrg = ffevaluate(mol, prm)
            energies = _formatEnergies(energies[:, 0])
            forces = forces[:, :, 0].squeeze()
            omm_energies, omm_forces = openmm_energy(prm, psf, coords)
            ediff = compareEnergies(energies, omm_energies, abstol=abstol)
            print('  ', force, 'Energy diff:', ediff, 'Force diff:', compareForces(forces, omm_forces))

        psf = parmed.charmm.CharmmPsfFile(psfFile)
        prm = parmed.charmm.CharmmParameterSet(*prmFiles)
        keepForces(prm, psf, mol)
        energies, forces, atmnrg = ffevaluate(mol, prm)
        energies = _formatEnergies(energies[:, 0])
        forces = forces[:, :, 0].squeeze()
        omm_energies, omm_forces = openmm_energy(prm, psf, coords)
        ediff = compareEnergies(energies, omm_energies, abstol=abstol)
        print('All forces. Total energy:', energies['total'], 'Energy diff:', ediff, 'Force diff:', compareForces(forces, omm_forces))

        # viewForces(mol, forces, omm_forces)




