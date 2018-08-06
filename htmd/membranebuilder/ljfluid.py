# (c) 2015-2018 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
from simtk import openmm
from simtk import unit
from simtk.openmm import app
import numpy as np


def _halton_sequence(p, n):
    """
    Halton deterministic sequence on [0,1].
    Parameters
    ----------
    p : int
       Prime number for sequence.
    n : int
       Sequence length to generate.
    Returns
    -------
    u : numpy.array of double
       Sequence on [0,1].
    Notes
    -----
    Code source: http://blue.math.buffalo.edu/sauer2py/
    More info: http://en.wikipedia.org/wiki/Halton_sequence
    Examples
    --------
    Generate some sequences with different prime number bases.
    >>> x = _halton_sequence(2,100)
    >>> y = _halton_sequence(3,100)
    >>> z = _halton_sequence(5,100)
    """
    eps = np.finfo(np.double).eps
    # largest number of digits (adding one for halton_sequence(2,64) corner case)
    b = np.zeros(int(np.ceil(np.log(n) / np.log(p))) + 1)
    u = np.empty(n)
    for j in range(n):
        i = 0
        b[0] += 1                       # add one to current integer
        while b[i] > p - 1 + eps:           # this loop does carrying in base p
            b[i] = 0
            i = i + 1
            b[i] += 1
        u[j] = 0
        for k in range(len(b)):         # add up reversed digits
            u[j] += b[k] * p**-(k + 1)
    return u


def _subrandom_particle_positions(nparticles, box_vectors, ndim, method='halton'):
    """Generate a deterministic list of subrandom particle positions."""
    # Create positions array.
    positions = np.zeros([nparticles, 3], np.float32)

    # Fill in each dimension.
    primes = [2, 3, 5]  # prime bases for Halton sequence
    for dim in range(ndim):
        x = _halton_sequence(primes[dim], nparticles)
        l = box_vectors[dim][dim]
        positions[:, dim] = x - 0.5

    np.random.shuffle(positions)
    return unit.Quantity(positions * l / l.unit, l.unit)


def distributeLipids(boxsize,
                     resnames,
                     sigmas,
                     cutoff,
                     mass=39.9 * unit.amu,  # argon
                     epsilon=0.238 * unit.kilocalories_per_mole,  # argon,
                     switch_width=3.4 * unit.angstrom,  # argon
                     ):
        nparticles = len(resnames)
                
        # Determine Lennard-Jones cutoff.
        cutoff = cutoff * unit.angstrom

        cutoff_type = openmm.NonbondedForce.CutoffPeriodic

        # Create an empty system object.
        system = openmm.System()

        # Periodic box vectors.
        a = unit.Quantity((boxsize[0] * unit.angstrom, 0 * unit.angstrom, 0 * unit.angstrom))
        b = unit.Quantity((0 * unit.angstrom, boxsize[1] * unit.angstrom, 0 * unit.angstrom))
        c = unit.Quantity((0 * unit.angstrom, 0 * unit.angstrom, boxsize[2] * unit.angstrom))
        system.setDefaultPeriodicBoxVectors(a, b, c)

        # Set up periodic nonbonded interactions with a cutoff.
        nb = openmm.NonbondedForce()
        nb.setNonbondedMethod(cutoff_type)
        nb.setCutoffDistance(cutoff)
        nb.setUseDispersionCorrection(True)

        nb.setUseSwitchingFunction(False)
        if (switch_width is not None):
            nb.setUseSwitchingFunction(True)
            nb.setSwitchingDistance(cutoff - switch_width)

        for s in sigmas:
            system.addParticle(mass)
            nb.addParticle(0.0 * unit.elementary_charge, s * unit.angstrom, epsilon)

        positions = _subrandom_particle_positions(nparticles, system.getDefaultPeriodicBoxVectors(), 2)

        # Add the nonbonded force.
        system.addForce(nb)

        # Add a restraining potential to keep atoms in z=0
        energy_expression = 'k * (z^2)'
        force = openmm.CustomExternalForce(energy_expression)
        force.addGlobalParameter('k', 10)
        for particle_index in range(nparticles):
            force.addParticle(particle_index, [])
        system.addForce(force)

        # Create topology.
        topology = app.Topology()
        chain = topology.addChain()
        elems = ['Ar', 'Cl', 'Na']
        _, idx = np.unique(resnames, return_inverse=True)
        for i in idx:
            element = app.Element.getBySymbol(elems[i])
            residue = topology.addResidue(elems[i], chain)
            topology.addAtom(elems[i], element, residue)

        topology.setUnitCellDimensions(unit.Quantity(boxsize, unit.angstrom)) 
            
        # Simulate it
        from simtk.openmm import VerletIntegrator
        from simtk.openmm.app import Simulation
        from simtk.unit import picoseconds, picosecond, angstrom
        nsteps = 10000

        integrator = VerletIntegrator(0.002 * picoseconds)
        simulation = Simulation(topology, system, integrator)
        simulation.context.setPositions(positions)
        simulation.minimizeEnergy()
        # simulation.reporters.append(DCDReporter('output.dcd', 1))
        # simulation.reporters.append(StateDataReporter(stdout, 1000, potentialEnergy=True, totalEnergy=True, step=True, separator='   '))
        simulation.step(nsteps)

        state = simulation.context.getState(getPositions=True, enforcePeriodicBox=True)
        allfinalpos = state.getPositions(asNumpy=True).value_in_unit(angstrom)

        return allfinalpos


if __name__ == '__main__':
    pass
