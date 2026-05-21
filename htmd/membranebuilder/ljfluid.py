# (c) 2015-2022 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
try:
    import openmm
except ImportError:
    raise ImportError(
        "openmm is not installed. Please install it using `pip install openmm`."
    )
from openmm import unit
from openmm import app
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
        b[0] += 1  # add one to current integer
        while b[i] > p - 1 + eps:  # this loop does carrying in base p
            b[i] = 0
            i = i + 1
            b[i] += 1
        u[j] = 0
        for k in range(len(b)):  # add up reversed digits
            u[j] += b[k] * p ** -(k + 1)
    return u


def _generate_halton_positions(n, box_lengths, ndim):
    """Halton positions in the centered box ``[-Lx/2, Lx/2] x ...`` (3D, with
    unused dims zeroed)."""
    positions = np.zeros([n, 3], np.float32)
    primes = [2, 3, 5]
    for dim in range(ndim):
        x = _halton_sequence(primes[dim], n)
        positions[:, dim] = (x - 0.5) * box_lengths[dim]
    return positions


def _subrandom_particle_positions(
    nparticles,
    box_vectors,
    ndim,
    method="halton",
    forbidden_xy=None,
    forbidden_radii=None,
):
    """Generate a deterministic list of subrandom particle positions.

    If ``forbidden_xy`` (shape ``(K, 2)``) and ``forbidden_radii`` (shape
    ``(K,)``) are provided, candidates whose XY position falls inside any
    per-atom disk are excluded under the periodic box's minimum-image
    convention -- the downstream LJ-fluid sim is periodic, so a candidate
    whose *direct* XY distance to an obstacle is large but whose *min-image*
    distance is small must be rejected here, otherwise the candidate starts
    inside an obstacle in the periodic frame and the WCA repulsion at sub-A
    separations destabilises the integrator (especially in single precision,
    e.g. WebGPU). The Halton sequence is regenerated at growing length until
    ``nparticles`` non-forbidden positions are available, with a 10x cap.
    """
    box_lengths = [
        box_vectors[d][d].value_in_unit(unit.angstrom) for d in range(3)
    ]

    if forbidden_xy is None or len(forbidden_xy) == 0:
        positions = _generate_halton_positions(nparticles, box_lengths, ndim)
        np.random.shuffle(positions)
        return unit.Quantity(positions, unit.angstrom)

    forbidden_xy = np.asarray(forbidden_xy)
    forbidden_radii = np.asarray(forbidden_radii)
    radii2 = forbidden_radii * forbidden_radii
    box_xy = np.asarray(box_lengths[:2], dtype=float)

    n_needed = int(np.ceil(nparticles * 1.3))
    max_factor = 10
    cap = nparticles * max_factor
    while True:
        positions = _generate_halton_positions(n_needed, box_lengths, ndim)
        diffs = positions[:, None, :2] - forbidden_xy[None, :, :]
        # Wrap to the minimum-image vector under XY periodicity so wrap-around
        # candidates (e.g. lipid near +Lx/2, obstacle near -Lx/2) are also
        # excluded; without this, those candidates start in deep WCA clash.
        diffs -= box_xy * np.round(diffs / box_xy)
        dists2 = np.sum(diffs * diffs, axis=2)
        in_forbidden = (dists2 < radii2[None, :]).any(axis=1)
        kept = positions[~in_forbidden]
        if len(kept) >= nparticles:
            np.random.shuffle(kept)
            return unit.Quantity(kept[:nparticles], unit.angstrom)
        if n_needed >= cap:
            raise RuntimeError(
                f"Could not generate {nparticles} Halton positions outside "
                f"the forbidden region after {n_needed} candidates "
                f"({len(kept)} kept). Increase the box size or check that "
                f"the protein footprint does not exceed the box area."
            )
        n_needed = min(n_needed * 2, cap)


def distributeLipids(
    boxsize,
    resnames,
    sigmas,
    cutoff,
    mass=39.9 * unit.amu,  # argon
    epsilon=0.238 * unit.kilocalories_per_mole,  # argon,
    switch_width=3.4 * unit.angstrom,  # argon
    forbidden_xy=None,
    forbidden_radii=None,
    platform_name=None,
):
    from moleculekit.periodictable import periodictable

    nparticles = len(resnames)

    # Determine Lennard-Jones cutoff.
    cutoff = cutoff * unit.angstrom

    cutoff_type = openmm.NonbondedForce.CutoffPeriodic

    # Create an empty system object.
    system = openmm.System()

    # Periodic box vectors.
    a = unit.Quantity(
        (boxsize[0] * unit.angstrom, 0 * unit.angstrom, 0 * unit.angstrom)
    )
    b = unit.Quantity(
        (0 * unit.angstrom, boxsize[1] * unit.angstrom, 0 * unit.angstrom)
    )
    c = unit.Quantity(
        (0 * unit.angstrom, 0 * unit.angstrom, boxsize[2] * unit.angstrom)
    )
    system.setDefaultPeriodicBoxVectors(a, b, c)

    # Set up periodic nonbonded interactions with a cutoff.
    nb = openmm.NonbondedForce()
    nb.setNonbondedMethod(cutoff_type)
    nb.setCutoffDistance(cutoff)
    nb.setUseDispersionCorrection(True)

    nb.setUseSwitchingFunction(False)
    if switch_width is not None:
        nb.setUseSwitchingFunction(True)
        nb.setSwitchingDistance(cutoff - switch_width)

    for s in sigmas:
        system.addParticle(mass)
        nb.addParticle(0.0 * unit.elementary_charge, s * unit.angstrom, epsilon)

    positions = _subrandom_particle_positions(
        nparticles,
        system.getDefaultPeriodicBoxVectors(),
        2,
        forbidden_xy=forbidden_xy,
        forbidden_radii=forbidden_radii,
    ).value_in_unit(unit.angstrom)

    # Obstacles are added as frozen particles plus a CustomNonbondedForce
    # restricted to lipid-obstacle pairs via an interaction group. Excluding
    # obstacle-obstacle pairs avoids astronomical clashes between protein
    # atoms collapsed onto the same z=0 plane.
    n_obstacles = 0 if forbidden_xy is None else len(forbidden_xy)
    obstacle_indices = []
    if n_obstacles > 0:
        forbidden_xy_arr = np.asarray(forbidden_xy)
        forbidden_radii_arr = np.asarray(forbidden_radii)
        # Treat the per-atom radius as the obstacle's WCA sigma.
        for i in range(n_obstacles):
            idx = system.addParticle(0)  # frozen
            obstacle_indices.append(idx)
            # Inert in the standard NonbondedForce so its index aligns with
            # the System; obstacle-obstacle and obstacle-lipid LJ here are
            # both zero (we use the CustomNonbondedForce below for the
            # obstacle-lipid repulsion).
            nb.addParticle(
                0.0 * unit.elementary_charge,
                1.0 * unit.angstrom,
                0.0 * unit.kilocalories_per_mole,
            )

        # WCA repulsion between lipids and obstacles. Wall strength is
        # matched to the lipid-lipid LJ epsilon so the two potentials live
        # on the same scale.
        wca = openmm.CustomNonbondedForce(
            "step(1.122462048309373*sig - r) * 4*eps*("
            "(sig/r)^12 - (sig/r)^6 + 0.25);"
            "sig = 0.5*(sigma1+sigma2)"
        )
        wca.addGlobalParameter("eps", epsilon.value_in_unit(unit.kilojoules_per_mole))
        wca.addPerParticleParameter("sigma")
        wca.setNonbondedMethod(openmm.CustomNonbondedForce.CutoffPeriodic)
        wca.setCutoffDistance(cutoff)
        wca.setUseLongRangeCorrection(False)
        # Lipid sigmas (in nm)
        for s in sigmas:
            wca.addParticle([s * 0.1])
        # Obstacle sigmas (radii treated as sigma, in nm)
        for r in forbidden_radii_arr:
            wca.addParticle([float(r) * 0.1])
        wca.addInteractionGroup(
            list(range(nparticles)), obstacle_indices
        )
        system.addForce(wca)

    # Add the nonbonded force.
    system.addForce(nb)

    # Append obstacle XY positions at z=0 so positions array aligns with System.
    if n_obstacles > 0:
        obstacle_positions = np.zeros((n_obstacles, 3), dtype=np.float32)
        obstacle_positions[:, :2] = forbidden_xy_arr
        positions = np.vstack([positions, obstacle_positions])

    # Add a restraining potential to keep lipids in z=0 (obstacles are frozen).
    energy_expression = "k * (z^2)"
    force = openmm.CustomExternalForce(energy_expression)
    force.addGlobalParameter("k", 10)
    for particle_index in range(nparticles):
        force.addParticle(particle_index, [])
    system.addForce(force)

    # Create topology.
    topology = app.Topology()
    chain = topology.addChain()
    elems = list(periodictable.keys())  # Just create residues of random dummy elements
    _, idx = np.unique(resnames, return_inverse=True)
    for i in idx:
        element = app.Element.getBySymbol(elems[i])
        residue = topology.addResidue(elems[i], chain)
        topology.addAtom(elems[i], element, residue)
    if n_obstacles > 0:
        obs_chain = topology.addChain()
        carbon = app.Element.getBySymbol("C")
        for _ in range(n_obstacles):
            res = topology.addResidue("OBS", obs_chain)
            topology.addAtom("X", carbon, res)

    topology.setUnitCellDimensions(unit.Quantity(boxsize, unit.angstrom))

    # Simulate it
    from openmm import VerletIntegrator, Platform
    from openmm.app import Simulation
    from openmm.unit import picoseconds, angstrom

    nsteps = 10000

    integrator = VerletIntegrator(0.002 * picoseconds)
    if platform_name is not None:
        simulation = Simulation(
            topology,
            system,
            integrator,
            Platform.getPlatformByName(platform_name),
        )
    else:
        simulation = Simulation(topology, system, integrator)
    simulation.context.setPositions(positions * angstrom)
    simulation.minimizeEnergy()
    simulation.step(nsteps)

    state = simulation.context.getState(getPositions=True, enforcePeriodicBox=False)
    allfinalpos = state.getPositions(asNumpy=True).value_in_unit(angstrom)

    box_xy = np.array([boxsize[0], boxsize[1]])
    if n_obstacles > 0:
        # Obstacles have mass=0 so they cannot have moved. The shift of
        # their XY centroid relative to the input is purely OpenMM's
        # internal wrapping; subtract it to recover the caller's frame.
        initial_obs_com = obstacle_positions[:, :2].mean(axis=0)
        final_obs_com = allfinalpos[nparticles:nparticles + n_obstacles, :2].mean(axis=0)
        allfinalpos[:, :2] -= final_obs_com - initial_obs_com
        anchor_xy = initial_obs_com
    else:
        anchor_xy = np.zeros(2)

    # Wrap each particle's XY into the periodic image closest to the
    # box center so lipids that drifted across a boundary during the
    # LJ sim end up in the same image as the obstacles.
    delta = allfinalpos[:, :2] - anchor_xy
    delta -= box_xy * np.round(delta / box_xy)
    allfinalpos[:, :2] = anchor_xy + delta

    # Return both the final positions (after LJ relaxation) and the initial
    # ones (Halton filter output, guaranteed outside the forbidden disks).
    return allfinalpos[:nparticles], np.asarray(positions[:nparticles])


if __name__ == "__main__":
    pass
