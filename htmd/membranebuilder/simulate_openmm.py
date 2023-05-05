# (c) 2015-2022 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import simtk.openmm as mm
from sys import stdout
from simtk import unit
from moleculekit.molecule import Molecule
from simtk.openmm import app


def _printEnergies(sim, text):
    state = sim.context.getState(getEnergy=True)
    pot = state.getPotentialEnergy()
    kin = state.getKineticEnergy()
    print(f"{text} {pot + kin} {pot} {kin}")


def _getPlatform(plat, device):
    platform = mm.Platform.getPlatformByName(plat)
    prop = None
    if plat == "CUDA":
        prop = {"CudaPrecision": "single", "CudaDeviceIndex": str(device)}
    return platform, prop


def equilibrateSystem(
    pdbfile,
    prmtopfile,
    outpdb,
    numsteps=30000,
    minimplatform="CPU",
    equilplatform="CUDA",
    device=0,
    temp=300,
    minimize=500000,
    minimizetol=100,
):
    pdb = Molecule(pdbfile)
    watcoo = pdb.get("coords", "water")
    celld = (watcoo.max(axis=0) - watcoo.min(axis=0)).squeeze()

    pdb = app.PDBFile(pdbfile)
    structure = app.AmberPrmtopFile(prmtopfile)

    a = unit.Quantity((celld[0] * unit.angstrom, 0 * unit.angstrom, 0 * unit.angstrom))
    b = unit.Quantity((0 * unit.angstrom, celld[1] * unit.angstrom, 0 * unit.angstrom))
    c = unit.Quantity((0 * unit.angstrom, 0 * unit.angstrom, celld[2] * unit.angstrom))
    structure.box_vectors = (a, b, c)
    system = structure.createSystem(
        nonbondedMethod=app.PME,
        nonbondedCutoff=1 * unit.nanometer,
        constraints=app.HBonds,
    )
    system.setDefaultPeriodicBoxVectors(a, b, c)

    system.addForce(mm.MonteCarloBarostat(1 * unit.atmospheres, temp * unit.kelvin, 25))

    platform, prop = _getPlatform(minimplatform, device)
    integrator = mm.LangevinIntegrator(
        temp * unit.kelvin, 1 / unit.picosecond, 0.002 * unit.picoseconds
    )
    simulation = app.Simulation(structure.topology, system, integrator, platform, prop)
    simulation.context.setPositions(pdb.positions)

    # Perform minimization
    _printEnergies(simulation, "Energy before minimization")
    simulation.minimizeEnergy(
        tolerance=minimizetol * unit.kilojoule / unit.mole, maxIterations=minimize
    )
    _printEnergies(simulation, "Energy after minimization")

    # Copy coordinates from miminization context to equilibration context
    state = simulation.context.getState(getPositions=True)
    pos = state.getPositions()

    # Set up the equilibration simulation
    platform, prop = _getPlatform(equilplatform, device)
    integrator = mm.LangevinIntegrator(
        temp * unit.kelvin, 1 / unit.picosecond, 0.002 * unit.picoseconds
    )
    simulation = app.Simulation(structure.topology, system, integrator, platform, prop)
    simulation.context.setPositions(pos)

    # Generate initial velocities
    simulation.context.setVelocitiesToTemperature(temp)

    from htmd.membranebuilder.pdbreporter import PDBReporter

    simulation.reporters.append(PDBReporter(outpdb, numsteps, enforcePeriodicBox=False))
    simulation.reporters.append(
        app.StateDataReporter(
            stdout,
            int(numsteps / 10),
            step=True,
            potentialEnergy=True,
            kineticEnergy=True,
            totalEnergy=True,
            temperature=True,
            progress=True,
            speed=True,
            remainingTime=True,
            totalSteps=numsteps,
            separator="\t",
        )
    )
    simulation.step(numsteps)
