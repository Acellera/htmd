# (c) 2015-2018 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import simtk.openmm as mm
from sys import stdout
from simtk import unit
from moleculekit.molecule import Molecule
from simtk.openmm import app


defaultCharmmFiles = [
    'top_all36_prot.rtf',
    'par_all36_prot.prm',
    'top_all36_na.rtf',
    'par_all36_na.prm',
    'top_all36_carb.rtf',
    'par_all36_carb.prm',
    'top_all36_lipid.rtf',
    'par_all36_lipid.prm',
    'top_all36_cgenff.rtf',
    'par_all36_cgenff.prm',
    'toppar_all36_prot_retinol.str',
    'toppar_all36_na_rna_modified.str',
    'toppar_all36_carb_glycopeptide.str',
    'toppar_all36_prot_fluoro_alkanes.str',
    'toppar_all36_prot_na_combined.str',
    'toppar_all36_prot_heme.str',
    'toppar_all36_lipid_bacterial.str',
    'toppar_all36_lipid_miscellaneous.str',
    'toppar_all36_lipid_cholesterol.str',
    'toppar_all36_lipid_yeast.str',
    'toppar_all36_lipid_sphingo.str',
    'toppar_all36_lipid_inositol.str',
    'toppar_all36_lipid_cardiolipin.str',
    'toppar_all36_lipid_detergent.str',
    'toppar_all36_lipid_lps.str',
    'toppar_water_ions.str',
    'toppar_dum_noble_gases.str',
    'toppar_all36_na_nad_ppi.str',
    'toppar_all36_carb_glycolipid.str',
    'toppar_all36_carb_imlab.str'
]


def _readCharmmParameters(folder, listfile=None):
    from glob import glob
    import os
    if listfile is None:
        files = glob(os.path.join(folder, '*'))
    elif isinstance(listfile, list):
        files = [os.path.join(folder, f) for f in listfile]
    return app.CharmmParameterSet(*files)


def _printEnergies(sim, text):
    state = sim.context.getState(getEnergy=True)
    pot = state.getPotentialEnergy()
    kin = state.getKineticEnergy()
    print("{} {} {} {}".format(text, pot + kin, pot, kin))


def _getPlatform(plat, device):
    platform = mm.Platform.getPlatformByName(plat)
    prop = None
    if plat == 'CUDA':
        prop = {'CudaPrecision': 'single', 'CudaDeviceIndex': str(device)}
    return platform, prop


def equilibrateSystem(pdbfile, psffile, outpdb, numsteps=30000, minimplatform='CPU', equilplatform='CUDA', device=0,
                      temp=300, minimize=500000, minimizetol=100, charmmfolder=None):
    pdb = Molecule(pdbfile)
    watcoo = pdb.get('coords', 'water')
    celld = unit.Quantity((watcoo.max(axis=0) - watcoo.min(axis=0)).squeeze(), unit=unit.angstrom)

    psf = app.CharmmPsfFile(psffile)
    pdb = app.PDBFile(pdbfile)
    psf.setBox(celld[0], celld[1], celld[2])

    params = _readCharmmParameters(charmmfolder, defaultCharmmFiles)

    system = psf.createSystem(params, nonbondedMethod=app.PME,
                              nonbondedCutoff=1*unit.nanometer,
                              constraints=app.HBonds)
    system.addForce(mm.MonteCarloBarostat(1*unit.atmospheres, temp*unit.kelvin, 25))

    platform, prop = _getPlatform(minimplatform, device)
    integrator = mm.LangevinIntegrator(temp*unit.kelvin, 1/unit.picosecond, 0.002*unit.picoseconds)
    simulation = app.Simulation(psf.topology, system, integrator, platform, prop)
    simulation.context.setPositions(pdb.positions)

    # Perform minimization
    _printEnergies(simulation, 'Energy before minimization')
    simulation.minimizeEnergy(tolerance=minimizetol*unit.kilojoule/unit.mole, maxIterations=minimize)
    _printEnergies(simulation, 'Energy after minimization')

    # Copy coordinates from miminization context to equilibration context
    state = simulation.context.getState(getPositions=True)
    pos = state.getPositions()

    # Set up the equilibration simulation
    platform, prop = _getPlatform(equilplatform, device)
    integrator = mm.LangevinIntegrator(temp*unit.kelvin, 1/unit.picosecond, 0.002*unit.picoseconds)
    simulation = app.Simulation(psf.topology, system, integrator, platform, prop)
    simulation.context.setPositions(pos)

    # Generate initial velocities
    simulation.context.setVelocitiesToTemperature(temp)

    from htmd.membranebuilder.pdbreporter import PDBReporter
    simulation.reporters.append(PDBReporter(outpdb, numsteps, enforcePeriodicBox=False))
    simulation.reporters.append(app.StateDataReporter(stdout, int(numsteps/10), step=True,
        potentialEnergy=True, kineticEnergy=True, totalEnergy=True, temperature=True,
        progress=True, speed=True, remainingTime=True, totalSteps=numsteps, separator='\t'))
    simulation.step(numsteps)

