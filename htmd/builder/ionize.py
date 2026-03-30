# (c) 2015-2022 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import numpy as np
from scipy.spatial import cKDTree
from moleculekit.molecule import Molecule
from moleculekit.util import sequenceID
import logging

logger = logging.getLogger(__name__)

# Info about all ions, with charge, name, amber and charmm names (may
# need residues too). Only the first column is currently used.
_ions = {
    "NA": (1, "sodium", "Na+", "SOD"),
    "MG": (2, "magnesium", "Mg2+", "MG"),
    "ZN": (2, "zinc", "Zn2+", "ZN"),
    "K": (1, "potassium", "K+", "POT"),
    "CS": (1, "cesium", "Cs+", "CES"),
    "CA": (2, "calcium", "Ca2+", "CAL"),
    "CL": (-1, "chloride", "Cl-", "CLA"),
}

_ion_aliases = {}
for _key, (_charge, _fullname, _ambername, _charmmname) in _ions.items():
    for _alias in (_key, _fullname, _ambername, _charmmname):
        _ion_aliases[_alias.upper()] = _key


def ionize(
    mol,
    netcharge,
    nwater,
    neutralize=True,
    saltconc=None,
    cation=None,
    anion=None,
    ff=None,
):
    if ff:
        logger.warning(
            "The ff option is no longer needed and is deprecated. "
            "Please use standard names (e.g. NA, CL, K) for `cation` and `anion` arguments."
        )

    if cation is None:
        cation = "NA"
    else:
        cation = _resolve_ion(cation)
    if anion is None:
        anion = "CL"
    else:
        anion = _resolve_ion(anion)

    if saltconc is not None:
        neutralize = False

    cationcharge = _ionGetCharge(cation)
    anioncharge = _ionGetCharge(anion)

    roundnetcharge = int(np.round(netcharge))
    done = False
    ncation = 0
    nanion = 0

    # First we neutralize the system
    if netcharge > 0:
        if anioncharge == -1:
            nanion = roundnetcharge
            ncation = 0
        elif anioncharge == -2:
            if np.mod(roundnetcharge, 2) == 0:  # Total charge is even
                nanion = int(roundnetcharge / 2)
                ncation = 0
            else:  # The charge is odd
                if cationcharge == 1:
                    nanion = 1 + int(np.floor(roundnetcharge / 2))
                    ncation = 1
                else:
                    raise NameError(
                        "Failed to neutralize because the cation and anion charges cannot reach 0"
                    )
    elif netcharge < 0:
        roundnetchargeabs = abs(roundnetcharge)
        if cationcharge == 1:
            nanion = 0
            ncation = roundnetchargeabs
        elif cationcharge == 2:
            if np.mod(roundnetchargeabs, 2) == 0:
                nanion = 0
                ncation = int(roundnetchargeabs / 2)
            else:
                if anioncharge == -1:
                    nanion = 1
                    ncation = 1 + int(np.floor(roundnetchargeabs / 2))
                else:
                    raise NameError(
                        "Failed to neutralize because the cation and anion charges cannot reach 0"
                    )
    else:
        logger.info("The system is already neutral.")
        done = True

    if not done:
        newtotalcharge = nanion * anioncharge + ncation * cationcharge + roundnetcharge
        if newtotalcharge != 0:
            raise NameError("Could not neutralize the system")
    if neutralize:
        return anion, cation, nanion, ncation

    # If we want a specific salt concentration, guess chemical formula
    if cationcharge == 1 and anioncharge == -1:  # e.g., NaCl, Kcl, ...
        cationstoich = 1
        anionstoich = 1
    elif cationcharge == 2 and anioncharge == -1:  # e.g., MgCl2
        cationstoich = 1
        anionstoich = 2
    elif cationcharge == 1 and anioncharge == -2:
        cationstoich = 2
        anionstoich = 1
    elif cationcharge == 2 and anioncharge == -2:
        cationstoich = 1
        anionstoich = 1
    else:
        raise NameError("Unsupported ion charge; cannot guess chamical formula.")

    # Count existing salt concentration to not add on top of existing conc
    # TODO: This will break on multi-atom ions! Will require a bit fancier code to calculate individual molecules
    cation_names = [cation, _ions[cation][3]]  # AMBER and CHARMM names
    anion_names = [anion, _ions[anion][3]]
    ncations_exist = np.sum(np.isin(mol.resname, cation_names)) / cationstoich
    nanions_exist = np.sum(np.isin(mol.resname, anion_names)) / anionstoich
    exist = int(min(ncations_exist, nanions_exist))

    num = int(np.floor(0.5 + 0.0187 * saltconc * nwater))
    num -= exist
    if num < 0:
        raise RuntimeError(
            f"Salt concentration in system is already higher than the requested concentration {saltconc} M. Please remove some ions or change salt concentration and try again."
        )

    logger.info(
        f"Adding {nanion} anions + {ncation} cations for neutralizing and {(cationstoich + anionstoich) * num} ions for the given salt concentration {saltconc} M."
    )
    ncation += cationstoich * num
    nanion += anionstoich * num
    return anion, cation, anion, cation, nanion, ncation


def _resolve_ion(name):
    """Resolve any ion alias (canonical key, full name, AMBER name, CHARMM name) to the canonical key."""
    key = _ion_aliases.get(name.upper())
    if key is None:
        valid = ", ".join(
            f"{k} ({v[1]}, {v[2]}, {v[3]})" for k, v in _ions.items()
        )
        raise NameError(
            f"Ion '{name}' not recognized. Valid ions and their aliases: {valid}"
        )
    return key


def _ionGetCharge(ion):
    ion = _resolve_ion(ion)
    return _ions[ion][0]


def ionizePlace(
    solvent_mol,
    solute_mol,
    anion_resname,
    cation_resname,
    anion_name,
    cation_name,
    nanion,
    ncation,
    dfrom=5,
    dbetween=5,
    segname=None,
):
    """Place a given number of negative and positive ions in the solvent.

    Uses farthest point sampling to produce well-spaced ion positions and
    randomly assigns ion types for uniform spatial mixing.
    Replaces water molecules as long as they respect the given distance criteria.

    Parameters
    ----------
    solvent_mol : :class:`Molecule <moleculekit.molecule.Molecule>` object
        The solvent (water) molecule
    solute_mol : :class:`Molecule <moleculekit.molecule.Molecule>` object or None
        The solute molecule used for minimum-distance filtering. If None or
        empty, no minimum distance from solute is enforced.
    anion_resname : str
        Resname of the added anions
    cation_resname : str
        Resname of the added cations
    anion_name : str
        Name of the added anions
    cation_name : str
        Name of the added cations
    nanion : int
        Number of anions to add
    ncation : int
        Number of cations to add
    dfrom : float
        Min distance of ions from solute molecule
    dbetween : float
        Min distance between ions
    segname : str
        Segment name to add

    Returns
    -------
    mol : :class:`Molecule <moleculekit.molecule.Molecule>` object
        The solvent molecule with ions replacing some water molecules
    """

    newmol = solvent_mol.copy()

    logger.debug(f"Min distance of ions from molecule: {dfrom}A")
    logger.debug(f"Min distance between ions: {dbetween}A")
    logger.debug(f"Placing {nanion} anions and {ncation} cations.")

    if (nanion + ncation) == 0:
        return newmol

    nions = nanion + ncation

    betabackup = newmol.beta.copy()
    newmol.set("beta", sequenceID((newmol.resid, newmol.insertion, newmol.segid)))

    # Find water oxygens to replace with ions
    wat_oh = newmol.atomselect("noh and water", indexes=True)

    if len(wat_oh) == 0:
        raise NameError("No water molecules found in the solvent to replace with ions.")

    wat_coords = newmol.coords[wat_oh, :, 0]

    # Filter by minimum distance from solute using KD-tree
    if solute_mol is not None and solute_mol.numAtoms > 0:
        solute_tree = cKDTree(solute_mol.coords[:, :, 0])
        dists_to_solute, _ = solute_tree.query(wat_coords)
        valid_mask = dists_to_solute >= dfrom
        valid_local = np.where(valid_mask)[0]
    else:
        valid_local = np.arange(len(wat_oh))

    if len(valid_local) == 0:
        raise NameError(
            f"No waters could be found further than {dfrom} from other molecules to be replaced by ions. You might need to "
            "solvate with a bigger box or disable the ionize property when building."
        )

    valid_coords = wat_coords[valid_local]

    # Farthest point sampling: iteratively pick the water farthest from all
    # previously placed ions.  This is deterministic, requires no retries,
    # and naturally maximises the minimum inter-ion distance.
    min_dists = np.full(len(valid_coords), np.inf)
    selected = []
    used = np.zeros(len(valid_coords), dtype=bool)

    for i in range(nions):
        idx = np.argmax(min_dists)
        if i > 0 and min_dists[idx] < dbetween:
            raise NameError(
                f"Could not place all ions while maintaining minimum distance of {dbetween} A between ions. Try decreasing the between parameter, "
                "decreasing ion concentration or making a larger water box."
            )
        selected.append(idx)
        used[idx] = True

        # Update min distances from the newly placed ion
        dists_to_new = np.linalg.norm(valid_coords - valid_coords[idx], axis=1)
        min_dists = np.minimum(min_dists, dists_to_new)
        min_dists[used] = -np.inf

    # Map selected indices back to atom indices in newmol
    ionlist = [wat_oh[valid_local[s]] for s in selected]

    # Randomly assign ion types to FPS positions.  Since FPS guarantees
    # well-spaced positions regardless of type, random labelling gives
    # uniform spatial mixing independent of the anion/cation ratio.
    type_labels = np.array(["anion"] * nanion + ["cation"] * ncation)
    np.random.shuffle(type_labels)

    # Delete waters but keep their coordinates
    waterpos = np.atleast_2d(newmol.get("coords", ionlist))
    betasel = np.zeros(newmol.numAtoms, dtype=bool)
    for b in newmol.beta[ionlist]:
        betasel |= newmol.beta == b

    sel = np.where(betasel)[0]
    newmol.remove(sel, _logger=False)
    betabackup = np.delete(betabackup, sel)

    # Build all ions as a single Molecule and append once
    ions_mol = Molecule()
    ions_mol.empty(nions)
    ions_mol.chain[:] = "I"
    ions_mol.segid[:] = "I"
    ions_mol.coords = waterpos[:, :, np.newaxis]

    anion_mask = type_labels == "anion"
    ions_mol.name[anion_mask] = anion_name
    ions_mol.resname[anion_mask] = anion_resname
    ions_mol.name[~anion_mask] = cation_name
    ions_mol.resname[~anion_mask] = cation_resname

    base_resid = newmol.resid[-1] + 1
    ions_mol.resid[:] = np.arange(base_resid, base_resid + nions)

    newmol.append(ions_mol)

    # Restoring the original betas
    newmol.beta[: len(betabackup)] = betabackup
    return newmol


def _getSegname(mol, segname):
    defsegnames = [
        "ION",
        "ION1",
        "ION3",
        "ION4",
        "ION5",
        "ION6",
        "ION7",
        "ION8",
        "ION9",
        "IN1",
        "IN2",
        "IN3",
        "IN4",
        "IN5",
        "IN6",
        "IN7",
        "IN8",
        "IN9",
    ]

    # Make sure segname is not taken
    if segname is None:
        for i in range(len(defsegnames)):
            sel = mol.segid == defsegnames[i]
            if not np.any(sel):
                segname = defsegnames[i]
                break
        if segname is None:
            raise NameError("All default segnames taken! Provide your own.")
    else:
        sel = mol.segid == segname
        if np.any(sel):
            raise NameError("Segname already taken. Provide a different one.")
