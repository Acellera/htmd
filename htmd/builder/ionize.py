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

# Info about all ions, keyed by their residue name in AMBER/OpenMM ion libraries:
# charge, name, AMBER name and CHARMM name (None when CHARMM has no such ion).
_ions = {
    "NA": (1, "sodium", "Na+", "SOD"),
    "MG": (2, "magnesium", "Mg2+", "MG"),
    "ZN": (2, "zinc", "Zn2+", "ZN"),
    "K": (1, "potassium", "K+", "POT"),
    "CS": (1, "cesium", "Cs+", "CES"),
    "CA": (2, "calcium", "Ca2+", "CAL"),
    "CL": (-1, "chloride", "Cl-", "CLA"),
    "F": (-1, "fluoride", "F-", None),
    "BR": (-1, "bromide", "Br-", None),
    "IOD": (-1, "iodide", "I-", None),
}

_ion_aliases = {}
for _key, (_charge, _fullname, _ambername, _charmmname) in _ions.items():
    for _alias in (_key, _fullname, _ambername, _charmmname):
        if _alias is not None:
            _ion_aliases[_alias.upper()] = _key
_ion_aliases["I"] = "IOD"


def _ionElement(ion):
    """Element symbol of a canonical ion, from its AMBER name: "Mg2+" -> "Mg"."""
    return _ions[ion][2].rstrip("+-0123456789")


# Every supported cation/anion pair, keyed by formula: "MgCl2" -> MG, CL
_salt_formulas = {
    f"{_ionElement(cat)}{_ionElement(an)}{'2' if qc == 2 else ''}": (cat, an)
    for cat, (qc, *_) in _ions.items()
    if qc > 0
    for an, (qa, *_) in _ions.items()
    if qa < 0
}


def ionize(
    mol: Molecule,
    netcharge: float,
    nwater: int,
    neutralize: bool = True,
    saltconc: float | None = None,
    cation: str | None = None,
    anion: str | None = None,
    ff: str | None = None,
):
    """Calculate the number of ions required to neutralize or salt a system.

    Computes the anion and cation counts needed to neutralize the system's net
    charge, and optionally adds extra ions to reach a target salt concentration.
    Ion stoichiometry is inferred from the ion charges (e.g. MgCl2 for Mg2+/Cl-).

    Parameters
    ----------
    mol : Molecule
        The molecular system (used to count existing ions when *saltconc* is set).
    netcharge : float
        Net charge of the molecular system (e.g. from partial charges).
    nwater : int
        Number of water molecules in the system (used to estimate box volume for
        salt concentration calculations).
    neutralize : bool
        If True, return only the counts needed for neutralization. If *saltconc*
        is not None, this is automatically set to False.
    saltconc : float, optional
        Target salt concentration in molar. If given, additional ions are added on
        top of the neutralizing ions to reach this concentration.
    cation : str, optional
        Cation residue name. Accepts canonical keys (``"NA"``, ``"K"``, ``"MG"``,
        ``"CA"``, ``"CS"``), full names (``"sodium"``), AMBER names (``"Na+"``),
        or CHARMM names (``"SOD"``). If None, defaults to ``"NA"``.
    anion : str, optional
        Anion residue name. Accepts canonical keys (``"CL"``), full names
        (``"chloride"``), AMBER names (``"Cl-"``), or CHARMM names (``"CLA"``).
        If None, defaults to ``"CL"``.
    ff : str, optional
        Deprecated. No longer needed; pass standard ion names via *cation* / *anion*.

    Returns
    -------
    tuple
        When *neutralize* is True: ``(anion, cation, nanion, ncation)``.
        When *saltconc* is given: ``(anion, cation, anion, cation, nanion, ncation)``.

    Raises
    ------
    NameError
        If the ion charges cannot neutralize the system, if no waters are found,
        or if the requested salt concentration is already exceeded.
    """
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

    cationstoich, anionstoich, num = _saltUnits(
        mol.resname, cation, anion, saltconc, nwater
    )

    logger.info(
        f"Adding {nanion} anions + {ncation} cations for neutralizing and {(cationstoich + anionstoich) * num} ions for the given salt concentration {saltconc} M."
    )
    ncation += cationstoich * num
    nanion += anionstoich * num
    return anion, cation, anion, cation, nanion, ncation


def _saltUnits(resnames, cation, anion, saltconc, nwater):
    """Formula of the salt made of ``cation`` and ``anion`` and how many formula units
    to add for ``saltconc`` (M) in ``nwater`` waters, net of the units already present
    in ``resnames``. Returns ``(cationstoich, anionstoich, num)``."""
    cationcharge = _ionGetCharge(cation)
    anioncharge = _ionGetCharge(anion)
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
    cation_names = [n for n in (cation, _ions[cation][3]) if n]  # AMBER, CHARMM
    anion_names = [n for n in (anion, _ions[anion][3]) if n]
    ncations_exist = np.sum(np.isin(resnames, cation_names)) / cationstoich
    nanions_exist = np.sum(np.isin(resnames, anion_names)) / anionstoich
    exist = int(min(ncations_exist, nanions_exist))

    num = int(np.floor(0.5 + 0.0187 * saltconc * nwater))
    num -= exist
    if num < 0:
        raise RuntimeError(
            f"Salt concentration in system is already higher than the requested concentration {saltconc} M. Please remove some ions or change salt concentration and try again."
        )
    return cationstoich, anionstoich, num


def resolveSalts(
    salts: list | None = None,
    saltconc: float | None = 0,
    cation: str | None = None,
    anion: str | None = None,
) -> list:
    """Normalize the salts of a build into ``[(cation, anion, concentration), ...]``.

    Builders take either a list of salts or the single-salt shorthand
    ``saltconc`` / ``saltcation`` / ``saltanion``, not both.

    Parameters
    ----------
    salts : list, optional
        ``[(formula, concentration), ...]`` with concentrations in molar, e.g.
        ``[("KCl", 0.15), ("MgCl2", 0.001)]``. A formula pairs Na, K, Cs, Mg, Ca or
        Zn with Cl, F, Br or I (e.g. ``"NaBr"``, ``"CaI2"``), in any case. The
        first salt also neutralizes the system.
    saltconc : float
        Concentration (M) of the single salt, when ``salts`` is not given.
    cation : str, optional
        Cation of the single salt. Defaults to ``"NA"``. Deprecated in favour of
        ``salts``.
    anion : str, optional
        Anion of the single salt. Defaults to ``"CL"``. Deprecated in favour of
        ``salts``.

    Returns
    -------
    salts : list
        ``[(cation, anion, concentration), ...]`` with canonical ion keys.
    """
    if not salts:
        if cation or anion:
            logger.warning(
                "saltcation / saltanion are deprecated. Pass the salt as a formula "
                'instead, e.g. salts=[("KCl", 0.15)].'
            )
        cation, anion = _resolve_ion(cation or "NA"), _resolve_ion(anion or "CL")
        if _ionGetCharge(cation) <= 0 or _ionGetCharge(anion) >= 0:
            raise ValueError(f"{cation} is not a cation or {anion} not an anion.")
        return [(cation, anion, float(saltconc or 0))]
    if saltconc or cation or anion:
        raise ValueError(
            "Give either salts or saltconc / saltcation / saltanion, not both."
        )
    resolved = []
    for formula, conc in salts:
        ions = next(
            (v for k, v in _salt_formulas.items() if k.upper() == formula.upper()),
            None,
        )
        if ions is None:
            raise ValueError(
                f"Unknown salt '{formula}'. Supported: {', '.join(_salt_formulas)}."
            )
        if any(ions == (c, a) for c, a, _ in resolved):
            raise ValueError(f"Salt '{formula}' is given twice.")
        if conc < 0:
            raise ValueError(f"Salt '{formula}' has a negative concentration.")
        resolved.append((*ions, float(conc)))
    return resolved


def ionizeSalts(resnames, netcharge: float, nwater: int, salts: list) -> dict:
    """Number of each ion to add to neutralize a system and reach several salt
    concentrations.

    The first salt's ions neutralize the system, as in :func:`ionize`. Every salt
    then adds the formula units for its concentration on top, like a single
    ``saltconc``.

    Parameters
    ----------
    resnames : array
        Residue names of the system's atoms, to count ions already present.
    netcharge : float
        Net charge of the system.
    nwater : int
        Number of water molecules, which sets the volume for the concentrations.
    salts : list
        ``[(cation, anion, concentration), ...]`` with canonical ion names, as
        returned by :func:`resolveSalts`.

    Returns
    -------
    counts : dict
        Number of ions to add, keyed by canonical ion name.
    """
    counts = {}
    for i, (cation, anion, conc) in enumerate(salts):
        nanion = ncation = 0
        if i == 0:
            _, _, nanion, ncation = ionize(
                None, netcharge, nwater, neutralize=True, cation=cation, anion=anion
            )
        cationstoich, anionstoich, num = _saltUnits(
            resnames, cation, anion, conc, nwater
        )
        if conc > 0 and num == 0:
            logger.warning(
                f"{conc} M of {cation}/{anion} is less than one formula unit in "
                f"{nwater} waters, so none is added. Use a larger box to include it."
            )
        counts[cation] = counts.get(cation, 0) + ncation + cationstoich * num
        counts[anion] = counts.get(anion, 0) + nanion + anionstoich * num
    logger.info(f"Adding ions: {counts}")
    return counts


def _resolve_ion(name):
    """Resolve any ion alias (canonical key, full name, AMBER name, CHARMM name) to the canonical key."""
    key = _ion_aliases.get(name.upper())
    if key is None:
        valid = ", ".join(f"{k} ({v[1]}, {v[2]}, {v[3]})" for k, v in _ions.items())
        raise NameError(
            f"Ion '{name}' not recognized. Valid ions and their aliases: {valid}"
        )
    return key


def _ionGetCharge(ion):
    ion = _resolve_ion(ion)
    return _ions[ion][0]


def ionizePlace(
    solvent_mol: Molecule,
    solute_mol: Molecule | None,
    anion_resname: str,
    cation_resname: str,
    anion_name: str,
    cation_name: str,
    nanion: int,
    ncation: int,
    dfrom: float = 5,
    dbetween: float = 5,
    segname: str | None = None,
) -> Molecule:
    """Place a given number of negative and positive ions in the solvent.

    Uses farthest point sampling to produce well-spaced ion positions and
    randomly assigns ion types for uniform spatial mixing.
    Replaces water molecules as long as they respect the given distance criteria.

    Parameters
    ----------
    solvent_mol : Molecule
        The solvent (water) molecule.
    solute_mol : Molecule or None
        The solute molecule used for minimum-distance filtering. If None or
        empty, no minimum distance from solute is enforced.
    anion_resname : str
        Resname of the added anions.
    cation_resname : str
        Resname of the added cations.
    anion_name : str
        Name of the added anions.
    cation_name : str
        Name of the added cations.
    nanion : int
        Number of anions to add.
    ncation : int
        Number of cations to add.
    dfrom : float
        Minimum distance of ions from solute molecule in Angstroms.
    dbetween : float
        Minimum distance between ions in Angstroms.
    segname : str, optional
        Segment name to assign to the ions.

    Returns
    -------
    mol : Molecule
        The solvent molecule with ions replacing some water molecules.
    """
    return _placeIons(
        solvent_mol,
        solute_mol,
        [(anion_resname, anion_name, nanion), (cation_resname, cation_name, ncation)],
        dfrom,
        dbetween,
    )


def placeIons(
    solvent_mol: Molecule,
    solute_mol: Molecule | None,
    counts: dict,
    dfrom: float = 5,
    dbetween: float = 5,
) -> Molecule:
    """Place ions of any number of types in the solvent, replacing waters.

    Like :func:`ionizePlace`, but for the ``{ion: count}`` of :func:`ionizeSalts`.
    Each ion's residue and atom name is its canonical name.

    Parameters
    ----------
    solvent_mol : Molecule
        The solvent (water) molecule.
    solute_mol : Molecule or None
        The solute molecule used for minimum-distance filtering.
    counts : dict
        Number of ions to add, keyed by canonical ion name.
    dfrom : float
        Minimum distance of ions from solute molecule in Angstroms.
    dbetween : float
        Minimum distance between ions in Angstroms.

    Returns
    -------
    mol : Molecule
        The solvent molecule with ions replacing some water molecules.
    """
    return _placeIons(solvent_mol, solute_mol, _ionSpecies(counts), dfrom, dbetween)


def _ionSpecies(counts):
    """``[(resname, name, count), ...]`` for ``counts``, anions first. With one salt
    this is the (anion, cation) order of :func:`ionizePlace`, so a single salt is
    labelled exactly as it always was. The atom name is the upper-case element, as
    in the AMBER and OpenMM ion libraries (``IOD`` holds atom ``I``)."""
    order = sorted(counts, key=lambda ion: _ionGetCharge(ion) > 0)
    return [(ion, _ionElement(ion).upper(), counts[ion]) for ion in order]


def _pickIonSites(wat_coords, solute_coords, nions, dfrom, dbetween):
    """Indices into ``wat_coords`` of ``nions`` well-spaced water oxygens at least
    ``dfrom`` from the solute and ``dbetween`` from each other."""
    # Filter by minimum distance from solute using KD-tree
    if solute_coords is not None and len(solute_coords) > 0:
        solute_tree = cKDTree(solute_coords)
        dists_to_solute, _ = solute_tree.query(wat_coords)
        valid_mask = dists_to_solute >= dfrom
        valid_local = np.where(valid_mask)[0]
    else:
        valid_local = np.arange(len(wat_coords))

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

    return valid_local[selected]


def _ionLabels(species):
    """Shuffled per-ion indices into ``species``. Since farthest point sampling
    guarantees well-spaced positions regardless of type, random labelling gives
    uniform spatial mixing independent of the ion ratio."""
    labels = np.concatenate(
        [np.full(n, i) for i, (_, _, n) in enumerate(species)]
    ).astype(int)
    np.random.shuffle(labels)
    return labels


def _placeIons(solvent_mol, solute_mol, species, dfrom, dbetween):
    newmol = solvent_mol.copy()

    nions = sum(n for _, _, n in species)
    logger.debug(f"Min distance of ions from molecule: {dfrom}A")
    logger.debug(f"Min distance between ions: {dbetween}A")
    logger.debug(f"Placing {nions} ions: {species}")

    if nions == 0:
        return newmol

    betabackup = newmol.beta.copy()
    newmol.set("beta", sequenceID((newmol.resid, newmol.insertion, newmol.segid)))

    # Find water oxygens to replace with ions
    wat_oh = newmol.atomselect("noh and water", indexes=True)

    if len(wat_oh) == 0:
        raise NameError("No water molecules found in the solvent to replace with ions.")

    solute_coords = None
    if solute_mol is not None and solute_mol.numAtoms > 0:
        solute_coords = solute_mol.coords[:, :, 0]
    picked = _pickIonSites(
        newmol.coords[wat_oh, :, 0], solute_coords, nions, dfrom, dbetween
    )
    ionlist = wat_oh[picked]
    type_labels = _ionLabels(species)

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

    for i, (resname, name, _) in enumerate(species):
        ions_mol.resname[type_labels == i] = resname
        ions_mol.name[type_labels == i] = name

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
