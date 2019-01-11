import numpy as np
import re
from htmd.parameterization.parameterset import _ATOM_TYPE_REG_EX
import logging

logger = logging.getLogger(__name__)


def _guessElement(name):
    import re
    name = re.sub('[0-9]*$', '', name)
    name = name.lower().capitalize()
    return name


def _guessMass(element):
    from htmdmol.molecule import vdw
    return vdw.massByElement(element)


# TODO: fix guessElement
def readRTF(filename):
    f = open(filename, 'r')
    lines = f.readlines()
    f.close()

    types = []
    mass_by_type = dict()
    element_by_type = dict()
    type_by_name = dict()
    type_by_index = []
    index_by_name = dict()
    names = []
    charge_by_name = dict()
    bonds = []
    improper_indices = []
    typeindex_by_type = dict()
    netcharge = 0.

    aidx = 0
    for l in lines:
        if l.startswith("MASS "):
            k = l.split()
            at = k[2]
            mass_by_type[at] = float(k[3])
            element_by_type[at] = k[4]
            typeindex_by_type[at] = int(k[1])
            types.append(at)
        elif l.startswith("RESI "):
            k = l.split()
            netcharge = float(k[2])
        elif l.startswith("ATOM "):
            k = l.split()
            names.append(k[1])
            index_by_name[k[1]] = aidx
            type_by_index.append(k[2])
            type_by_name[k[1]] = k[2]
            charge_by_name[k[1]] = float(k[3])
            aidx += 1
        elif l.startswith("BOND "):
            k = l.split()
            bonds.append([index_by_name[k[1]], index_by_name[k[2]]])
        elif l.startswith("IMPR "):
            k = l.split()
            improper_indices.append([index_by_name[k[1]], index_by_name[k[2]], index_by_name[k[3]], index_by_name[k[4]]])

    # if there weren't any "MASS" lines, we need to guess them
    typeindex = 4000
    for idx in range(len(names)):
        atype = type_by_index[idx]
        name = names[idx]
        if atype not in element_by_type:
            element_by_type[atype] = _guessElement(name)
            logger.info("Guessing element %s for atom %s type %s" % (element_by_type[atype], name, atype))
        if atype not in mass_by_type:
            mass_by_type[atype] = _guessMass(element_by_type[atype])

        if atype not in typeindex_by_type:
            typeindex_by_type[atype] = typeindex
            typeindex += 1
        if atype not in types:
            types.append(atype)

    names = np.array(names, dtype=object)
    type_by_index = np.array(type_by_index, dtype=object)
    element_by_idx = np.array([element_by_type[t].lower().capitalize() for t in type_by_index], dtype=object)
    charge_by_idx = np.array([charge_by_name[n] for n in names], dtype=np.float32)
    mass_by_idx = np.array([mass_by_type[t] for t in type_by_index], dtype=np.float32)

    improper_indices = np.array(improper_indices).astype(np.uint32)
    if improper_indices.ndim == 1:
        improper_indices = improper_indices[:, np.newaxis]

    for type_ in type_by_index:
        if re.match(_ATOM_TYPE_REG_EX, type_):
            raise ValueError('Atom type %s is incompatible. It cannot finish with "x" + number!'.format(type_))

    return names, element_by_idx, type_by_index, charge_by_idx, mass_by_idx, improper_indices


def readPREPI(mol, prepi):
    f = open(prepi, 'r')
    lines = f.readlines()
    f.close()
    f = lines

    # the prepi has the atoms re-ordered. Reorder the info based on the order in the mol
    index_by_name = {name: i for i, name in enumerate(mol.name)}

    types = []
    names = np.array(['' for _ in range(mol.numAtoms)], dtype=object)
    type_by_idx = np.array(['' for _ in range(mol.numAtoms)], dtype=object)
    charge_by_idx = np.zeros(mol.numAtoms, dtype=np.float32)

    if f[4].split()[1] != 'INT':
        raise ValueError('Invalid prepi format line 5')
    if f[5].strip() != "CORRECT     OMIT DU   BEG":
        raise ValueError('Invalid prepi format line 6')

    ctr = 10
    while f[ctr].strip() != '':
        ff = f[ctr].split()
        ff[1] = ff[1].upper()
        idx = index_by_name[ff[1]]
        names[idx] = ff[1]
        type_by_idx[idx] = ff[2]
        charge_by_idx[idx] = float(ff[10])
        ctr += 1

    # Read improper section
    with open(prepi) as file:
        text = file.read()
    improper_indices = []
    impropers = re.search('^IMPROPER\n(.+)\n\n', text, re.MULTILINE | re.DOTALL)  # extract improper section
    if impropers:
        impropers = impropers.group(1).split('\n')  # array of improper lines
        impropers = [improper.split() for improper in impropers]  # impropers by names
        for improper in impropers:
            idx = [index_by_name[name.upper()] for name in improper]  # conv atom name to indices
            improper_indices.append(idx)

    improper_indices = np.array(improper_indices).astype(np.uint32)
    if improper_indices.ndim == 1:
        improper_indices = improper_indices[:, np.newaxis]

    for type_ in type_by_idx:
        if re.match(_ATOM_TYPE_REG_EX, type_):
            raise ValueError('Atom type %s is incompatible. It cannot finish with "x" + number!'.format(type_))

    return names, type_by_idx, charge_by_idx, improper_indices


def readFRCMOD(atomtypes, frcmod):
    from periodictable import elements
    mass2element = {e.mass: e.symbol for e in list(elements._element.values())[1:]}

    # Read MASS section
    with open(frcmod) as file:
        text = file.read()
    section = re.search('^MASS\n(.+?)\n\n', text, re.MULTILINE | re.DOTALL)

    mass_by_atomtype = {line.split()[0]: float(line.split()[1]) for line in section.group(1).split('\n')}
    element_by_atomtype = {}
    for at in mass_by_atomtype:
        for m in mass2element:
            if np.isclose(mass_by_atomtype[at], m, atol=1e-1):
                element_by_atomtype[at] = mass2element[m]

    element_by_idx = np.array([element_by_atomtype[at] for at in atomtypes]).astype(np.object)
    mass_by_idx = np.array([mass_by_atomtype[at] for at in atomtypes]).astype(np.float32)

    return mass_by_idx, element_by_idx
