from htmd.parameterization.parameterset import findImproperType, getParameter
import os
import parmed
import numpy as np
import logging


logger = logging.getLogger(__name__)


def getSortedAndUniqueTypes(types, field):
    if field == 'atom_types':
        return np.unique(types)
    if field == 'bond_types':
        for i in range(types.shape[0]):
            types[i] = sorted(types[i])
    elif field == 'angle_types':
        for i in range(types.shape[0]):
            types[i][0], types[i][2] = sorted([types[i][0], types[i][2]])
    elif field == 'dihedral_types':
        for i in range(types.shape[0]):
            if types[i][0] > types[i][3]:
                types[i] = types[i][::-1]
    elif field == 'improper_types' or field == 'improper_periodic_types':
        for i in range(types.shape[0]):
            types[i][0], types[i][1], types[i][3] = sorted([types[i][0], types[i][1], types[i][3]])
    else:
        raise RuntimeError('Invalid field')
    return sorted(list({tuple(row) for row in types.tolist()}))


def getAtomTypeMapping(prm):
    # Make a type mapping for any name != 2 chars in length
    # Because Amber file formats are horrid
    atom_type_map = dict()
    idx = 97
    for atom_type in prm.atom_types:
        atom_type_alias = atom_type
        if len(atom_type) != 2:
            atom_type_alias = "z%c" % idx
            idx += 1
        atom_type_map[atom_type] = atom_type_alias
    return atom_type_map


def mapAtomTypesParameterSet(parameters, typemap):
    from copy import deepcopy, copy
    from parmed.parameters import ParameterSet

    newparams = ParameterSet()
    for type, val in parameters.atom_types.items():
        if type in typemap:
            newparams.atom_types[typemap[type]] = copy(val)

    for f in ('bond_types', 'angle_types', 'dihedral_types', 'improper_types', 'improper_periodic_types'):
        for key in parameters.__dict__[f]:
            newkey = np.array(key, dtype=object)
            newkey = tuple(np.vectorize(typemap.get)(newkey))
            newparams.__dict__[f][newkey] = copy(parameters.__dict__[f][key])
    return newparams


def writeFRCMOD(mol, parameters, filename, typemap=None):
    from htmd.version import version as htmdversion
    if typemap is not None:
        parameters = mapAtomTypesParameterSet(parameters, typemap)
        atomtypes = np.vectorize(typemap.get)(mol.atomtype)
    else:
        atomtypes = mol.atomtype

    f = open(filename, "w")
    f.write("Frcmod generated by HTMD parameterize version {}\n".format(htmdversion()))
    f.write("MASS\n")
    for at in np.unique(atomtypes):
        f.write('{:<2s}  {:>8.3f} {:>10.2f}\n'.format(at, parameters.atom_types[at].mass, 0.))

    f.write('\nBOND\n')
    types = getSortedAndUniqueTypes(atomtypes[mol.bonds], 'bond_types')
    for type in types:
        val = getParameter(type, parameters.bond_types)
        f.write('{:<2s}-{:<2s} {:>10.3f}{:>8.3f}\n'.format(type[0], type[1], val.k, val.req))

    f.write('\nANGL\n')
    types = getSortedAndUniqueTypes(atomtypes[mol.angles], 'angle_types')
    for type in types:
        val = getParameter(type, parameters.angle_types)
        f.write('{:<2s}-{:<2s}-{:<2s} {:>10.3f}{:>9.3f}\n'.format(type[0], type[1], type[2], val.k, val.theteq))

    f.write('\nDIHE\n')
    types = getSortedAndUniqueTypes(atomtypes[mol.dihedrals], 'dihedral_types')
    for type in types:
        val = getParameter(type, parameters.dihedral_types)

        toprint = []
        for term in val:
            if term.phi_k < 1e-7:
                continue
            toprint.append(term)

        # HACK: print at least one dihedral, even if the force constant is 0, otherwise "tleap" is not happy!
        if len(toprint) == 0:
            toprint.append(val[0])

        for i, term in enumerate(toprint):
            # All terms of the same dihedral except the last one should be negative. http://ambermd.org/formats.html#frcmod
            per = term.per
            if i != len(toprint)-1:
                per = -per
            fmt = '{:<2s}-{:<2s}-{:<2s}-{:<2s} {:>4d}{:>15.8f}{:>9.3f}{:>6.1f}    SCEE={:1.1f} SCNB={:1.1f}\n'
            f.write(fmt.format(type[0], type[1], type[2], type[3], 1, term.phi_k, term.phase, per, term.scee, term.scnb))

    f.write('\nIMPR\n')
    types = getSortedAndUniqueTypes(atomtypes[mol.impropers], 'improper_types')
    for type in types:
        type, field = findImproperType(type, parameters)
        val = parameters.__dict__[field][type]
        fmt = '{:<2s}-{:<2s}-{:<2s}-{:<2s}     {:>10.8f}{:>9.3f}{:>6.1f}\n'
        if field == 'improper_periodic_types':
            if val.phi_k < 1e-7:
                continue
            f.write(fmt.format(type[0], type[1], type[2], type[3], val.phi_k, val.phase, val.per))
        elif field == 'improper_types':
            if val.psi_k < 1e-7:
                continue
            f.write(fmt.format(type[0], type[1], type[2], type[3], val.psi_k, val.psi_eq, 1))

    f.write('\nNONB\n')
    # Have to iterate over the types in use, which include cloned types, and map them back
    # to original type (which has the same vdw params), because a copy of a copy won't be in self.nonbonded.
    types = getSortedAndUniqueTypes(atomtypes, 'atom_types')
    for type in types:
        val = parameters.atom_types[type]
        f.write('{:<2s}    {:>10.8f}   {:>10.8f}\n'.format(type, val.rmin, val.epsilon))
    f.write('\n')
    f.close()


def writePRM(mol, parameters, filename):
    from htmd.version import version as htmdversion
    from htmd.util import ensurelist

    # for type, val in parameters.atom_types.items():
    #     if val.epsilon_14 != 1.0:
    #         raise ValueError("Can't express 1-4 electrostatic scaling in Charmm file format")

    f = open(filename, "w")
    print("* prm file built by HTMD parameterize version {}".format(htmdversion()), file=f)
    print("*\n", file=f)

    print("BONDS", file=f)
    types = getSortedAndUniqueTypes(mol.atomtype[mol.bonds], 'bond_types')
    for type in types:
        val = getParameter(type, parameters.bond_types)
        print("%-6s %-6s %8.2f %8.4f" % (type[0], type[1], val.k, val.req), file=f)

    print("\nANGLES", file=f)
    types = getSortedAndUniqueTypes(mol.atomtype[mol.angles], 'angle_types')
    for type in types:
        val = getParameter(type, parameters.angle_types)
        print("%-6s %-6s %-6s %8.2f %8.2f" % (type[0], type[1], type[2], val.k, val.theteq), file=f)

    print("\nDIHEDRALS", file=f)
    types = getSortedAndUniqueTypes(mol.atomtype[mol.dihedrals], 'dihedral_types')
    for type in types:
        val = getParameter(type, parameters.dihedral_types)
        for term in val:
            print("%-6s %-6s %-6s %-6s %12.8f %d %12.8f" % (type[0], type[1], type[2], type[3], term.phi_k, term.per, term.phase), file=f)

    print("\nIMPROPER", file=f)
    types = getSortedAndUniqueTypes(mol.atomtype[mol.impropers], 'improper_types')
    for type in types:
        type, field = findImproperType(type, parameters)
        val = parameters.__dict__[field][type]
        if field == 'improper_periodic_types':
            for term in ensurelist(val):
                print("%-6s %-6s %-6s %-6s %12.8f %d %12.8f" % (type[0], type[1], type[2], type[3], term.phi_k, term.per, term.phase), file=f)
        elif field == 'improper_types':
            print("%-6s %-6s %-6s %-6s %12.8f %d %12.8f" % (type[0], type[1], type[2], type[3], val.psi_k, 0, val.psi_eq), file=f)

    print("\nNONBONDED nbxmod  5 atom cdiel shift vatom vdistance vswitch -", file=f)
    print("cutnb 14.0 ctofnb 12.0 ctonnb 10.0 eps 1.0 e14fac 1.0 wmin 1.5", file=f)
    types = getSortedAndUniqueTypes(mol.atomtype, 'atom_types')
    for type in types:
        val = parameters.atom_types[type]
        if val.epsilon_14 != val.epsilon:
            print("%-6s 0.0000 %8.4f %8.4f 0.0000 %8.4f %8.4f" % (type, val.epsilon, val.rmin, val.epsilon_14, val.rmin_14), file=f)
        else:
            print("%-6s 0.0000 %8.4f %8.4f" % (type, val.epsilon, val.rmin), file=f)
    f.close()


def writeRTF(mol, parameters, netcharge, filename):
    import periodictable
    from htmd.version import version as htmdversion

    f = open(filename, "w")
    print("* Charmm RTF built by HTMD parameterize version {}".format(htmdversion()), file=f)
    print("* ", file=f)
    print("  22     0", file=f)
    types = getSortedAndUniqueTypes(mol.atomtype, 'atom_types')
    for type in types:
        val = parameters.atom_types[type]
        print("MASS %5d %s %8.5f %s" % (val.number, type, val.mass, periodictable.elements[val.atomic_number]), file=f)
    print("\nAUTO ANGLES DIHE\n", file=f)
    assert np.unique(mol.resname).size == 1 # Check if all atoms have the same residue name
    print("RESI %s %8.5f" % (mol.resname[0], netcharge), file=f)
    print("GROUP", file=f)

    maxnamelen = max(4, np.max([len(na) for na in mol.name]))  # Minimum field size of 4
    maxatomtypelen = max(6, np.max([len(at) for at in mol.atomtype]))  # Minimum field size of 6
    for n, a, c in zip(mol.name, mol.atomtype, mol.charge):
        print("ATOM {0:>{1}s} {2:>{3}s} {4:8.6f}".format(n, maxnamelen, a, maxatomtypelen, c), file=f)
    for a in mol.bonds:
        print("BOND {:>4s} {:>4s}".format(*sorted([mol.name[a[0]], mol.name[a[1]]])), file=f)
    for a in mol.impropers:
        print("IMPR %4s %4s %4s %4s" % (mol.name[a[0]], mol.name[a[1]], mol.name[a[2]], mol.name[a[3]]),
              file=f)
    print("PATCH FIRST NONE LAST NONE", file=f)
    print("\nEND", file=f)
    f.close()


def writeParameters(outdir, mol, parameters, method, netcharge, original_coords=None):
    os.makedirs(outdir, exist_ok=True)

    typemap = None
    extensions = ('mol2', 'pdb', 'coor')

    if method == 'CGenFF_2b6':
        extensions += ('psf', 'rtf', 'prm')

        # TODO: remove?
        f = open(os.path.join(outdir, "input.namd"), "w")
        tmp = '''parameters mol.prm
paraTypeCharmm on
coordinates mol.pdb
bincoordinates mol.coor
temperature 0
timestep 0
1-4scaling 1.0
exclude scaled1-4
outputname .out
outputenergies 1
structure mol.psf
cutoff 20.
switching off
stepsPerCycle 1
rigidbonds none
cellBasisVector1 50. 0. 0.
cellBasisVector2 0. 50. 0.
cellBasisVector3 0. 0. 50.
run 0'''
        print(tmp, file=f)
        f.close()

    elif method in ('GAFF', 'GAFF2'):
        # types need to be remapped because Amber FRCMOD format limits the type to characters
        # writeFrcmod does this on the fly and returns a mapping that needs to be applied to the mol
        # TODO: get rid of this mapping
        frcFile = os.path.join(outdir, 'mol.frcmod')
        typemap = getAtomTypeMapping(parameters)
        writeFRCMOD(mol, parameters, frcFile, typemap=typemap)
        logger.info('Write FRCMOD file: %s' % frcFile)

        tleapFile = os.path.join(outdir, 'tleap.in')
        with open(tleapFile, 'w') as file_:
            file_.write('loadAmberParams mol.frcmod\n')
            file_.write('A = loadMol2 mol.mol2\n')
            file_.write('saveAmberParm A structure.prmtop mol.crd\n')
            file_.write('quit\n')
        logger.info('Write tleap input file: %s' % tleapFile)

        # TODO: remove?
        f = open(os.path.join(outdir, "input.namd"), "w")
        tmp = '''parmfile structure.prmtop
amber on
coordinates mol.pdb
bincoordinates mol.coor
temperature 0
timestep 0
1-4scaling 0.83333333
exclude scaled1-4
outputname .out
outputenergies 1
cutoff 20.
switching off
stepsPerCycle 1
rigidbonds none
cellBasisVector1 50. 0. 0.
cellBasisVector2 0. 50. 0.
cellBasisVector3 0. 0. 50.
run 0'''
        print(tmp, file=f)
        f.close()

    else:
        raise NotImplementedError

    def remapAtomTypes(mol):
        tmpmol = mol
        if typemap is not None:
            tmpmol = mol.copy()
            tmpmol.atomtype[:] = [typemap[atomtype] for atomtype in mol.atomtype]
        return tmpmol

    tmpmol = remapAtomTypes(mol)

    for ext in extensions:
        file_ = os.path.join(outdir, "mol." + ext)
        if ext == 'prm':
            writePRM(mol, parameters, file_)
        elif ext == 'rtf':
            writeRTF(mol, parameters, netcharge, file_)
        else:
            tmpmol.write(file_)
        logger.info('Write %s file: %s' % (ext.upper(), file_))

    if original_coords is not None:
        molFile = os.path.join(outdir, 'mol-orig.mol2')
        tmpmol.coords = original_coords
        tmpmol.write(molFile)
        logger.info('Write MOL2 file (with original coordinates): {}'.format(molFile))
