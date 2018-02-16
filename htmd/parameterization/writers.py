from htmd.parameterization.util import _qm_method_name
from htmd.parameterization.fftype import FFTypeMethod
import os
import parmed
import numpy as np
import logging


logger = logging.getLogger(__name__)


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


def mapAtomTypesParameterSet(prm, typemap):
    from copy import deepcopy, copy
    from parmed.parameters import ParameterSet

    newprm = ParameterSet()
    for type, val in prm.atom_types.items():
        if type in typemap:
            newprm.atom_types[typemap[type]] = copy(val)

    for f in ('bond_types', 'angle_types', 'dihedral_types', 'improper_types', 'improper_periodic_types'):
        for key in prm.__dict__[f]:
            newkey = np.array(key)
            newkey = tuple(np.vectorize(typemap.get)(newkey))
            newprm.__dict__[f][newkey] = copy(prm.__dict__[f][key])


def writeFRCMOD(prm, typemap, outfile):
    from htmd.version import version as htmdversion
    prm = mapAtomTypesParameterSet(prm, typemap)
    parmed.amber.AmberParameterSet.write(prm, outfile, title='FRCMOD built by HTMD parameterize version {}'.format(htmdversion()))


def writePRM(prm, filename):
    from htmd.version import version as htmdversion

    for type, val in prm.dihedral_types.items():
        if val.epsilon_14 != 1.0:
            raise ValueError("Can't express 1-4 electrostatic scaling in Charmm file format")

    f = open(filename, "w")
    print("* prm file built by HTMD parameterize version {}".format(htmdversion()), file=f)
    print("*\n", file=f)

    print("BONDS", file=f)
    for type, val in prm.bond_types.items():
        print("%-6s %-6s %8.2f %8.4f" % (type[0], type[1], val.k, val.req), file=f)

    print("\nANGLES", file=f)
    for type, val in prm.angle_types.items():
        print("%-6s %-6s %-6s %8.2f %8.2f" % (type[0], type[1], type[2], val.k, val.theteq), file=f)

    print("\nDIHEDRALS", file=f)
    for type, val in prm.dihedral_types.items():
        print("%-6s %-6s %-6s %-6s %12.8f %d %12.8f" % (type[0], type[1], type[2], type[3], val.phi_k, val.per, val.phase), file=f)

    print("\nIMPROPER", file=f)
    for type, val in prm.improper_types.items():
        print("%-6s %-6s %-6s %-6s %12.8f %d %12.8f" % (type[0], type[1], type[2], type[3], val.psi_k, 0, val.psi_eq), file=f)
    for type, val in prm.improper_periodic_types.items():
        print("%-6s %-6s %-6s %-6s %12.8f %d %12.8f" % (type[0], type[1], type[2], type[3], val.phi_k, val.per, val.phase), file=f)

    print("\nNONBONDED nbxmod  5 atom cdiel shift vatom vdistance vswitch -", file=f)
    print("cutnb 14.0 ctofnb 12.0 ctonnb 10.0 eps 1.0 e14fac 1.0 wmin 1.5", file=f)
    for type, val in prm.atom_types.items():
        if val.epsilon_14 is not None:
            # Charmm prm stores rmin/2
            print("%-6s 0.0000 %8.4f %8.4f 0.0000 %8.4f %8.4f" % (type, val.epsilon, val.rmin, val.epsilon_14, val.rmin_14), file=f)
        else:
            print("%-6s 0.0000 %8.4f %8.4f" % (type, val.epsilon, val.rmin), file=f)
    f.close()


def writeRTF(mol, prm, netcharge, filename):
    import periodictable
    from htmd.version import version as htmdversion

    f = open(filename, "w")
    print("* Charmm RTF built by HTMD parameterize version {}".format(htmdversion()), file=f)
    print("* ", file=f)
    print("  22     0", file=f)
    for type, val in prm.atom_types.items():
        print("MASS %5d %s %8.5f %s" % (typeindex_by_type[type], type, val.mass, periodictable.elements[val.atomic_number]), file=f)
    print("\nAUTO ANGLES DIHE\n", file=f)
    print("RESI  MOL %8.5f" % netcharge, file=f)
    print("GROUP", file=f)
    for n, a, c in (mol.name, mol.atomtype, mol.charge):
        print("ATOM %4s %6s %8.6f" % (n, a, c), file=f)
    for a in mol.bonds:
        print("BOND %4s %4s" % (mol.names[a[0]], mol.names[a[1]]), file=f)
    for a in mol.impropers:
        print("IMPR %4s %4s %4s %4s" % (mol.names[a[0]], mol.names[a[1]], mol.names[a[2]], mol.names[a[3]]),
              file=f)
    print("PATCH FIRST NONE LAST NONE", file=f)
    print("\nEND", file=f)
    f.close()


def writeParameters(mol, prm, qm, method, outdir, original_molecule=None):

    paramDir = os.path.join(outdir, 'parameters', method.name, _qm_method_name(qm))
    os.makedirs(paramDir, exist_ok=True)

    typemap = None
    extensions = ('mol2', 'pdb', 'coor')

    if method == FFTypeMethod.CGenFF_2b6:
        extensions += ('psf', 'rtf', 'prm')

        # TODO: remove?
        f = open(os.path.join(paramDir, "input.namd"), "w")
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

    elif method in (FFTypeMethod.GAFF, FFTypeMethod.GAFF2):
        # types need to be remapped because Amber FRCMOD format limits the type to characters
        # writeFrcmod does this on the fly and returns a mapping that needs to be applied to the mol
        # TODO: get rid of this mapping
        frcFile = os.path.join(paramDir, 'mol.frcmod')
        typemap = getAtomTypeMapping(prm)
        writeFRCMOD(prm, typemap, frcFile)
        logger.info('Write FRCMOD file: %s' % frcFile)

        tleapFile = os.path.join(paramDir, 'tleap.in')
        with open(tleapFile, 'w') as file_:
            file_.write('loadAmberParams mol.frcmod\n')
            file_.write('A = loadMol2 mol.mol2\n')
            file_.write('saveAmberParm A structure.prmtop mol.crd\n')
            file_.write('quit\n')
        logger.info('Write tleap input file: %s' % tleapFile)

        # TODO: remove?
        f = open(os.path.join(paramDir, "input.namd"), "w")
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
        file_ = os.path.join(paramDir, "mol." + ext)
        if ext == 'prm':
            writePRM()
        elif ext == 'rtf':
            writeRTF()
        else:
            tmpmol.write(file_)
        logger.info('Write %s file: %s' % (ext.upper(), file_))

    if original_molecule:
        molFile = os.path.join(paramDir, 'mol-orig.mol2')
        tmpmol = remapAtomTypes(original_molecule)
        tmpmol.write(molFile)
        logger.info('Write MOL2 file (with original coordinates): {}'.format(molFile))