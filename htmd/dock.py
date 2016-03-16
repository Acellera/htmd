# (c) 2015-2016 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import numpy as np
from tempfile import NamedTemporaryFile
from os import path
from htmd.home import home
import os
import shutil
from subprocess import call, check_output
from htmd.util import tempname
from htmd.molecule.molecule import Molecule


def dock(protein, ligand, center=None, extent=None):
    """ Molecular docking, using Vina

    If centre and extent are not provided, docking will be performed over the whole protein

    Parameters
    ----------
    protein
        Molecule object representing the receptor
    ligand
        Molecule object representing the ligand to dock
    center
        3-vec centre of of the search bounding box (optional)
    extent
        3-vec linear extent of the search bounding box (optional)

    Returns
    -------
    poses
        Molecule object representing the N<10 best poses
    scores
        3x num_poses matrix containing kcal, rmsd lb, rmsd ub

    Examples
    --------
    >>> poses = dock(protein, ligand)
    >>> poses = dock(protein, ligand, center=[ 10., 5., 12. ], extent=[ 15., 15., 15. ] )

    """
    if np.size(protein.coords, 2) != 1 or np.size(ligand.coords, 2) != 1:
        raise NameError('Protein and ligand Molecules should be single frames')

    buffer = 10.  # Angstrom buffer around protein for whole-protein docking
    c_min = np.min(protein.coords, 0).reshape((1, 3))[0]
    c_max = np.max(protein.coords, 0).reshape((1, 3))[0]

    if center is None:
        center = (buffer + (c_max + c_min)) / 2
    if extent is None:
        extent = (c_max - c_min) + buffer

    # babel -i pdb protein.pdb  -o pdbqt protein.pdbqt -xr
    # babel -i pdb ligand.pdb   -o pdbqt ligand.pdbqt -xhn
    # vina --ligand ligand.pdbqt --receptor protein.pdbqt --center_x 0. --center_y 0. --center_z 0. --size_x 60. --size_y 60. --size_z 60 --exhaustiveness 10
    # babel -m -i pdbqt ligand_out.pdbqt -o pdb out_.pdb -xhn

    protein_pdb = tempname(suffix=".pdb")
    ligand_pdb = tempname(suffix=".pdb")
    output_pdb = tempname(suffix="_.pdb")
    output_prefix = path.splitext(output_pdb)[0]

    protein_pdbqt = tempname(suffix=".pdbqt")
    ligand_pdbqt = tempname(suffix=".pdbqt")
    output_pdbqt = tempname(suffix=".pdbqt")

    protein.write(protein_pdb)
    ligand.write(ligand_pdb)

    try:
        import platform
        suffix=""
        if platform.system() == "Windows":
          suffix=".exe"
        vinaexe = shutil.which( platform.system() + "-vina" + suffix, mode=os.X_OK )
    except:
        raise NameError('Could not find vina, or no execute permissions are given')
    try:
        babelexe = shutil.which('htmd_babel', mode=os.X_OK)
    except:
        raise NameError('Could not find babel, or no execute permissions are given')
    #babelexe = path.join(home(), 'bin', 'htmd_babel')
    #vinaexe = path.join(home(), 'bin', 'htmd_vina')

    #if not os.access(babelexe, os.X_OK):
    #    raise NameError('Could not find ' + babelexe + ' or no execute permissions are given')
    #if not os.access(vinaexe, os.X_OK):
    #    raise NameError('Could not find ' + vinaexe + ' or no execute permissions are given')

    #from IPython.core.debugger import Tracer
    #Tracer()()
    call([babelexe, '-i', 'pdb', protein_pdb, '-o', 'pdbqt', protein_pdbqt, '-xr'])
    call([babelexe, '-i', 'pdb', ligand_pdb, '-o', 'pdbqt', ligand_pdbqt, '-xhn'])

    if not path.isfile(ligand_pdbqt):
        raise NameError('Ligand could not be converted to PDBQT')
    if not path.isfile(protein_pdbqt):
        raise NameError('Protein could not be converted to PDBQT')

    call([vinaexe, '--receptor', protein_pdbqt, '--ligand', ligand_pdbqt, '--out', output_pdbqt,
          '--center_x', str(center[0]), '--center_y', str(center[1]), '--center_z', str(center[2]),
          '--size_x', str(extent[0]), '--size_y', str(extent[1]), '--size_z', str(extent[2])])

    call([babelexe, '-m', '-i', 'pdbqt', output_pdbqt, '-o', 'pdb', output_pdb, '-xhn'])

    scoring = np.zeros((0,3))
    coords = []
    idx = 1
    name = '{}{}.pdb'.format(output_prefix, idx)
    while path.isfile(name):
        # First get the scoring
        kcal = float(check_output('grep "VINA RESULT" ' + name + ' | awk \'{print $4}\'', shell=True).decode('ascii').strip())
        rmsdlb = float(check_output('grep "VINA RESULT" ' + name + ' | awk \'{print $5}\'', shell=True).decode('ascii').strip())
        rmsdub = float(check_output('grep "VINA RESULT" ' + name + ' | awk \'{print $6}\'', shell=True).decode('ascii').strip())
        scoring = np.append(scoring, np.array([[float(kcal), float(rmsdlb), float(rmsdub)]]), axis=0)
        next_pose = Molecule(name)
        os.remove(name)
        c = next_pose.coords
        co = c.copy()
        natoms = len(ligand.name)

        for idx_i in range(natoms):
            for idx_j in range(natoms):
                if ligand.name[idx_i] == next_pose.name[idx_j]:
                    co[idx_i, :, :] = c[idx_j, :, :]

        '''if idx == 1:
            coords = co
        else:
            coords = np.append(coords, co, axis=2)'''
        coords.append(co)

        idx += 1
        name = '{}{}.pdb'.format(output_prefix, idx)

    poses = []
    for i, c in enumerate(coords):
        l = ligand.copy()
        l.viewname = 'Pose {}'.format(i)
        l.coords = c
        poses.append(l)

    os.remove(protein_pdb)
    os.remove(ligand_pdb)
    os.remove(protein_pdbqt)
    os.remove(ligand_pdbqt)
    os.remove(output_pdbqt)

    return poses, scoring

if __name__ == "__main__":
    from htmd.molecule.molecule import Molecule
    from os import path
    from htmd.home import home
    protein = Molecule(path.join(home(), 'data', 'docking', 'protein.pdb'))
    ligand = Molecule(path.join(home(), 'data', 'docking', 'ligand.pdb'))
    poses, scoring = dock(protein, ligand)

