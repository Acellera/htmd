from htmd.molecule.util import boundingBox
import numpy as np
import ctypes
import htmd
import os

_order = ('hydrophobic', 'aromatic', 'hbond_acceptor', 'hbond_donor', 'positive_ionizable', 'negative_ionizable', 'metal')
libdir = htmd.home(libDir=True)
occupancylib = ctypes.cdll.LoadLibrary(os.path.join(libdir, "occupancy_ext.so"))


def getVoxelDescriptors(mol, buffer, voxelsize=1):
    """ Calculate descriptors of atom properties for voxels in a grid bounding the Molecule object.

    Constructs a bounding box around Molecule with some buffer space. Then it

    Parameters
    ----------
    mol :
        A Molecule object.
    buffer : float
        The buffer space to add to the bounding box. This adds zeros to the grid around the protein so that properties
        which are at the edge of the box can be found in the center of one. Should be usually set to boxsize/2.
    voxelsize : float
        The voxel size in A

    Returns
    -------
    features : np.ndarray
        A list of boxes containing voxels and their properties
    centers : np.ndarray
        A list of the centers of all boxes
    """
    properties = _getAtomtypePropertiesPDBQT(mol)

    # Calculate the bbox and the number of voxels
    [bbm, bbM] = boundingBox(mol)
    bbm -= buffer
    bbM += buffer
    N = np.ceil((bbM - bbm) / voxelsize).astype(int)

    # Calculate for each channel the atom sigmas
    multisigmas = np.zeros([mol.numAtoms, len(_order) + 1])  # add 1 channel for occupancy
    sigmas = _getSigmas(mol)
    for i, p in enumerate(_order):
        multisigmas[properties[p], i] = sigmas[properties[p]]
    multisigmas[:, -1] = sigmas  # occupancy

    # Calculate occupancies and centers
    occs, centers = _getGridDescriptors(mol, bbm, N, multisigmas, voxelsize)

    occs = np.swapaxes(occs, 2, 3)
    occs = np.swapaxes(occs, 1, 2)
    occs = np.swapaxes(occs, 0, 1)

    return occs, centers


def _getAtomtypePropertiesPDBQT(mol):
    """ Matches PDBQT atom types to specific properties
    ('hydrophobic', 'aromatic', 'hbond_acceptor', 'hbond_donor', 'positive_ionizable', 'negative_ionizable', 'metal')

    Parameters
    ----------
    mol :
        A Molecule object

    Returns
    -------
    properties : dict
        A dictionary of atom masks for the Molecule showing for each property (key) which atoms (value) belong to it.
    """
    """ AutoDock 4 atom types http://autodock.scripps.edu/faqs-help/faq/faqsection_view?section=Scientific%20Questions
    #        --     ----  -----  -------  --------  ---  ---  -  --  -- --
    atom_par H      2.00  0.020   0.0000   0.00051  0.0  0.0  0  -1  -1  3    # Non H-bonding Hydrogen
    atom_par HD     2.00  0.020   0.0000   0.00051  0.0  0.0  2  -1  -1  3    # Donor 1 H-bond Hydrogen
    atom_par HS     2.00  0.020   0.0000   0.00051  0.0  0.0  1  -1  -1  3    # Donor S Spherical Hydrogen
    atom_par C      4.00  0.150  33.5103  -0.00143  0.0  0.0  0  -1  -1  0    # Non H-bonding Aliphatic Carbon
    atom_par A      4.00  0.150  33.5103  -0.00052  0.0  0.0  0  -1  -1  0    # Non H-bonding Aromatic Carbon
    atom_par N      3.50  0.160  22.4493  -0.00162  0.0  0.0  0  -1  -1  1    # Non H-bonding Nitrogen
    atom_par NA     3.50  0.160  22.4493  -0.00162  1.9  5.0  4  -1  -1  1    # Acceptor 1 H-bond Nitrogen
    atom_par NS     3.50  0.160  22.4493  -0.00162  1.9  5.0  3  -1  -1  1    # Acceptor S Spherical Nitrogen
    atom_par OA     3.20  0.200  17.1573  -0.00251  1.9  5.0  5  -1  -1  2    # Acceptor 2 H-bonds Oxygen
    atom_par OS     3.20  0.200  17.1573  -0.00251  1.9  5.0  3  -1  -1  2    # Acceptor S Spherical Oxygen
    atom_par F      3.09  0.080  15.4480  -0.00110  0.0  0.0  0  -1  -1  4    # Non H-bonding Fluorine
    atom_par Mg     1.30  0.875   1.5600  -0.00110  0.0  0.0  0  -1  -1  4    # Non H-bonding Magnesium
    atom_par MG     1.30  0.875   1.5600  -0.00110  0.0  0.0  0  -1  -1  4    # Non H-bonding Magnesium
    atom_par P      4.20  0.200  38.7924  -0.00110  0.0  0.0  0  -1  -1  5    # Non H-bonding Phosphorus
    atom_par SA     4.00  0.200  33.5103  -0.00214  2.5  1.0  5  -1  -1  6    # Acceptor 2 H-bonds Sulphur
    atom_par S      4.00  0.200  33.5103  -0.00214  0.0  0.0  0  -1  -1  6    # Non H-bonding Sulphur
    atom_par Cl     4.09  0.276  35.8235  -0.00110  0.0  0.0  0  -1  -1  4    # Non H-bonding Chlorine
    atom_par CL     4.09  0.276  35.8235  -0.00110  0.0  0.0  0  -1  -1  4    # Non H-bonding Chlorine
    atom_par Ca     1.98  0.550   2.7700  -0.00110  0.0  0.0  0  -1  -1  4    # Non H-bonding Calcium
    atom_par CA     1.98  0.550   2.7700  -0.00110  0.0  0.0  0  -1  -1  4    # Non H-bonding Calcium
    atom_par Mn     1.30  0.875   2.1400  -0.00110  0.0  0.0  0  -1  -1  4    # Non H-bonding Manganese
    atom_par MN     1.30  0.875   2.1400  -0.00110  0.0  0.0  0  -1  -1  4    # Non H-bonding Manganese
    atom_par Fe     1.30  0.010   1.8400  -0.00110  0.0  0.0  0  -1  -1  4    # Non H-bonding Iron
    atom_par FE     1.30  0.010   1.8400  -0.00110  0.0  0.0  0  -1  -1  4    # Non H-bonding Iron
    atom_par Zn     1.48  0.550   1.7000  -0.00110  0.0  0.0  0  -1  -1  4    # Non H-bonding Zinc
    atom_par ZN     1.48  0.550   1.7000  -0.00110  0.0  0.0  0  -1  -1  4    # Non H-bonding Zinc
    atom_par Br     4.33  0.389  42.5661  -0.00110  0.0  0.0  0  -1  -1  4    # Non H-bonding Bromine
    atom_par BR     4.33  0.389  42.5661  -0.00110  0.0  0.0  0  -1  -1  4    # Non H-bonding Bromine
    atom_par I      4.72  0.550  55.0585  -0.00110  0.0  0.0  0  -1  -1  4    # Non H-bonding Iodine
    """
    props = dict()
    elements = np.array([el.upper() for el in mol.element])

    props['hydrophobic'] = (elements == 'C') | (elements == 'A')
    props['aromatic'] = elements == 'A'
    props['hbond_acceptor'] = (elements == 'NA') | (elements == 'NS') | (elements == 'OA') | (elements == 'OS') | (
    elements == 'SA')
    props['hbond_donor'] = _findDonors(mol, mol._getBonds())
    props['positive_ionizable'] = mol.charge > 0
    props['negative_ionizable'] = mol.charge < 0
    props['metal'] = (elements == 'MG') | (elements == 'ZN') | (elements == 'MN') | \
                     (elements == 'CA') | (elements == 'FE')

    return props


def _findDonors(mol, bonds):
    """ Finds atoms connected to HD and HS atoms to find the heavy atom donor (O or N)

    Parameters
    ----------
    mol :
        A Molecule object
    bonds : np.ndarray
        An array of atom bonds

    Returns
    -------
    donors : np.ndarray
        Boolean array indicating the donor heavy atoms in Mol
    """
    donors = np.zeros(mol.numAtoms, dtype=bool)
    hydrogens = np.where((mol.element == 'HD') | (mol.element == 'HS'))[0]
    for h in hydrogens:
        partners = bonds[bonds[:, 0] == h, 1]
        partners = np.hstack((partners, bonds[bonds[:, 1] == h, 0]))
        for p in partners:
            if mol.name[p][0] == 'N' or mol.name[p][0] == 'O':
                donors[p] = True
    return donors


def _getSigmas(mol):
    """ Returns the VDW radii of each atom in the molecule

    Parameters
    ----------
    mol :
        A Molecule object

    Returns
    -------
    sigmas : np.ndarray
        The radius of each atom in Molecule
    """
    # NOTE: this should improved. Element or name? More than first characeter? FF based?
    # e.g. Cl gets transformed wrongly to C. Also Mn and Zn are not present in dictionary
    sigmas = {"H": 1.0, "C": 1.5, "N": 1.4, "O": 1.3, "F": 1.2, "P": 2.0, "S": 1.9,
              "Cl": 2.5, "D": 4.0, "A": 1.5}
    sig = np.zeros(mol.numAtoms)
    for a in range(mol.numAtoms):  # TODO: STEFAN - Should be doable faster, although it's not really time-consuming
        elem = mol.element[a][0]
        if elem in sigmas:
            sigma = sigmas[elem]
        else:
            print('unknown element -', mol.element[a], '- at atom index ', a)
            sigma = 1.3
        sig[a] = sigma
    return sig


def _getGridDescriptors(mol, llc, N, channelsigmas, resolution):
    """ Calls the C code to calculate the voxels values for each property.

    Parameters
    ----------
    mol
    llc
    N
    channelsigmas
    resolution

    Returns
    -------
    occupancies
    centers
    """
    center = np.array(llc + 0.5 * resolution * np.array(N))
    xrange = [llc[0] + resolution * x for x in range(0, N[0])]
    yrange = [llc[1] + resolution * x for x in range(0, N[1])]
    zrange = [llc[2] + resolution * x for x in range(0, N[2])]
    centers = np.zeros((N[0], N[1], N[2], 3))
    for i, x in enumerate(xrange):
        for j, y in enumerate(yrange):
            for k, z in enumerate(zrange):
                centers[i, j, k, :] = np.array([x, y, z])
    centers1D = centers.reshape(np.prod(N), 3)
    nchannels = channelsigmas.shape[1]
    occus = np.zeros((centers1D.shape[0], nchannels))
    coords = np.squeeze(mol.coords[:, :, 0])

    occupancylib.descriptor_ext(centers1D.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                       coords.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                       channelsigmas.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                       occus.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                       center.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                       ctypes.c_double(np.max(0.5 * np.array(N) * resolution)),  # max distance
                       ctypes.c_int(occus.shape[0]),  # n of centers
                       ctypes.c_int(coords.shape[0]),  # n of atoms
                       ctypes.c_int(nchannels))  # n of channels
    return occus.reshape((N[0], N[1], N[2], nchannels)), centers

if __name__ == '__main__':
    from htmd.molecule.molecule import Molecule
    from htmd.home import home
    import os
    import numpy as np
    resOcc, resCent = getVoxelDescriptors(Molecule('3PTB'), buffer=8, voxelsize=1)
    refOcc = np.load(os.path.join(home(), 'data', 'test-voxeldescriptors', '3PTB_occ.npy'))
    refCent = np.load(os.path.join(home(), 'data', 'test-voxeldescriptors', '3PTB_center.npy'))
    assert np.allclose(resOcc, refOcc)
    assert np.allclose(resCent, refCent)
