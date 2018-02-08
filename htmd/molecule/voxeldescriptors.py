# (c) 2015-2017 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
from htmd.molecule.util import boundingBox
from htmd.molecule import vdw
import numpy as np
import ctypes
import htmd.home
import os
from numba import cuda, jit
import numba
from math import sqrt, exp

_order = ('hydrophobic', 'aromatic', 'hbond_acceptor', 'hbond_donor', 'positive_ionizable',
          'negative_ionizable', 'metal', 'occupancies')
libdir = htmd.home.home(libDir=True)
occupancylib = ctypes.cdll.LoadLibrary(os.path.join(libdir, "occupancy_ext.so"))


def getVoxelDescriptors(mol, usercenters=None, voxelsize=1, buffer=0, channels=None, method='C'):
    """ Calculate descriptors of atom properties for voxels in a grid bounding the Molecule object.

    Constructs a bounding box around Molecule with some buffer space. Then it computes
    pharmacophoric-like descriptors on a defined grid.

    Parameters
    ----------
    mol :
        A Molecule object. Needs to be read from Autodock 4 .pdbqt format
    usercenters : np.ndarray
        A 2D array specifying the centers of the voxels. If None is given, it will discretize the bounding box of the
        Molecule plus any buffer space requested into voxels of voxelsize.
    voxelsize : float
        The voxel size in A
    buffer : float
        The buffer space to add to the bounding box. This adds zeros to the grid around the protein so that properties
        which are at the edge of the box can be found in the center of one. Should be usually set to localimagesize/2.
    channels : np.ndarray
        A 2D array of size (mol.numAtoms, nchannels) where nchannels is the number of channels we want to have.
        Each column i then has True (or a float) in the rows of the atoms which belong to channel i and False (or 0) 
        otherwise. Such boolean arrays can be obtained for example by using mol.atomselect.
        If the array is boolean, each atom will get assigned its VdW radius. If the array is float, these floats will 
        be used as the corresponding atom radii. Make sure the numpy array is of dtype=bool if passing boolean values.
        If no channels are given, the default ('hydrophobic', 'aromatic', 'hbond_acceptor', 'hbond_donor', 
        'positive_ionizable', 'negative_ionizable', 'metal', 'occupancies') channels will be used.

    Returns
    -------
    features : np.ndarray
        A 2D array of size (centers, channels) containing the effect of each channel in the voxel with that center. 
    centers : np.ndarray
        A list of the centers of all boxes
    N : np.ndarray
        Is returned only when no user centers are passed. It corresponds to the number of centers in each of the x,y,z
        dimensions
    method : str
        Voxel descriptors can be calculated either with our C implementation or CUDA or NUMBA implementations.

    Examples
    --------
    >>> mol = Molecule('3PTB')
    >>> mol.filter('protein')
    >>> features, centers, N = getVoxelDescriptors(mol, buffer=8)
    >>> # Features can be reshaped to a 4D array (3D for each grid center in xyz, 1D for the properties) like this:
    >>> features = features.reshape(N[0], N[1], N[2], features.shape[1])
    >>> # The user can provide his own centers
    >>> features, centers = getVoxelDescriptors(mol, usercenters=[[0, 0, 0], [16, 24, -5]], buffer=8)
    """
    if channels is None:
        channels = _getAtomtypePropertiesPDBQT(mol)

    if channels.dtype == bool:
        # Calculate for each channel the atom sigmas
        sigmas = _getRadii(mol)
        channels = sigmas[:, np.newaxis] * channels.astype(float)

    N = None
    if usercenters is None:
        # Calculate the bbox and the number of voxels
        [bbm, bbM] = boundingBox(mol)
        bbm -= buffer
        bbM += buffer
        N = np.ceil((bbM - bbm) / voxelsize).astype(int) + 1

        # Calculate grid centers
        centers = _getGridCenters(bbm, N, voxelsize)
        centers = centers.reshape(np.prod(N), 3)
    else:
        centers = usercenters

    # Calculate features
    if method.upper() == 'C':
        features = _getOccupancyC(mol.coords[:, :, mol.frame], centers, channels)
    elif method.upper() == 'CUDA':
        features = _getOccupancyCUDA(mol.coords[:, :, mol.frame], centers, channels)
    elif method.upper() == 'NUMBA':
        features = _getOccupancyNUMBA(mol.coords[:, :, mol.frame], centers, channels, 5)


    if N is None:
        return features, centers
    else:
        return features, centers, N


def getPointDescriptors(mol, point, size, resolution=1):
    """ Compute descriptors around a specific point in space.

    Parameters
    ----------
    mol:
        A Molecule object. Needs to be read from Autodock 4 .pdbqt format
    point: array-like
        Central box point where descriptors are computed.
    size: array-like
        Size of the resulting box.
    resolution: float
        Resolution of the internal grid.

    Returns
    -------
    features: array-like
        Computed features around queried point.

    """
    size = np.array(size)
    bbm = point - (size * resolution) / 2 + resolution / 2  # Position centers + half res.
    inbox = _getGridCenters(bbm, size, resolution)
    inbox = inbox.reshape(np.prod(size), 3)
    features, _ = getVoxelDescriptors(mol, usercenters=inbox)
    features = features.reshape((size[0], size[1], size[2], features.shape[1]))
    return features


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
    from collections import OrderedDict
    props = OrderedDict()
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
    props['occupancies'] = (elements != 'H') & (elements != 'HS') & (elements != 'HD')

    channels = np.zeros((len(elements), len(props)), dtype=bool)
    for i, p in enumerate(_order):
        channels[:, i] = props[p]
    return channels


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


def _getRadii(mol):
    """ Gets vdW radius for each elem in mol.element. Source VMD.
    
    Parameters
    ----------
    mol :
        A Molecule object. Needs to be read from Autodock 4 .pdbqt format

    Returns
    -------
    radii : np.ndarray
        vdW radius for each element in mol.
    """

    mappings = {  # Mapping pdbqt representation to element.
        'HD': 'H',
        'HS': 'H',
        'A': 'C',
        'NA': 'N',
        'NS': 'N',
        'OA': 'O',
        'OS': 'O',
        'MG': 'Mg',
        'SA': 'S',
        'CL': 'Cl',
        'CA': 'Ca',
        'MN': 'Mn',
        'FE': 'Fe',
        'ZN': 'Zn',
        'BR': 'Br'
    }

    for el in ['H', 'C', 'N', 'O', 'F', 'Mg', 'P', 'S', 'Cl', 'Ca', 'Fe', 'Zn', 'Br', 'I']:
        mappings[el] = el

    res = np.zeros(mol.numAtoms)
    for a in range(mol.numAtoms):
        elem = mol.element[a]

        if elem not in mappings:
            raise ValueError('PDBQT element {} does not exist in mappings.'.format(elem))
        elem = mappings[elem]
        if elem in vdw.radiidict:
            rad = vdw.radiusByElement(elem)
        else:
            print('Unknown element -', mol.element[a], '- at atom index ', a)
            rad = 1.5
        res[a] = rad
    return res


def _getGridCenters(llc, N, resolution):
    xrange = [llc[0] + resolution * x for x in range(0, N[0])]
    yrange = [llc[1] + resolution * x for x in range(0, N[1])]
    zrange = [llc[2] + resolution * x for x in range(0, N[2])]
    centers = np.zeros((N[0], N[1], N[2], 3))
    for i, x in enumerate(xrange):
        for j, y in enumerate(yrange):
            for k, z in enumerate(zrange):
                centers[i, j, k, :] = np.array([x, y, z])
    return centers


def _getOccupancyC(coords, centers, channelsigmas):
    """ Calls the C code to calculate the voxels values for each property."""
    centers = centers.astype(np.float64)
    channelsigmas = channelsigmas.astype(np.float64)

    nchannels = channelsigmas.shape[1]
    occus = np.zeros((centers.shape[0], nchannels))

    occupancylib.descriptor_ext(centers.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                       coords.ctypes.data_as(ctypes.POINTER(ctypes.c_float)),
                       channelsigmas.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                       occus.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                       ctypes.c_int(occus.shape[0]),  # n of centers
                       ctypes.c_int(coords.shape[0]),  # n of atoms
                       ctypes.c_int(nchannels))  # n of channels
    return occus


@jit('float64[:,:](float32[:,:], float64[:,:], float64[:,:], float64)', nopython=True)
def _getOccupancyNUMBA(coords, centers, channelsigmas, trunc):
    ncenters = centers.shape[0]
    natoms = coords.shape[0]
    nchannels = channelsigmas.shape[1]
    trunc = trunc * trunc  # Since we calculate the d**2 we need to get trunc**2
    occus = np.zeros((ncenters, nchannels))
    for a in range(natoms):
        coo = coords[a, :]
        atomsigmas = channelsigmas[a, :]

        for c in range(ncenters):
            cent = centers[c, :]
            d = np.zeros(3)
            d[0] = coo[0] - cent[0]
            d[1] = coo[1] - cent[1]
            d[2] = coo[2] - cent[2]
            d2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2]
            if d2 < trunc:
                d1 = 1 / sqrt(d2)
                for h in range(nchannels):
                    if atomsigmas[h] == 0:
                        continue
                    x = atomsigmas[h] * d1
                    x3 = x * x * x
                    x12 = x3 * x3 * x3 * x3
                    value = 1 - exp(-x12)
                    occus[c, h] = max(occus[c, h], value)
    return occus


def _getOccupancyCUDA(coords, centers, channelsigmas, trunc=5, device=0):
    cuda.select_device(device)
    occus = np.zeros((centers.shape[0], channelsigmas.shape[1]))
    threadsperblock = 1024
    natomblocks = int(np.ceil(coords.shape[0] / threadsperblock))
    blockspergrid = (centers.shape[0], natomblocks)

    centers = cuda.to_device(centers)
    coords = cuda.to_device(coords)
    channelsigmas = cuda.to_device(channelsigmas)
    _getOccupancyCUDAkernel[blockspergrid, threadsperblock](occus, coords, centers, channelsigmas, trunc * trunc)

    return occus

@cuda.jit
def _getOccupancyCUDAkernel(occus, coords, centers, channelsigmas, trunc):
    centeridx = cuda.blockIdx.x
    blockidx = cuda.blockIdx.y
    atomidx = (cuda.threadIdx.x + (cuda.blockDim.x * blockidx))

    if atomidx >= coords.shape[0] or centeridx >= centers.shape[0]:
        return

    # TODO: Can remove this. Barely any speedup
    centcoor = cuda.shared.array(shape=(3), dtype=numba.float32)
    centcoor[0] = centers[centeridx, 0]
    centcoor[1] = centers[centeridx, 1]
    centcoor[2] = centers[centeridx, 2]
    cuda.syncthreads()

    dx = coords[atomidx, 0] - centcoor[0]
    dy = coords[atomidx, 1] - centcoor[1]
    dz = coords[atomidx, 2] - centcoor[2]
    d2 = dx * dx + dy * dy + dz * dz
    if d2 >= trunc:
        return

    d1 = 1 / sqrt(d2)
    for h in range(channelsigmas.shape[1]):
        if channelsigmas[atomidx, h] == 0:
            continue
        x = channelsigmas[atomidx, h] * d1
        value = 1 - exp(-(x ** 12))
        cuda.atomic.max(occus, (centeridx, h), value)




if __name__ == '__main__':
    from htmd.molecule.molecule import Molecule
    from htmd.home import home
    import os
    import numpy as np
    testf = os.path.join(home(), 'data', 'test-voxeldescriptors')
    resOcc, resCent, N = getVoxelDescriptors(Molecule(os.path.join(testf, '3ptb.pdbqt')), buffer=8, voxelsize=1)
    resOcc = resOcc.reshape(N[0], N[1], N[2], resOcc.shape[1])
    refOcc = np.load(os.path.join(testf, '3PTB_occ.npy'))
    refCent = np.load(os.path.join(testf, '3PTB_center.npy'))
    assert np.allclose(resOcc, refOcc)
    assert np.allclose(resCent, refCent)

    import numpy as np
    from htmd.molecule.voxeldescriptors import _getOccupancyC, _getOccupancyCUDA
    centers = np.load(os.path.join(testf, '3PTB_centers_inp.npy'))
    coords = np.load(os.path.join(testf, '3PTB_coords_inp.npy'))
    sigmas = np.load(os.path.join(testf, '3PTB_channels_inp.npy'))
    centers = centers[::10, :].copy()

    res_C = _getOccupancyC(coords, centers, sigmas)
    # res_cuda = _getOccupancyCUDA(coords, centers, sigmas, 5)
    res_numba = _getOccupancyNUMBA(coords, centers, sigmas, 5)

    # assert np.abs(res_C - res_cuda).max() < 1e-4
    assert np.abs(res_C - res_numba).max() < 1e-4






