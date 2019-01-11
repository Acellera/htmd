# (c) 2015-2018 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import ctypes
import os
import unittest
from math import exp, sqrt
from unittest import TestCase

import htmd.home
import numba
import numpy as np
from htmd.home import home
from htmd.molecule import vdw
from htmd.molecule.atomtyper import getFeatures, molPDBQT, prepProtForFeats
from htmd.molecule.molecule import Molecule
from htmd.molecule.util import boundingBox
from numba import cuda, jit

_order = ('hydrophobic', 'aromatic', 'hbond_acceptor', 'hbond_donor', 'positive_ionizable',
          'negative_ionizable', 'metal', 'occupancies')
libdir = htmd.home.home(libDir=True)
occupancylib = ctypes.cdll.LoadLibrary(
    os.path.join(libdir, "occupancy_ext.so"))


def getVoxelDescriptors(mol, usercenters=None, voxelsize=1, buffer=0, channels=None, aromaticNitrogen=False, method='C'):
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
    mol = molPDBQT(prepProtForFeats(mol), aromaticNitrogen=aromaticNitrogen)

    if channels is None:
        channels = getFeatures(mol)

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
        features = _getOccupancyC(
            mol.coords[:, :, mol.frame], centers, channels)
    elif method.upper() == 'CUDA':
        features = _getOccupancyCUDA(
            mol.coords[:, :, mol.frame], centers, channels)
    elif method.upper() == 'NUMBA':
        features = _getOccupancyNUMBA(
            mol.coords[:, :, mol.frame], centers, channels, 5)

    if N is None:
        return features, centers
    else:
        return features, centers, N


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

    for el in vdw.radiidict.keys():
        mappings[el] = el

    res = np.zeros(mol.numAtoms)
    for a in range(mol.numAtoms):
        elem = mol.element[a]

        if elem not in mappings:
            raise ValueError(
                'PDBQT element {} does not exist in mappings.'.format(elem))
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
                                coords.ctypes.data_as(
                                    ctypes.POINTER(ctypes.c_float)),
                                channelsigmas.ctypes.data_as(
                                    ctypes.POINTER(ctypes.c_double)),
                                occus.ctypes.data_as(
                                    ctypes.POINTER(ctypes.c_double)),
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


def _memsetArray(array, val=0, threadsperblock=256):
    from math import ceil
    totalelem = np.prod(array.shape)
    nblocks = ceil(totalelem / threadsperblock)
    _memsetArrayCUDAkernel[nblocks, threadsperblock](
        array.reshape(totalelem), val)


@cuda.jit
def _memsetArrayCUDAkernel(array, val):
    threadidx = (cuda.threadIdx.x + (cuda.blockDim.x * cuda.blockIdx.x))
    if threadidx >= array.shape[0]:
        return
    array[threadidx] = val


def _getOccupancyCUDA(coords, centers, channelsigmas, trunc=5, device=0, resD=None, asnumpy=True, threadsperblock=256):
    # cuda.select_device(device)
    if resD is None:
        resD = cuda.device_array(
            (centers.shape[0], channelsigmas.shape[1]), dtype=np.float32)
    _memsetArray(resD, val=0)

    natomblocks = int(np.ceil(coords.shape[0] / threadsperblock))
    blockspergrid = (centers.shape[0], natomblocks)

    centers = cuda.to_device(centers)
    coords = cuda.to_device(coords)
    channelsigmas = cuda.to_device(channelsigmas)
    _getOccupancyCUDAkernel[blockspergrid, threadsperblock](
        resD, coords, centers, channelsigmas, trunc * trunc)

    if asnumpy:
        return resD.copy_to_host()


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


class TestVoxel(TestCase):
    @classmethod
    def setUpClass(self):
        self.testf = os.path.join(home(), 'data', 'test-voxeldescriptors')
        self.mol = Molecule(os.path.join(self.testf, '3PTB.pdb'))
        self.refOcc = np.load(os.path.join(self.testf, '3PTB_occ.npy'))
        self.refCent = np.load(os.path.join(self.testf, '3PTB_centers.npy'))

    def test_featC(self):
        resOcc, resCent, N = getVoxelDescriptors(self.mol, method='C')
        assert np.allclose(resOcc, self.refOcc)
        assert np.allclose(resCent, self.refCent)

    def test_featNUMBA(self):
        resOcc, resCent, N = getVoxelDescriptors(self.mol, method='NUMBA')
        assert np.allclose(resOcc, self.refOcc)
        assert np.allclose(resCent, self.refCent)

    # def test_featCUDA(self):
    #     resOcc, resCent, N = getVoxelDescriptors(self.mol, method='CUDA')
    #     assert np.allclose(resOcc, self.refOcc)
    #     assert np.allclose(resCent, self.refCent)


if __name__ == '__main__':
    unittest.main(verbosity=2)
