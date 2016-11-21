from htmd import *
from htmd.builder.builder import sequenceID
from htmd.molecule.molecule import _residueNameTable
import ctypes as ct

mylib = ct.cdll.LoadLibrary("/shared/sdoerr/Work/htmdacellera/htmd/lib/Linux/tmalign.so")
tmalign = mylib.tmalign

ref = Molecule('./data/ntl9_2hbb.pdb')
mol = Molecule('./data/filtered.pdb')
# Taken from /workspace3/casp/folding/crystal/ss/1ns/ntl9/batches/0/filtered/e30s1_e28s7p0f1/e30s1_e28s7p0f1-SDOERR_CASP22S_crystal_ss_1ns_ntl9_0-0-1-RND1648_9.filtered.xtc
mol.read('./data/traj.xtc')

#ref = Molecule('PDB2.pdb')
#mol = Molecule('PDB1.pdb')

# def inefficientPP(array2D, ct_type):
#     # It's more efficient to pass single pointers and index in C. However this is quite convenient.
#     return (ct.POINTER(ct_type) * array2D.shape[0])(*(r.ctypes.data_as(ct.POINTER(ct_type)) for r in array2D))

def calculateVariables(mol):
    ASCIIlimit = 123
    res = sequenceID((mol.resid, mol.insertion, mol.segid, mol.chain))
    mol.segid[:] = 'X'
    mol.chain[:] = 'X'
    caidx = mol.name == 'CA'

    cares = res[caidx]
    res = np.unique(res)
    reslen  = len(res)

    res = res.astype(np.int32).ctypes.data_as(ct.POINTER(ct.c_int))

    # Calculate the protein sequence
    seq = ''.join([_residueNameTable[x] for x in mol.resname[caidx]])
    seq = ct.c_char_p(seq.encode('utf-8'))

    # Keep only CA coordinates
    coords = mol.coords[caidx, :, :].copy()
    return reslen, res, seq, coords

def callTMAlign(ref, mol):
    mol = mol.copy()
    ref = ref.copy()
    ref.filter('protein', _logger=False)
    mol.filter('protein', _logger=False)

    reslenMOL, residMOL, seqMOL, coordsMOL = calculateVariables(mol)
    reslenREF, residREF, seqREF, coordsREF = calculateVariables(ref)

    # void tmalign(int xlen, int ylen, int* xresno, int* yresno, char* seqx, char* seqy, float* xcoor, float* ycoor, int nframes, double &TM1, double &TM2, double &rmsd)
    resTM1 = (ct.c_double * mol.numFrames)()
    resTM2 = (ct.c_double * mol.numFrames)()
    resRMSD = (ct.c_double * mol.numFrames)()
    retval = fn = tmalign(ct.c_int(reslenREF), 
                          ct.c_int(reslenMOL), 
                          residREF,
                          residMOL,
                          seqREF,
                          seqMOL,
                          coordsREF.ctypes.data_as(ct.POINTER(ct.c_float)),
                          coordsMOL.ctypes.data_as(ct.POINTER(ct.c_float)),
                          ct.c_int(mol.numFrames),
                          ct.byref(resTM1),
                          ct.byref(resTM2),
                          ct.byref(resRMSD))
    return np.ctypeslib.as_array(resTM1), np.ctypeslib.as_array(resTM2), np.ctypeslib.as_array(resRMSD)
    
import time
t = time.time()
resTM1, resTM2, resRMSD = callTMAlign(ref, mol)
print('time: {}'.format(time.time() - t))
print(resTM1, resTM2, resRMSD)


