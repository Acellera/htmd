from htmd.util import tempname
from subprocess import call
import numpy as np
from htmd.molecule.molecule import Molecule


def loopModeller(mol, segid, seq, startresid, movstart=None, movend=None, modellerexe='/shared/sdoerr/Software/modeller-9.16/build/bin/mod9.16'):
    """ Uses the Modeller software to predict missing loops in a Molecule.

    Parameters
    ----------
    mol : :class:`Molecule <htmd.molecule.molecule.Molecule>` object
        A Molecule object.
    segid : str
        The name of the segment containing the gap.
    seq : str
        The sequence of residues to be added by the loop modeller.
    startresid : int
        The resid of the residue before the gap.
    movstart : int
        The resid after which the residues will be allowed to move in modeller.
    movend : int
        The resid before which the residues will be allowed to move in modeller.
    modellerexe : str
        The path to the modeller executable.

    Returns
    -------
    newmol : :class:`Molecule <htmd.molecule.molecule.Molecule>` object
        A new Molecule object containing the protein with the modelled loop.

    Examples
    --------
    >>> newmol = loopModeller(mol, 'B', 'MLSGSK', 263)
    """
    if movstart is None:
        movstart = startresid
    if movend is None:
        movend = movstart + len(seq) + 1

    segatm = mol.atomselect('segid {}'.format(segid))
    segres = np.unique(mol.get('resid', sel=segatm))
    chain = np.unique(mol.get('chain', sel=segatm))
    if len(chain) != 1:
        raise RuntimeError('More than one chain detected in segment {}'.format(segid))
    pos = np.where(segres == startresid)[0][0]
    if (segres[pos+1] - segres[pos]) != (len(seq) + 1):
        raise RuntimeError('Sequence is {} characters long while sequence gap ({}-{}) is {} long. Cannot insert loop.'.format(
            len(seq), segres[pos], segres[pos+1], (segres[pos+1] - segres[pos])))

    chain = chain[0]
    segresidstart = np.min(segres)
    segresidend = np.max(segres)

    currseq = mol.sequence()[segid]
    minuses = ['-'] * len(seq)
    minuses = ''.join(minuses)
    inspos = np.where(segres == startresid)[0] + 1
    gapseq = currseq[:inspos] + minuses + currseq[inspos:] + '*'
    fullseq = currseq[:inspos] + seq + currseq[inspos:] + '*'

    # Get the sequence of the 1qg8 PDB file, and write to an alignment file
    pdbfile = tempname(suffix='.pdb')
    mol.write(pdbfile)

    alifile = tempname(suffix='.ali')
    f = open(alifile, 'w')
    f.write('>P1;prot\nstructure:{}:{}:{}:{}:{}::::\n'.format(pdbfile, segresidstart, chain, segresidend, chain))
    f.write(gapseq)
    f.write('\n>P1;prot_fill\nsequence:::::::::\n')
    f.write(fullseq)
    f.close()

    script = '''
from modeller import *
from modeller.automodel import *    # Load the automodel class

log.verbose()
env = environ()

class MyModel(automodel):
    def select_atoms(self):
        return selection(self.residue_range('{}', '{}'))

a = MyModel(env, alnfile = '{}', knowns = 'prot', sequence = 'prot_fill')
a.starting_model= 1
a.ending_model  = 1

a.make()
    '''.format(movstart+1-segresidstart, movend-segresidstart, alifile)
    pyfile = tempname(suffix='.py')
    f = open(pyfile, 'w')
    f.write(script)
    f.close()

    print(pdbfile)
    print(alifile)
    print(pyfile)

    call('{} {}'.format(modellerexe, pyfile), shell=True)
    newmol = Molecule('./prot_fill.B99990001.pdb')
    return newmol


'''
def loopmodelerFALC(mol, segid, seq, pos, inspos, outname='loop'):
    currseq = mol.sequence()[segid]
    #from IPython.core.debugger import Tracer
    #Tracer()()
    fullseq = currseq[:inspos] + seq + currseq[inspos:]

    seqfile = tempname(suffix='.fasta')
    f = open(seqfile, 'w')
    f.write(fullseq)
    f.close()

    ulrfile = tempname(suffix='.ulr')  # Unreliable local region file: http://galaxy.seoklab.org/softwares/falc.html
    f = open(ulrfile, 'w')
    f.write('1 L {}-{} 0.000'.format(pos, pos+len(seq)+1))
    f.close()

    pdbfile = tempname(suffix='.pdb')
    mol.write(pdbfile)

    call('/shared/sdoerr/Software/FALC/bin/falc.py -title {} -pdb {} -fa {} -ulr {} -n {}'.format(outname, pdbfile, seqfile, ulrfile, 1))
'''