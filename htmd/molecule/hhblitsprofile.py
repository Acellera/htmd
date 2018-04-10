# (c) 2015-2018 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import tempfile
import os
import re
import numpy as np
import pandas


def getProfile(sequence, hhblits, hhblitsdb, ncpu=6, niter=4):
    """ Calculates the sequence profile of a protein sequence using HHBlits

    File description and unit conversions taken from section 6 https://hpc.nih.gov/apps/hhsuite-userguide.pdf

    Parameters
    ----------
    sequence : str
        A string encoding the one letter sequence of the protein
    hhblits : str
        The path to the hhblits executable
    hhblitsdb : str
        The path to the hhblits database that we want to search against
    ncpu : int
        Number of CPUs to use for the search
    niter : int
        The number of hhblits iterations. The higher the value the more remote homologues it will find
    
    Returns
    -------
    df : pandas.DataFrame
        A pandas dataframe containing all the information read from the file
    pssm : np.ndarray
        A Nx20 numpy array where N the number of residues of the protein. Contains the transition probabilities to all 20 residues.

    Examples
    --------
    >>> hhb = '~/hhsuite-2.0.16-linux-x86_64/bin/hhblits'
    >>> hhbdb = '~/hhsuite-2.0.16-linux-x86_64/databases/uniprot20_2016_02/uniprot20_2016_02'
    >>> seq = 'MKVIFLKDVKGMGKKGEIKNVADGYANNFLFKQGLAIEA'
    >>> df, prof = getProfile(seq, hhb, hhbdb)
    """
    regex = re.compile('^\w\s\d+')
    with tempfile.TemporaryDirectory() as tmpdir:
        with open(tmpdir + '/input.fasta', 'w') as fp:
            fp.write('>\n{}\n'.format(sequence))

        os.system('{} -i {}/input.fasta -d {} -cpu {} -M first -n {} -ohhm {}/output.hhm -v 0'.format(hhblits, tmpdir, hhblitsdb, ncpu, niter, tmpdir))
        
        data = []
        seq = []
        with open(tmpdir + '/output.hhm', 'r') as fp:
            starting = 0
            lines = fp.readlines()
            for i in range(len(lines)):
                if lines[i].startswith('NULL'):
                    pieces = lines[i].split()
                    seq.append([pieces[0]])
                    data.append([2 ** (-int(x) / 1000) if x != '*' else 0 for x in pieces[1:21]] + [0] * 10)
                if lines[i].startswith('HMM    A	C	D'):
                    col_desc = lines[i].split()[1:] + lines[i + 1].split()
                    starting = 1
                if starting > 0:
                    starting += 1
                if starting >= 4 and regex.match(lines[i]):
                    pieces = lines[i].split()
                    seq.append([pieces[0]])
                    d = [2 ** (-int(x) / 1000) if x != '*' else 0 for x in pieces[2:22]]
                    pieces = lines[i + 1].split()
                    d += [2 ** (-int(x) / 1000) if x != '*' else 0 for x in pieces[:7]]
                    d += [0.001 * int(x) for x in pieces[7:10]]
                    data.append(d)

        df = pandas.DataFrame(np.hstack((np.vstack(seq), np.vstack(data))), columns=['seq'] + col_desc)
        return df, df.as_matrix()[1:, 1:21].astype(float)


if __name__ == '__main__':
    pass
