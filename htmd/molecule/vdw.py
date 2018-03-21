# (c) 2015-2018 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#


"""
All values below are taken from PeriodicTable.C from VMD source code

van der Waals radii are taken from A. Bondi,
J. Phys. Chem., 68, 441 - 452, 1964,
except the value for H, which is taken from R.S. Rowland & R. Taylor,
J.Phys.Chem., 100, 7384 - 7391, 1996. Radii that are not available in
either of these publications have RvdW = 2.00 �.
The radii for Ions (Na, K, Cl, Ca, Mg, and Cs are based on the CHARMM27
Rmin/2 parameters for (SOD, POT, CLA, CAL, MG, CES) by default.
"""

massdict = {'Ac': 227.0, 'Ag': 107.8682, 'Al': 26.981538, 'Am': 243.0, 'Ar': 39.948, 'As': 74.9216, 'At': 210.0,
            'Au': 196.96655, 'B': 10.811, 'Ba': 137.327, 'Be': 9.012182, 'Bh': 264.0, 'Bi': 208.98038, 'Bk': 247.0,
            'Br': 79.904, 'C': 12.0107, 'Ca': 40.078, 'Cd': 112.411, 'Ce': 140.116, 'Cf': 251.0, 'Cl': 35.453,
            'Cm': 247.0, 'Co': 58.9332, 'Cr': 51.9961, 'Cs': 132.90545, 'Cu': 63.546, 'Db': 262.0, 'Ds': 271.0,
            'Dy': 162.5, 'Er': 167.259, 'Es': 252.0, 'Eu': 151.964, 'F': 18.9984032, 'Fe': 55.845, 'Fm': 257.0,
            'Fr': 223.0, 'Ga': 69.723, 'Gd': 157.25, 'Ge': 72.64, 'H': 1.00794, 'He': 4.0026, 'Hf': 178.49,
            'Hg': 200.59, 'Ho': 164.93032, 'Hs': 269.0, 'I': 126.90447, 'In': 114.818, 'Ir': 192.217, 'K': 39.0983,
            'Kr': 83.798, 'La': 138.9055, 'Li': 6.941, 'Lr': 262.0, 'Lu': 174.967, 'Md': 258.0, 'Mg': 24.305,
            'Mn': 54.938049, 'Mo': 95.94, 'Mt': 268.0, 'N': 14.0067, 'Na': 22.98977, 'Nb': 92.90638, 'Nd': 144.24,
            'Ne': 20.1797, 'Ni': 58.6934, 'No': 259.0, 'Np': 237.0, 'O': 15.9994, 'Os': 190.23, 'P': 30.973761,
            'Pa': 231.03588, 'Pb': 207.2, 'Pd': 106.42, 'Pm': 145.0, 'Po': 209.0, 'Pr': 140.90765, 'Pt': 195.078,
            'Pu': 244.0, 'Ra': 226.0, 'Rb': 85.4678, 'Re': 186.207, 'Rf': 261.0, 'Rg': 272.0, 'Rh': 102.9055,
            'Rn': 222.0, 'Ru': 101.07, 'S': 32.065, 'Sb': 121.76, 'Sc': 44.95591, 'Se': 78.96, 'Sg': 266.0,
            'Si': 28.0855, 'Sm': 150.36, 'Sn': 118.71, 'Sr': 87.62, 'Ta': 180.9479, 'Tb': 158.92534, 'Tc': 98.0,
            'Te': 127.6, 'Th': 232.0381, 'Ti': 47.867, 'Tl': 204.3833, 'Tm': 168.93421, 'U': 238.02891, 'V': 50.9415,
            'W': 183.84, 'X': 0.0, 'Xe': 131.293, 'Y': 88.90585, 'Yb': 173.04, 'Zn': 65.409, 'Zr': 91.224}

radiidict = {'Ac': 2.0, 'Ag': 1.72, 'Al': 2.0, 'Am': 2.0, 'Ar': 1.88, 'As': 1.85, 'At': 2.0, 'Au': 1.66, 'B': 2.0,
             'Ba': 2.0, 'Be': 2.0, 'Bh': 2.0, 'Bi': 2.0, 'Bk': 2.0, 'Br': 1.85, 'C': 1.7, 'Ca': 1.37, 'Cd': 1.58,
             'Ce': 2.0, 'Cf': 2.0, 'Cl': 2.27, 'Cm': 2.0, 'Co': 2.0, 'Cr': 2.0, 'Cs': 2.1, 'Cu': 1.4, 'Db': 2.0,
             'Ds': 2.0, 'Dy': 2.0, 'Er': 2.0, 'Es': 2.0, 'Eu': 2.0, 'F': 1.47, 'Fe': 2.0, 'Fm': 2.0, 'Fr': 2.0,
             'Ga': 1.07, 'Gd': 2.0, 'Ge': 2.0, 'H': 1.2, 'He': 1.4, 'Hf': 2.0, 'Hg': 1.55, 'Ho': 2.0, 'Hs': 2.0,
             'I': 1.98, 'In': 1.93, 'Ir': 2.0, 'K': 1.76, 'Kr': 2.02, 'La': 2.0, 'Li': 1.82, 'Lr': 2.0, 'Lu': 2.0,
             'Md': 2.0, 'Mg': 1.18, 'Mn': 2.0, 'Mo': 2.0, 'Mt': 2.0, 'N': 1.55, 'Na': 1.36, 'Nb': 2.0, 'Nd': 2.0,
             'Ne': 1.54, 'Ni': 1.63, 'No': 2.0, 'Np': 2.0, 'O': 1.52, 'Os': 2.0, 'P': 1.8, 'Pa': 2.0, 'Pb': 2.02,
             'Pd': 1.63, 'Pm': 2.0, 'Po': 2.0, 'Pr': 2.0, 'Pt': 1.72, 'Pu': 2.0, 'Ra': 2.0, 'Rb': 2.0, 'Re': 2.0,
             'Rf': 2.0, 'Rg': 2.0, 'Rh': 2.0, 'Rn': 2.0, 'Ru': 2.0, 'S': 1.8, 'Sb': 2.0, 'Sc': 2.0, 'Se': 1.9,
             'Sg': 2.0, 'Si': 2.1, 'Sm': 2.0, 'Sn': 2.17, 'Sr': 2.0, 'Ta': 2.0, 'Tb': 2.0, 'Tc': 2.0, 'Te': 2.06,
             'Th': 2.0, 'Ti': 2.0, 'Tl': 1.96, 'Tm': 2.0, 'U': 1.86, 'V': 2.0, 'W': 2.0, 'X': 1.5, 'Xe': 2.16, 'Y': 2.0,
             'Yb': 2.0, 'Zn': 1.39, 'Zr': 2.0}


def massByElement(element):
    try:
        return massdict[element]
    except:
        return 0.


def radiusByElement(element):
    try:
        return radiidict[element]
    except:
        return 1.5
