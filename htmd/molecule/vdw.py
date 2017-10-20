# (c) 2015-2017 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
# All values taken from VMD
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
elements = [
    "X", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
    "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc",
    "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge",
    "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc",
    "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe",
    "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb",
    "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os",
    "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr",
    "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf",
    "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt",
    "Ds", "Rg"
]  # periodic table of elements for translation of ordinal to atom type

masses = [
    0.00000, 1.00794, 4.00260, 6.941, 9.012182, 10.811,
    12.0107, 14.0067, 15.9994, 18.9984032, 20.1797,
    22.989770, 24.3050, 26.981538, 28.0855, 30.973761,
    32.065, 35.453, 39.948, 39.0983, 40.078, 44.955910,
    47.867, 50.9415, 51.9961, 54.938049, 55.845, 58.9332,
    58.6934, 63.546, 65.409, 69.723, 72.64, 74.92160,
    78.96, 79.904, 83.798, 85.4678, 87.62, 88.90585,
    91.224, 92.90638, 95.94, 98.0, 101.07, 102.90550,
    106.42, 107.8682, 112.411, 114.818, 118.710, 121.760,
    127.60, 126.90447, 131.293, 132.90545, 137.327,
    138.9055, 140.116, 140.90765, 144.24, 145.0, 150.36,
    151.964, 157.25, 158.92534, 162.500, 164.93032,
    167.259, 168.93421, 173.04, 174.967, 178.49, 180.9479,
    183.84, 186.207, 190.23, 192.217, 195.078, 196.96655,
    200.59, 204.3833, 207.2, 208.98038, 209.0, 210.0, 222.0,
    223.0, 226.0, 227.0, 232.0381, 231.03588, 238.02891,
    237.0, 244.0, 243.0, 247.0, 247.0, 251.0, 252.0, 257.0,
    258.0, 259.0, 262.0, 261.0, 262.0, 266.0, 264.0, 269.0,
    268.0, 271.0, 272.0
]  # corresponding table of masses

radii = [
    1.50, 1.20, 1.40, 1.82, 2.00, 2.00,
    1.70, 1.55, 1.52, 1.47, 1.54,
    1.36, 1.18, 2.00, 2.10, 1.80,
    1.80, 2.27, 1.88, 1.76, 1.37, 2.00,
    2.00, 2.00, 2.00, 2.00, 2.00, 2.00,
    1.63, 1.40, 1.39, 1.07, 2.00, 1.85,
    1.90, 1.85, 2.02, 2.00, 2.00, 2.00,
    2.00, 2.00, 2.00, 2.00, 2.00, 2.00,
    1.63, 1.72, 1.58, 1.93, 2.17, 2.00,
    2.06, 1.98, 2.16, 2.10, 2.00,
    2.00, 2.00, 2.00, 2.00, 2.00, 2.00,
    2.00, 2.00, 2.00, 2.00, 2.00,
    2.00, 2.00, 2.00, 2.00, 2.00, 2.00,
    2.00, 2.00, 2.00, 2.00, 1.72, 1.66,
    1.55, 1.96, 2.02, 2.00, 2.00, 2.00, 2.00,
    2.00, 2.00, 2.00, 2.00, 2.00, 1.86,
    2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00,
    2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00,
    2.00, 2.00, 2.00
]  # corresponding table of VDW radii


def massByElement(element):
    try:
        return masses[elements.index(element)]
    except:
        return 0.


def radiusByElement(element):
    try:
        return radii[elements.index(element)]
    except:
        return 1.5
