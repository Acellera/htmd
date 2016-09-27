ff-nucleic-OL15 force field is ff99bsc0 plus the following Olomouc (OL) dihedral corrections:

DNA:
(i) chiOL4 - correction for DNA glycosidic torsion chi
(ii) ezOL1 - correction for eps/zeta torsions
(iii) bOL1 - correction for beta torsion

RNA:
(iv) chiOL3 - correction for RNA glycosidic torsion chi


Parameters may be read into tLEaP using command:
tleap -s -f leap-ff-nucleic-OL15.in

Files are prepared for use with AMBER 14. Note that all downloaded files should be places in the current directory (from which tLEaP was launched) or you should define the path to these in leap-ff-nucleic-OL15.in

The following atom types were changed in DNA nucleotides (no changes in RNA nucleotides):
atom name C5' is changed from CI to CJ type in all nucleotides
atom name C3' is changed from CT to C7 type in all nucleotides
atom name C8 is changed from CK to C2 type in Adenines
atom names C5 and C6 are changed from CM to C1 type in Cytosines
(For RNA the atom types are the same as in ff10, ff12SB and ff14SB)



Please cite these corrections as:

(i) Krepl, M.; Zgarbova, M.; Stadlbauer, P.; Otyepka, M.; Banas, P.; Koca, J.; Cheatham, T. E.; Jurecka, P.; Sponer, J., Reference Simulations of Noncanonical Nucleic Acids with Different chi Variants of the AMBER Force Field: Quadruplex DNA, Quadruplex RNA, and Z-DNA. J. Chem. Theory Comput. 2012, 8 (7), 2506-2520. 
chiOL3 (RNA):

(ii) Zgarbova, M.; Luque, F. J.; Sponer, J.; Cheatham, T. E.; Otyepka, M.; Jurecka, P., Toward Improved Description of DNA Backbone: Revisiting Epsilon and Zeta Torsion Force Field Parameters. J. Chem. Theory Comput. 2013, 9 (5), 2339-2354.

(iii) Zgarbova, M.; Spner, J.; Otyepka, M., Cheatham, T. E.; Galindo-Murillo, R.; Jurecka, P., Refinement of the Sugar-Phosphate Backbone Torsion Beta for AMBER Force Fields Improves the Description of Z- and B-DNA. J. Chem. Theory Comput. 2015, DOI:10.1021/acs.jctc.5b00716. 

(iv) Zgarbova, M.; Otyepka, M.; Sponer, J.; Mladek, A.; Banas, P.; Cheatham, T. E.; Jurecka, P., Refinement of the Cornell et al. Nucleic Acids Force Field Based on Reference Quantum Chemical Calculations of Glycosidic Torsion Profiles. J. Chem. Theory Comput. 2011, 7 (9), 2886-2902.
Banas, P.; Hollas, D.; Zgarbova, M.; Jurecka, P.; Orozco, M.; Cheatham, T. E.; Sponer, J.; Otyepka, M., Performance of Molecular Mechanics Force Fields for RNA Simulations: Stability of UUCG and GNRA Hairpins. J. Chem. Theory Comput. 2010, 6 (12), 3836-3849.

If you have any questions please contact us (Petr Jurecka: petr.jurecka@upol.cz; Marie Zgarbova: marie.zgarbova@upol.cz)
