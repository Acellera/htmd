from htmd import *

cases="""
1z98
2w2e
2z73
3cll
3jbr
3lut
3spc
3spg
3syq
4csk
4hkr
4kfm
4pe5
4tlm
5an8
""".split()

for p in cases:
    m=Molecule(p)
    mo,rd=prepareProtein(m)
