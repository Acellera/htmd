# Small Molecule (Sub)Package

Package for efficient handling of small molecules.


## Voxelization 

example:

```python
from htmd.smallmol.smallmol import SmallMol, SmallMolStack
ms = SmallMolStack("/path/to/my/file/x.sdf")  # Generate Molecule stack 
my_gen = ms.voxel_generator(batch_size=32, n_jobs=10)  # Make a Voxel generator
for voxellist in my_gen:  # Iterate over the elements
    pass  # Do something with it.
```

## Tautomer generation
 
 
Generate and filter tautomers:
```python
from htmd.smallmol.smallmol import SmallMol
from htmd.smallmol.tautomers import TautomerCanonicalizer

# input smile string / .mol2 file / rdkit molecule
my_mol = SmallMol("/path/to/some/file/mf.mol2")

canon = TautomerCanonicalizer()
tautomers, scores, depictatoms, details = canon.compute_score(my_mol, returndetails=True)

# Return filtered list of SmallMol objects
tautomers_mols = canon.filter_tautomers(tautomers, scores)
```

Visalize tauotomers:
```python
from htmd.smallmol.smallmol import SmallMol
from htmd.smallmol.tautomers import TautomerCanonicalizer, MolFragmenter
from rdkit.Chem.Draw import rdMolDraw2D
from IPython.display import SVG
drawer = rdMolDraw2D.MolDraw2DSVG(400,200)

my_mol = SmallMol("/path/to/some/file/mf.mol2")
mf = MolFragmenter(my_mol)
mf.fragment()
mf.compute_paths()
path = mf.get_paths[0]
atoms = [a for fr in path for a in fr.atoms]

drawer.DrawMolecule(my_mol._mol, highlightAtoms=atoms)
drawer.FinishDrawing()
svg = drawer.GetDrawingText().replace('svg:','')
img = SVG(svg)
```
