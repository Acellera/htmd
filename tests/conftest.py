import os

# Serve Molecule("XXXX") from tests/pdb instead of RCSB to avoid flaky network fetches
os.environ["LOCAL_PDB_REPO"] = os.path.join(os.path.dirname(__file__), "pdb")
