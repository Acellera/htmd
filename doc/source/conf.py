# Configuration file for the Sphinx documentation builder.
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
from pathlib import Path

from acellera_docs_theme import apply

# Drop the timestamp from HTMD / moleculekit log lines in rendered tutorial output.
os.environ.setdefault("HTMD_LOG_FORMAT", "%(name)s - %(levelname)s - %(message)s")
os.environ.setdefault("MOLECULEKIT_LOG_FORMAT", "%(name)s - %(levelname)s - %(message)s")

# -- Project information -----------------------------------------------------

project = "HTMD"
copyright = "2026, Acellera"
author = "Acellera"

# -- General configuration ---------------------------------------------------

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.doctest",
    "sphinx.ext.intersphinx",
    "sphinx.ext.mathjax",
    "sphinx.ext.ifconfig",
    "sphinx.ext.napoleon",
    "sphinx_autodoc_typehints",
    "sphinx_design",
    "sphinx_copybutton",
    "sphinxcontrib.mermaid",
    "myst_nb",
    "sphinxarg.ext",
]

source_suffix = {
    ".rst": "restructuredtext",
    ".md": "myst-nb",
    ".ipynb": "myst-nb",
}

exclude_patterns = ["build", "**/.ipynb_checkpoints"]

# -- MyST / MyST-NB ----------------------------------------------------------

myst_enable_extensions = [
    "colon_fence",
    "deflist",
    "substitution",
    "tasklist",
]
myst_heading_anchors = 3

nb_execution_mode = "auto"  # execute .md tutorials on first build; .ipynb tutorials ship with cached outputs
nb_execution_excludepatterns = ["**/*.ipynb"]  # legacy notebooks carry their own outputs (and pop viewers if rerun)
nb_execution_in_temp = True  # run each notebook in a temp dir so build artefacts don't pollute source/
nb_execution_timeout = 600  # antechamber + tLeap inside the system-prep tutorials can take minutes
nb_merge_streams = True

# -- Acellera unified branding ----------------------------------------------

apply(
    globals(),
    project_name="HTMD",
    github_repo="Acellera/htmd",
)

# -- LLM full-corpus artifact ------------------------------------------------


def _emit_llms_full_txt(app, exception):
    """build-finished hook: concatenate every rendered page source into llms-full.txt."""
    if exception is not None:
        return
    srcdir = Path(app.srcdir)
    output = srcdir / "llms-full.txt"
    parts = []
    for path in sorted(srcdir.rglob("*.md")):
        if "build" in path.parts:
            continue
        rel = path.relative_to(srcdir)
        parts.append(f"# === {rel} ===\n\n{path.read_text(encoding='utf-8')}\n")
    for path in sorted(srcdir.rglob("*.rst")):
        if "build" in path.parts:
            continue
        rel = path.relative_to(srcdir)
        parts.append(f"# === {rel} ===\n\n{path.read_text(encoding='utf-8')}\n")
    output.write_text("\n".join(parts), encoding="utf-8")


def setup(app):
    app.connect("build-finished", _emit_llms_full_txt)
    return {"version": "1.0", "parallel_read_safe": True}
