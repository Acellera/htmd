# Configuration file for the Sphinx documentation builder.
# https://www.sphinx-doc.org/en/master/usage/configuration.html

from pathlib import Path

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
}

templates_path = ["_templates"]
exclude_patterns = ["build", "**/.ipynb_checkpoints"]

# -- MyST / MyST-NB ----------------------------------------------------------

myst_enable_extensions = [
    "colon_fence",
    "deflist",
    "substitution",
    "tasklist",
]
myst_heading_anchors = 3

nb_execution_mode = "off"  # tutorials are pre-rendered; flip to "cache" when adopting MyST-NB
nb_merge_streams = True

# -- HTML output -------------------------------------------------------------

html_theme = "pydata_sphinx_theme"
html_static_path = ["_static"]
html_css_files = ["custom.css"]
html_show_sourcelink = True
html_logo = "_static/img/acellera_new_web.png"
html_favicon = "_static/img/acellera-logo-16x16.png"
html_context = {"default_mode": "light"}

html_theme_options = {
    "logo": {
        "image_light": "_static/img/acellera_new_web.png",
        "image_dark": "_static/img/acellera_new_web.png",
        "text": "HTMD",
    },
    "header_links_before_dropdown": 5,
    "show_toc_level": 2,
    "navigation_depth": 3,
    "use_edit_page_button": False,
    "navigation_with_keys": False,
    "footer_start": ["copyright"],
    "footer_end": [],
    "icon_links": [
        {
            "name": "Acellera",
            "url": "https://www.acellera.com",
            "icon": "_static/img/acellera-logo-white.png",
            "type": "local",
        },
        {
            "name": "Twitter",
            "url": "https://twitter.com/acellera",
            "icon": "fab fa-twitter",
            "type": "fontawesome",
        },
        {
            "name": "GitHub",
            "url": "https://github.com/Acellera/htmd",
            "icon": "fab fa-github-square",
            "type": "fontawesome",
        },
        {
            "name": "LinkedIn",
            "url": "https://www.linkedin.com/company/acellera/",
            "icon": "fab fa-linkedin",
            "type": "fontawesome",
        },
        {
            "name": "Youtube",
            "url": "https://www.youtube.com/user/acelleralive",
            "icon": "fab fa-youtube",
            "type": "fontawesome",
        },
    ],
}

html_sidebars = {
    "**": ["sidebar-nav-bs.html"],
}

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
