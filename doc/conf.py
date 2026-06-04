# Sphinx configuration for the liblfun documentation.
#
# TOOLCHAIN CHOICE (recorded here per bead lfunctions-ofr)
# -------------------------------------------------------
# Doxygen + Breathe + Sphinx (with MyST for Markdown sources).
#
# Why this stack:
#   * include/glfunc.h is the public API and is already commented, so we keep it
#     as the single source of truth for the per-function reference. Doxygen
#     parses the header into XML; Breathe surfaces that XML through Sphinx's C
#     domain. Hand-written narrative (MyST Markdown under doc/*.md) sits on top
#     and links into the generated reference.
#   * Doxygen is packaged everywhere we build (apt, Read the Docs native
#     support, GitHub Actions) and emits structured C-domain XML, which is more
#     robust than scraping a header with a libclang Python binding whose version
#     must track the system libclang (the hawkmoth / sphinx-c-autodoc route).
#
# Doxygen runs automatically below (run_doxygen), so a plain `sphinx-build` and a
# Read the Docs build both regenerate the XML; the doc/Makefile build only has to
# invoke sphinx-build. The Doxygen settings live in doc/Doxyfile.
#
# NOTE for the API-guide author (bead lfunctions-rrt): the reference is wired up
# in doc/api.md via Breathe's `doxygenfile` directive against the project named
# "liblfun" below. Fill that page; no change is needed here.

import os
import shutil
import subprocess

# -- Paths -------------------------------------------------------------------
# conf.py lives in doc/; the repository root is its parent. Doxygen is run from
# the repo root so its relative INPUT/OUTPUT paths match the Makefile target.
_here = os.path.abspath(os.path.dirname(__file__))
_repo_root = os.path.abspath(os.path.join(_here, ".."))

# -- Project information -----------------------------------------------------
project = "liblfun"
author = "Bober, Booker, Costa, Lee, Platt, and Sutherland"
copyright = "2016-2026 Edgar Costa and contributors"

# Version is intentionally light: liblfun ships no public version macro yet, so
# track the Autoconf package version (AC_INIT in configure.ac).
release = "0.1"
version = "0.1"

# -- General configuration ---------------------------------------------------
extensions = [
    "myst_parser",       # Markdown (MyST) sources, so doc/*.md just work
    "breathe",           # bridge Doxygen XML -> Sphinx C domain
    "sphinx.ext.intersphinx",
]

# Accept both Markdown and reStructuredText sources.
source_suffix = {
    ".md": "markdown",
    ".rst": "restructuredtext",
}

master_doc = "index"
root_doc = "index"

exclude_patterns = ["_build", "Thumbs.db", ".DS_Store", ".venv"]

# Treat the C header as C (not C++) when Sphinx resolves cross-references.
primary_domain = "c"
highlight_language = "c"

# MyST niceties used by the narrative pages (and by sibling-PR topic pages).
myst_enable_extensions = [
    "colon_fence",       # ::: fenced directives, friendlier than ```{...}
    "deflist",
    "linkify",           # bare URLs become links
    "smartquotes",
]
myst_heading_anchors = 3  # auto heading anchors so pages can deep-link

# -- Breathe (Doxygen bridge) ------------------------------------------------
# XML is generated under doc/_build/doxygen/xml by run_doxygen() below.
_doxygen_xml = os.path.join(_here, "_build", "doxygen", "xml")
breathe_projects = {"liblfun": _doxygen_xml}
breathe_default_project = "liblfun"
breathe_domain_by_extension = {"h": "c"}
breathe_show_define_initializer = True  # show the values of the ERR_* / tri-state macros


def run_doxygen(app=None):
    """Generate the Doxygen XML that Breathe consumes.

    Run from the repository root so doc/Doxyfile's relative paths resolve the
    same way regardless of where sphinx-build is invoked from (doc/Makefile, or
    Read the Docs). If the
    doxygen binary is missing we emit a clear, actionable message instead of
    failing with a cryptic Breathe "project ... not found" error later.
    """
    if shutil.which("doxygen") is None:
        msg = (
            "doxygen not found on PATH: the API reference cannot be generated. "
            "Install it (Debian/Ubuntu: `apt-get install doxygen`; macOS: "
            "`brew install doxygen`) and rebuild."
        )
        if app is not None:
            from sphinx.util import logging as sphinx_logging

            sphinx_logging.getLogger(__name__).warning(msg)
        else:
            print("WARNING:", msg)
        return
    subprocess.run(
        ["doxygen", os.path.join("doc", "Doxyfile")],
        cwd=_repo_root,
        check=True,
    )


def setup(app):
    # Regenerate the XML at the start of every build (local and Read the Docs).
    app.connect("builder-inited", run_doxygen)
    return {"parallel_read_safe": True, "parallel_write_safe": True}


# Read the Docs builds run conf.py directly (no Makefile), so make sure the XML
# exists before the build proper begins there too.
if os.environ.get("READTHEDOCS") == "True":
    run_doxygen()

# -- HTML output -------------------------------------------------------------
html_theme = "furo"
html_title = "liblfun documentation"
html_static_path = ["_static"]
html_css_files = ["custom.css"]

# -- Cross-project references ------------------------------------------------
intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
}
