# Minimal makefile for the liblfun documentation.
#
# Sphinx + Breathe over the Doxygen XML of ../include/glfunc.h. This build is
# self-contained under doc/: it does not depend on the library build, and the
# library Makefile does not depend on it. conf.py runs Doxygen itself, so the
# only tools needed are sphinx-build and the doxygen binary on PATH.
#
# Set up the toolchain once, then build from this directory:
#   python3 -m venv .venv
#   .venv/bin/pip install -r requirements.txt
#   sudo apt-get install doxygen        # or: brew install doxygen
#   make html                           # HTML in _build/html/index.html
#
# Override the tool or options from the command line, e.g.
#   make html SPHINXBUILD=/path/to/sphinx-build
#   make html SPHINXOPTS=               # do not treat warnings as errors

# Prefer a sphinx-build under .venv if a contributor created one here; otherwise
# use whatever is on PATH (an activated venv, a system package, Read the Docs).
SPHINXBUILD  ?= $(if $(wildcard .venv/bin/sphinx-build),.venv/bin/sphinx-build,sphinx-build)
# -W --keep-going turns Sphinx warnings into errors (collecting them all first),
# so a broken cross-reference or an orphaned page fails the build.
SPHINXOPTS   ?= -W --keep-going
SOURCEDIR     = .
BUILDDIR      = _build

# `make` with no target prints Sphinx's list of builders (html, clean, ...).
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help Makefile

# Catch-all: route any target (html, clean, linkcheck, ...) through sphinx-build.
%: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)
