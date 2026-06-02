# Makefile for Sphinx documentation
#

SPHINXOPTS    ?=
SPHINXBUILD   ?= sphinx-build
SPHINXAPIDOC  ?= sphinx-apidoc
PAPER         =
BUILDDIR      = build
SOURCEDIR     = source
TUTORIALSDIR  = $(SOURCEDIR)/tutorials
HOWTODIR      = $(SOURCEDIR)/how-to
REFERENCEDIR  = $(SOURCEDIR)/reference
HTMDTAG       = $(shell python gethtmdtag.py tag)
HTMDBRANCH    = $(shell python gethtmdtag.py branch)

# Notebook mapping: each list is the basename set that lands in a specific
# Diátaxis subdir. Sources all live in ../tutorials/ at repo root (so the
# tutorials test fixture stays put).
TUT_ADAPTIVE  = adaptive-sampling.ipynb adaptive-bandit.ipynb
TUT_ANALYSIS  =
HT_NOTEBOOKS  =

# User-friendly check for sphinx-build
ifeq ($(shell which $(SPHINXBUILD) >/dev/null 2>&1; echo $$?), 1)
$(error The '$(SPHINXBUILD)' command was not found. Make sure you have Sphinx installed, then set the SPHINXBUILD environment variable to point to the full path of the '$(SPHINXBUILD)' executable. Alternatively you can add the directory with the executable to your PATH. If you don\'t have Sphinx installed, grab it from http://sphinx-doc.org/)
endif

.PHONY: help clean html analysis-html tutorials rst publish

html: Makefile rst
	@$(SPHINXBUILD) -M html "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

# Bake the MSM analysis tutorials: lifts the default exclude on
# tutorials/analysis/*.md so wget downloads run and cells execute. Expect
# several hours and ~5 GB of disk.
analysis-html: Makefile rst
	BAKE_ANALYSIS=1 $(SPHINXBUILD) -M html "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

clean:
	rm -rf $(REFERENCEDIR)/htmd.* $(REFERENCEDIR)/modules.rst
	rm -rf $(TUTORIALSDIR)/*.ipynb $(TUTORIALSDIR)/system-prep $(TUTORIALSDIR)/adaptive $(TUTORIALSDIR)/analysis
	rm -rf $(TUTORIALSDIR)/img $(TUTORIALSDIR)/**/img
	rm -rf $(HOWTODIR)/*.ipynb $(HOWTODIR)/img $(HOWTODIR)/wrapping_img
	rm -f $(SOURCEDIR)/version $(SOURCEDIR)/llms-full.txt
	rm -rf $(BUILDDIR)

publish:
	ssh software.acellera.com 'rm -rf /var/www/software.acellera.com/source/htmd/$(HTMDBRANCH)/*'
	scp -r $(SOURCEDIR)/* www.acellera.com:/var/www/software.acellera.com/source/htmd/$(HTMDBRANCH)/

# Copy .ipynb tutorials into their Diátaxis destinations. myst-nb consumes
# the .ipynb files directly (see source_suffix in conf.py); cached outputs
# render as-is because nb_execution_mode = "off".
tutorials:
	@mkdir -p $(TUTORIALSDIR)/adaptive $(TUTORIALSDIR)/analysis $(HOWTODIR)
	@for nb in $(TUT_ADAPTIVE); do cp ../tutorials/$$nb $(TUTORIALSDIR)/adaptive/$$nb; done
	@for nb in $(TUT_ANALYSIS); do cp ../tutorials/$$nb $(TUTORIALSDIR)/analysis/$$nb; done
	@for nb in $(HT_NOTEBOOKS); do cp ../tutorials/$$nb $(HOWTODIR)/$$nb; done
	@if [ -d ../tutorials/img ]; then \
		cp -R ../tutorials/img $(TUTORIALSDIR)/adaptive/; \
		cp -R ../tutorials/img $(TUTORIALSDIR)/analysis/; \
		cp -R ../tutorials/img $(HOWTODIR)/; \
	fi
	@if [ -d ../tutorials/wrapping_img ]; then cp -R ../tutorials/wrapping_img $(HOWTODIR)/; fi

rst: tutorials
	@mkdir -p $(REFERENCEDIR)
	rm -f $(REFERENCEDIR)/htmd.* $(REFERENCEDIR)/modules.rst
	$(SPHINXAPIDOC) --no-toc --force --separate -o $(REFERENCEDIR) ../htmd/ >/dev/null
	echo "$(HTMDTAG)" > $(SOURCEDIR)/version
	@echo
	@echo "Build finished. The RST files are on $(REFERENCEDIR)/"
