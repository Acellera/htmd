# Makefile for Sphinx documentation
#

# You can set these variables from the command line.
SPHINXOPTS    =
SPHINXBUILD   = sphinx-build
SPHINXAPIDOC  = sphinx-apidoc
IPYTHON       = jupyter
PAPER         =
BUILDDIR      = build
TUTORIALSDIR  = source/tutorials
HTMDTAG       = $(shell python gethtmdtag.py tag)
HTMDBRANCH    = $(shell python gethtmdtag.py branch)
IPYNB_FILES   = $(shell ls ../tutorials/*.ipynb)

# User-friendly check for sphinx-build
ifeq ($(shell which $(SPHINXBUILD) >/dev/null 2>&1; echo $$?), 1)
$(error The '$(SPHINXBUILD)' command was not found. Make sure you have Sphinx installed, then set the SPHINXBUILD environment variable to point to the full path of the '$(SPHINXBUILD)' executable. Alternatively you can add the directory with the executable to your PATH. If you don\'t have Sphinx installed, grab it from http://sphinx-doc.org/)
endif

# Internal variables.
PAPEROPT_a4     = -D latex_paper_size=a4
PAPEROPT_letter = -D latex_paper_size=letter
ALLSPHINXOPTS   = -d $(BUILDDIR)/doctrees $(PAPEROPT_$(PAPER)) $(SPHINXOPTS) source
# the i18n builder cannot share the environment and doctrees with the others
I18NSPHINXOPTS  = $(PAPEROPT_$(PAPER)) $(SPHINXOPTS) source

RST_FILES := $(addprefix $(TUTORIALSDIR)/,$(notdir $(IPYNB_FILES:.ipynb=.rst)))

.PHONY: help clean html dirhtml singlehtml pickle json htmlhelp qthelp devhelp epub latex latexpdf text man changes linkcheck doctest coverage gettext static

help:
	@echo "Please use \`make <target>' where <target> is one of"
	@echo "  html       to make standalone HTML files"
	@echo "  dirhtml    to make HTML files named index.html in directories"
	@echo "  singlehtml to make a single large HTML file"
	@echo "  pickle     to make pickle files"
	@echo "  json       to make JSON files"
	@echo "  htmlhelp   to make HTML files and a HTML help project"
	@echo "  qthelp     to make HTML files and a qthelp project"
	@echo "  applehelp  to make an Apple Help Book"
	@echo "  devhelp    to make HTML files and a Devhelp project"
	@echo "  epub       to make an epub"
	@echo "  latex      to make LaTeX files, you can set PAPER=a4 or PAPER=letter"
	@echo "  latexpdf   to make LaTeX files and run them through pdflatex"
	@echo "  latexpdfja to make LaTeX files and run them through platex/dvipdfmx"
	@echo "  text       to make text files"
	@echo "  man        to make manual pages"
	@echo "  texinfo    to make Texinfo files"
	@echo "  info       to make Texinfo files and run them through makeinfo"
	@echo "  gettext    to make PO message catalogs"
	@echo "  changes    to make an overview of all changed/added/deprecated items"
	@echo "  xml        to make Docutils-native XML files"
	@echo "  pseudoxml  to make pseudoxml-XML files for display purposes"
	@echo "  linkcheck  to check all external links for integrity"
	@echo "  doctest    to run all doctests embedded in the documentation (if enabled)"
	@echo "  coverage   to run coverage check of the documentation (if enabled)"

clean:
	rm -f source/htmd.*
	rm -rf $(TUTORIALSDIR)
	rm source/version

publish:
	ssh www.acellera.com 'rm -rf ~/software.acellera.com/source/htmd/$(HTMDBRANCH)/*'
	scp -r source/* www.acellera.com:~/software.acellera.com/source/htmd/$(HTMDBRANCH)/
	ssh www.acellera.com 'rm -rf ~/software.acellera.com/source/parameterize/*'
	scp -r source-parameterize/* www.acellera.com:~/software.acellera.com/source/parameterize/
	ssh www.acellera.com 'chmod -R g+rwX software.acellera.com/'

# Rule for building .rst from .ipynb
$(TUTORIALSDIR)/%.rst: ../tutorials/%.ipynb
	mkdir -p $(TUTORIALSDIR)
	$(IPYTHON) nbconvert $< --to rst --output-dir=$(TUTORIALSDIR)/

tutorials: $(RST_FILES)

rst: tutorials
	rm -f source/htmd.*
	$(SPHINXAPIDOC) --force --separate -o source ../htmd/ >/dev/null
	echo "$(HTMDTAG)" > source/version
	@echo
	@echo "Build finished. The RST files are on source/"

