RSCRIPT ?= Rscript


all: doc-doxygen check build install test

build:
	$(RSCRIPT) .scripts/build_package.R

install: doc-doxygen
	$(RSCRIPT) .scripts/install_package.R


check:
	$(RSCRIPT) .scripts/check_package.R

doc-doxygen:
	$(RSCRIPT) .scripts/build_documentation.R

doc-html: install
	$(RSCRIPT) .scripts/build_documentation_html.R

doc: doc-html

test:
	$(RSCRIPT) .scripts/test_package.R
