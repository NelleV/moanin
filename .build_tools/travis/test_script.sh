#!/bin/bash
# This script is meant to be called by the "script" step defined in
# .travis.yml. See http://docs.travis-ci.com/ for more details.
# The behavior of the script is controlled by environment variabled defined
# in the .travis.yml in the top level folder of the project.

# License: 3-clause BSD

set -e

install_github_dependencies() {
    # install_github keeps failing because of the Github API rate limit. So
    # we're just going to do this by hand…
    pushd /tmp
    wget -P . https://github.com/NelleV/timecoursedata/archive/master.zip
    unzip /tmp/master.zip
    cd timecoursedata-master
    make install-extra
    make install
    popd
}

run_tests() {
    # This runs simple "normal" + bioconductor checks on the package
    # make check
    make test
    make doc
    make check
}

# Start by installing
install_github_dependencies
make install-extra
make install

# Now get data and run the tests and build the manuscript
run_tests
