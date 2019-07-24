#!/bin/bash
# This script is meant to be called by the "after_success" step defined in
# .travis.yml. See https://docs.travis-ci.com/ for more details.

# License: 3-clause BSD

# Ok, now that we succeeded, I want this to be published somewhere not too
# public. How can I do that?

set -e


# commands taken shamelessly from https://gist.github.com/willprice/e07efd73fb7f13f917ea
setup_git() {
  git config --global user.email "travis@travis-ci.org"
  git config --global user.name "Travis CI"
}

commit_website_files() {
  # Use force to add the files and commit them.
  git add man/
  git commit --message "Travis building documentation: $TRAVIS_BUILD_NUMBER"  || echo "No changes to commit"
}

upload_files() {
  git push origin master 
}

if [ "$TRAVIS_PULL_REQUEST" = "false" ]; then
    # cp -r scripts/reports/*.html /tmp/html_outputs
    # cp -r scripts/reports/*.pdf /tmp/html_outputs
    # cd /tmp/html_outputs
    setup_git
    commit_website_files
    upload_files
fi;
