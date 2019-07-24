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
  git add . *.html
  git commit --message "Travis build: $TRAVIS_BUILD_NUMBER"
}

upload_files() {
  git remote add origin-pages https://${GH_TOKEN}@github.com/NelleV/2019timecourse-rnaseq-pipeline_html_outputs.git > /dev/null 2>&1
  git push --quiet --set-upstream origin-pages master 
}

if [ "$TRAVIS_PULL_REQUEST" = "false" ]; then
    # git clone https://${GH_TOKEN}@github.com/NelleV/2019timecourse-rnaseq-pipeline_html_outputs.git /tmp/html_outputs
    # cp -r scripts/reports/*.html /tmp/html_outputs
    # cp -r scripts/reports/*.pdf /tmp/html_outputs
    cd /tmp/html_outputs
    # setup_git
    # commit_website_files
    # upload_files
fi;
