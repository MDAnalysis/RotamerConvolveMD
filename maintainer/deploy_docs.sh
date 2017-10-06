#!/bin/bash
# Deploying docs from travis-ci.
# See https://github.com/MDAnalysis/mdanalysis/issues/386
# Script based on https://github.com/steveklabnik/automatically_update_github_pages_with_travis_example

# Run this script from the top-level of the checked out git
# repository. A github OAuth token must be available in the evironment
# variable GH_TOKEN and is set up through the .travis.yml
# env:global:secure parameter (encrypted with travis-ci's public key)/
#
# Additional environment variables set in .travis.yml
#  GH_REPOSITORY     repo to full from and push to
#  GH_DOC_BRANCH     branch from which the docs are built
#  GIT_CI_USER       name of the user to push docs as
#  GIT_CI_EMAIL      email of the user to push docs as
#  MDA_DOCDIR        path to the docdir from top of repo
#
# NOTE: If any of these environment variables are not set or 
#       empty then the script will exit with an error (-o nounset).

set -o errexit -o nounset

function die () {
    local msg="$1" err=${2:-1}
    echo "ERROR: $msg [$err]"
    exit $err
}

GIT_CI_USER=TravisCI
GIT_CI_EMAIL=travisCI@mdanalysis.org
MDA_DOCDIR=doc/sphinx/build/html
MDA_BUILDDIR=doc/sphinx

GH_DOC_BRANCH=master
GH_REPOSITORY=MDAnalysis/RotamerConvolveMD.git

rev=$(git rev-parse --short HEAD)

# the following tests should be superfluous because of -o nounset
#test -n "${GH_TOKEN}" || die "GH_TOKEN is empty: need OAuth GitHub token to continue" 100
test -n "${GH_REPOSITORY}" || die "GH_REPOSITORY must be set in .travis.yml" 100
test -n "${MDA_DOCDIR}" || die "MDA_DOCDIR must be set in .travis.yml" 100

(cd ${MDA_BUILDDIR} && make html) || die "Failed to build docs in ${MDA_BUILDDIR}"

cd ${MDA_DOCDIR} || die "Failed to 'cd ${MDA_DOCDIR}'. Run from the top level of the repository"

test -d .git && rm -rf .git

git init
git config user.name "${GIT_CI_USER}"
git config user.email "${GIT_CI_EMAIL}"

git remote add upstream "git@github.com:${GH_REPOSITORY}"
git fetch --depth 50 upstream ${GH_DOC_BRANCH} gh-pages
git reset upstream/gh-pages

touch .
touch .nojekyll

git add -A .
git commit -m "rebuilt html docs from branch ${GH_DOC_BRANCH} with sphinx at ${rev}"
git push -q upstream HEAD:gh-pages


