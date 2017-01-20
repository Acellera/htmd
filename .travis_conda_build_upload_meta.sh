#!/usr/bin/env bash
# (c) 2015-2016 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#

export ANACONDA_TOKEN=$ANACONDA_TOKEN_BASIC

if [ "$ANACONDA_TOKEN" != "" ]; then
	ANACONDA_TOKEN="-t $ANACONDA_TOKEN"
fi

export MAJOR=$(echo $TRAVIS_BRANCH | sed 's/[\.-]/ /g'  | awk '{print $1}')
export MINOR=$(echo $TRAVIS_BRANCH | sed 's/[\.-]/ /g'  | awk '{print $2}')
export BUGFIX=$(echo $TRAVIS_BRANCH | sed 's/[\.-]/ /g'  | awk '{print $3}')

echo "$MINOR" | egrep -qe "[02468]$"
if [ "$?" == "0" ]; then
 export PACKAGE_NAME=htmd-stable
else
 export PACKAGE_NAME=htmd-latest
fi

export BUGFIX
export MINOR_VERSION="${MAJOR}.${MINOR}"
export BUGFIX_VERSION="${MAJOR}.${MINOR}.${BUGFIX}"

conda build package/htmd-meta --no-include-recipe

export CHANNEL=acellera
echo "Uploading to channel: $CHANNEL ; PACKAGE: $PACKAGE_NAME (based on $TRAVIS_BRANCH version $MAJOR.$MINOR.$BUGFIX)"

if [ "$CROSS_COMPILE" == "1" ]; then
    conda convert -f -p win-64 $HOME/miniconda/conda-bld/linux-64/$PACKAGE_NAME-[0-9]*.tar.bz2
    anaconda $ANACONDA_TOKEN upload win-64/$PACKAGE_NAME-[0-9]*.tar.bz2 -u $CHANNEL -p $PACKAGE_NAME
else
	anaconda $ANACONDA_TOKEN upload $HOME/miniconda/conda-bld/*-64/$PACKAGE_NAME-*.tar.bz2 -u $CHANNEL -p $PACKAGE_NAME
fi
