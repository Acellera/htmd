#!/usr/bin/env bash

export MAJOR=$(git describe --tags | sed 's/[\.-]/ /g'  | awk '{print $1}')
export MINOR=$(git describe --tags | sed 's/[\.-]/ /g'  | awk '{print $2}')
export BUGFIX=$(git describe --tags | sed 's/[\.-]/ /g'  | awk '{print $3}')

echo "$MINOR" | egrep -qe "[02468]$"
if [ "$?" == "0" ]; then
 export PACKAGE_NAME=htmd-stable
else
 export PACKAGE_NAME=htmd-devel
fi

export BUGFIX
export MINOR_VERSION="${MAJOR}.${MINOR}"
export BUGFIX_VERSION="${MAJOR}.${MINOR}.${BUGFIX}"

echo "Building package $PACKAGE_NAME version $MINOR_VERSION bugfix $BUGFIX"

conda build . --no-include-recipe

if [ "$ANACONDA_TOKEN" != "" ]; then
	anaconda -t $ANACONDA_TOKEN upload $HOME/miniconda/conda-bld/*-64/$PACKAGE_NAME-*.tar.bz2 -u acellera
fi
