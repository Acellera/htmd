#!/usr/bin/env bash

export ANACONDA_TOKEN=$ANACONDA_TOKEN_BASIC

export MAJOR=$(echo $TRAVIS_BRANCH | sed 's/[\.-]/ /g'  | awk '{print $1}')
export MINOR=$(echo $TRAVIS_BRANCH | sed 's/[\.-]/ /g'  | awk '{print $2}')
export BUGFIX=$(echo $TRAVIS_BRANCH | sed 's/[\.-]/ /g'  | awk '{print $3}')

echo "$MINOR" | egrep -qe "[02468]$"
if [ "$?" == "0" ]; then
 export META_PACKAGE_NAME=htmd-stable
else
 export META_PACKAGE_NAME=htmd-latest
fi

export BUGFIX
export MINOR_VERSION="${MAJOR}.${MINOR}"
export BUGFIX_VERSION="${MAJOR}.${MINOR}.${BUGFIX}"

conda build --python $TRAVIS_PYTHON_VERSION package/htmd-meta

export CHANNEL=acellera
echo "Uploading to channel: $CHANNEL ; META_PACKAGE: $META_PACKAGE_NAME (based on $PACKAGE_NAME version $MAJOR.$MINOR.$BUGFIX)"

if [ "$CROSS_COMPILE" == "1" ]; then
    conda convert -f -p win-64 $HOME/miniconda/conda-bld/linux-64/$META_PACKAGE_NAME-[0-9]*.tar.bz2
    anaconda -t $ANACONDA_TOKEN upload win-64/$META_PACKAGE_NAME-[0-9]*.tar.bz2 -u $CHANNEL -p $META_PACKAGE_NAME
else
	anaconda -t $ANACONDA_TOKEN upload $HOME/miniconda/conda-bld/*-64/$META_PACKAGE_NAME-*.tar.bz2 -u $CHANNEL -p $META_PACKAGE_NAME
fi