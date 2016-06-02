#!/bin/bash
set -ev

# TG set command above is to fail immediately, https://docs.travis-ci.com/user/customizing-the-build/#Implementing-Complex-Build-Steps

export ANACONDA_TOKEN=$ANACONDA_TOKEN_BASIC


export CHANNEL=acellera
echo "Uploading to channel: $CHANNEL : PACKAGE $PACKAGE_NAME"

if [ "$CROSS_COMPILE" == "1" ]; then
    conda convert -f -p win-64 $HOME/miniconda/conda-bld/linux-64/$PACKAGE_NAME-[0-9]*.tar.bz2
    anaconda -t $ANACONDA_TOKEN upload win-64/$PACKAGE_NAME-[0-9]*.tar.bz2 -u $CHANNEL
else
		anaconda -t $ANACONDA_TOKEN upload  $HOME/miniconda/conda-bld/*-64/$PACKAGE_NAME-[0-9]*.tar.bz2 -u $CHANNEL
fi



