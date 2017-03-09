#!/bin/bash
# (c) 2015-2017 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#

export ANACONDA_TOKEN=$ANACONDA_TOKEN_BASIC
export CHANNEL=acellera

for PACKAGE_NAME in htmd; do
    echo "Uploading to channel: $CHANNEL : PACKAGE $PACKAGE_NAME"

    if [ "$OSNAME" == "Windows" ]; then
        conda convert -f -p win-64 $HOME/miniconda/conda-bld/linux-64/$PACKAGE_NAME-[0-9]*.tar.bz2
        anaconda -t $ANACONDA_TOKEN upload win-64/$PACKAGE_NAME-[0-9]*.tar.bz2 -u $CHANNEL
    elif [ "$OSNAME" == "Linux" ]; then
        anaconda -t $ANACONDA_TOKEN upload  $HOME/miniconda/conda-bld/linux-64/$PACKAGE_NAME-[0-9]*.tar.bz2 -u $CHANNEL
    elif [ "$OSNAME" == "Darwin" ]; then
        anaconda -t $ANACONDA_TOKEN upload  $HOME/miniconda/conda-bld/osx-64/$PACKAGE_NAME-[0-9]*.tar.bz2 -u $CHANNEL
    fi
done


if [ "$MAKE_NOARCH" == "1" ]; then
    # Loop in case we need to add more noarch packages
    for PACKAGE_NAME in htmd-data htmd-deps; do
        echo "Uploading to channel: $CHANNEL : PACKAGE $PACKAGE_NAME"
        anaconda -t $ANACONDA_TOKEN upload  $HOME/miniconda/conda-bld/*/$PACKAGE_NAME-[0-9]*.tar.bz2 -u $CHANNEL
    done
fi


