#!/bin/bash
# (c) 2015-2016 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#

export ANACONDA_TOKEN=$ANACONDA_TOKEN_BASIC

for PACKAGE_NAME in htmd htmd-deps htmd-data; do
export CHANNEL=acellera
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




