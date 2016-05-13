#!/bin/bash
set -ev

# TG set command above is to fail immediately, https://docs.travis-ci.com/user/customizing-the-build/#Implementing-Complex-Build-Steps

if [ "$TRAVIS_REPO_SLUG" == "multiscalelab/htmd" ]; then
    export TYPE=basic
    export ANACONDA_TOKEN=$ANACONDA_TOKEN_BASIC
elif [ "$TRAVIS_REPO_SLUG" == "Acellera/htmd" ]; then
    export TYPE=pro
    export ANACONDA_TOKEN=$ANACONDA_TOKEN_PRO
else
    echo "Unexpected TRAVIS_REPO_SLUG = $TRAVIS_REPO_SLUG"
    exit 0
fi

export CHANNEL=acellera
echo "Uploading to channel: $CHANNEL"
anaconda -t $ANACONDA_TOKEN upload  $HOME/miniconda/conda-bld/*-64/htmd-[0-9]*.tar.bz2 -u $CHANNEL

if [ "$CROSS_COMPILE" == "1" ]; then
    conda convert -f -p win-64 $HOME/miniconda/conda-bld/linux-64/htmd-[0-9]*.tar.bz2
    anaconda -t $ANACONDA_TOKEN upload win-64/htmd-[0-9]*.tar.bz2 -u $CHANNEL
fi



