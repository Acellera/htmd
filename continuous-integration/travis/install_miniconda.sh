#!/bin/bash
# taken from https://raw.githubusercontent.com/broadinstitute/viral-ngs/master/travis/install-conda.sh
set -e

# the miniconda directory may exist if it has been restored from cache
if [ -d "$MINICONDA_DIR" ] && [ -e "$MINICONDA_DIR/bin/conda" ]; then
    echo "Miniconda install already present from cache: $MINICONDA_DIR"
    export PATH="$MINICONDA_DIR/bin:$PATH"
    hash -r
else # if it does not exist, we need to install miniconda
    rm -rf "$MINICONDA_DIR" # remove the directory in case we have an empty cached directory

    if [[ "$TRAVIS_OS_NAME" == "osx" ]]; then
        wget --quiet https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh -O miniconda.sh;
    else
        wget --quiet https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh;
    fi

    bash miniconda.sh -b -p "$MINICONDA_DIR"
    source "$MINICONDA_DIR/etc/profile.d/conda.sh"
    hash -r

    conda config --set always_yes yes \
        --set changeps1 no \
        --set quiet yes \
        --set auto_update_conda false \
        --system # important to write to system cfg, otherwise we loose the changes upon cache reloading.
fi

# we want to have an up to date conda-build.
conda install -q conda-build=3
conda info -a # for debugging
conda create -q -n test python=$TRAVIS_PYTHON_VERSION
echo "python $TRAVIS_PYTHON_VERSION*" > $MINICONDA_DIR/envs/test/conda-meta/pinned