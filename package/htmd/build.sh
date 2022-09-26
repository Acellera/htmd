#!/usr/bin/env bash
# (c) 2015-2018 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#

echo "building"

export PATH=$PATH:/usr/bin/:/bin/
printenv

find htmd -type d -name __pycache__ -exec rm -rf {} \; -print || true

STARTDIR="$PWD"

echo "Installing into $PREFIX"

DIR="$SP_DIR"

if [ "$DIR" != "" ]; then
    mkdir -p "$DIR"
fi
if [ -e "$DIR" ]; then
    pwd
    ls
    cp -r htmd  $DIR/
    rm -rf $(find "$DIR" -name .git -type d)
    rm -rf $DIR/htmd/data
    rm -rf $DIR/htmd/.idea
    rm -rf $DIR/htmd/.ipynb_checkpoints
    rm -rf $DIR/htmd/Untitled*
    echo "def version():" > $DIR/htmd/version.py
    echo "    return \"$PKG_VERSION\"" >> $DIR/htmd/version.py
else
    echo "Error: PREFIX not defined"
    exit 1
fi

cd "$DIR/../../"

chmod -R a+rX "$PREFIX"
