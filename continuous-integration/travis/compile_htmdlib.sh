#!/usr/bin/env bash
# (c) 2015-2017 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#

# Attempt to clone (needs these environment variables to be set in Travis)
git clone https://$GITHUB_HTMDLIB_USERNAME:$GITHUB_HTMDLIB_PASSWORD@github.com/Acellera/htmdlib --depth 1

# If clone is successful
if [ "$?" == "0" ]; then
    # Remove compiled DSOs present in the repo
    rm htmd/lib/*/*.so
    # Compile the DSOs
    htmdlib/C/build.sh $PWD/htmd/lib/$OSNAME/
fi