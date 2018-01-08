#!/usr/bin/env bash
# (c) 2015-2017 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#

# Attempt to clone (needs these environment variables to be set in Travis; they are not available on PRs due to security reasons)
git clone https://$GITHUB_HTMDLIB_USERNAME:$GITHUB_HTMDLIB_PASSWORD@github.com/Acellera/htmdlib --depth 1

# If clone is successful
if [ "$?" == "0" ]; then
    # Check if we're on stable or latest, and choose the appropriate htmdlib branch
    export TAG_DESCRIBE=$(git --git-dir=htmd/.git describe)
    export MINOR=$(echo $TAG_DESCRIBE | sed 's/[\.-]/ /g'  | awk '{print $2}')
    echo "$MINOR" | egrep -qe "[02468]$"
    if [ "$?" == "0" ]; then
        echo "Checking out stable version of HTMDLIB"
        git --git-dir=htmdlib/.git checkout stable
    else
        echo "Checking out latest version of HTMDLIB"
        git --git-dir=htmdlib/.git checkout master
    fi
    # Remove compiled DSOs present in the repo
    rm htmd/lib/*/*.so
    # Compile the DSOs
    htmdlib/C/build.sh $PWD/htmd/lib/$OSNAME/
fi