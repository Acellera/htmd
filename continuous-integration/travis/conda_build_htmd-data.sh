#!/usr/bin/env bash
# (c) 2015-2018 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#

export MAJOR=0
export MINOR=1
export SHORTHASH=$(git --git-dir=./htmd/.git --work-tree=./htmd ls-tree HEAD htmd/data | awk '{print substr($3,0,7)}')

export HTMD_DATA_VERSION="${MAJOR}.${MINOR}.${SHORTHASH}"

conda build package/htmd-data --no-include-recipe
