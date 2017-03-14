#!/bin/bash
# (c) 2015-2017 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
#DIR=$(dirname $(readlink -f "$0"))
DIR="$PWD/package/htmd-deps/"

echo "DIR is [$DIR]"
echo "PWD is [$PWD]"

conda install --file $DIR/DEPENDENCIES -y

# Just in case this got installed (eg acecloud-client,acemd depend on it)
conda remove htmd -y

python $DIR/write_meta_yaml.py $DIR
exit $?
