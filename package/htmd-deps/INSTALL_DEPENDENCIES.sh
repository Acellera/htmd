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

echo "
package:
  name: htmd-deps
  version: {{ environ.get('BUILD_VERSION') }}

source:
   path: .

requirements:
  build:
    - python
    - requests


  run:
" > $DIR/meta.yaml

RET=0

conda list > /tmp/list
for T in $(for i in $(cat package/htmd-deps/DEPENDENCIES); do echo ${i%%=*}; done); do
	VER=$( egrep -e "^$T " /tmp/list |awk '{print $2}')
	echo "    - $T $VER" >> $DIR/meta.yaml
	if [ "$VER" == "" ]; then
			RET=1
			VER="[NOT FOUND]"
	fi
	echo "Package $T at version $VER"
done

exit $RET
