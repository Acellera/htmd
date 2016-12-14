#!/bin/bash
DIR=$(dirname $(readlink -f "$0"))
conda install $(cat $DIR/DEPENDENCIES) -y

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

conda list > /tmp/list
for T in $(cat $DIR/DEPENDENCIES); do
	VER=$( egrep -e "^$T " /tmp/list |awk '{print $2}')
	echo "Package $T at version $VER"
	echo "    - $T $VER" >> $DIR/meta.yaml
done
