#!/bin/bash
# (c) 2015-2017 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
if [ "$OSNAME" == "Linux" ]; then
	DIR=$(dirname $(readlink -f $(which python)))/../lib/python3.5/site-packages/htmd/lib
  for T in $DIR/Linux/*so; do
		echo $T
	  strings $T | egrep -e  'GLIBC_2.1[45]'
	done
	strings $DIR/Linux/*so | egrep -q -e  'GLIBC_2.1[45]'
	if [ "$?" == "0" ]; then
		echo "GLIBC 2.14 use detected in Linux dsos"
		exit 1
	fi
	exit 0
fi

