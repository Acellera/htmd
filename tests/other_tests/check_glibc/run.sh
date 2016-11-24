#!/bin/bash



strings ../../../htmd/lib/*/Linux/*so ../../../htmd/lib/Linux/*so | egrep -q -e  'GLIBC_2.1[45]'
if [ "$?" == "0" ]; then
	echo "GLIBC 2.14 use detected in Linux dsos"
	exit 1
fi
exit 0
