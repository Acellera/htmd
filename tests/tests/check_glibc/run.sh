#!/bin/bash

strings ../../../htmd/lib/*/Linux/*so | grep -q GLIBC_2.14
if [ "$?" == "0" ]; then
	echo "GLIBC 2.14 use detected in Linux dsos"
	exit 1
fi
exit 0
