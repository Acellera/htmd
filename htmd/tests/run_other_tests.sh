#!/bin/sh
# (c) 2015-2017 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#

R=0

for T in tests/other_tests/*/run.sh; do
	S="$PWD"
	DD=$(dirname $T)
	cd $DD
	./run.sh
	RET=$?
	echo "Test: $DD : $RET"
	if [ $RET != "0" ]; then
		R=1
	fi	
	cd "$S"
done

exit $R
