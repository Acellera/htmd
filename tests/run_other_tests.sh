#!/bin/sh

R=0

for T in tests/*/run.sh; do
	S="$PWD"
	DD=$(dirname $T)
	cd $DD
	./run.sh
	RET=$?
	echo "Test: $DD : $RET"
	if [ $RET != "0" ]; then
		R=1
	fi	
done

exit $R
