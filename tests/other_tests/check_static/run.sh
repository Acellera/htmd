#!/bin/bash

which ldd
if [ "$?" == "1" ]; then
	# No ldd -- infer that this is Travis OS X
	exit 0
fi

RET=0
for T in  $(dirname $(which parameterize))/../lib/*/*/htmd/parameterization/share/bin/* ; do
	ldd $T | grep -q "not a dynamic executable"
	if [ "$?" == "1" ]; then
		RET=1
	fi
done

exit $RET
