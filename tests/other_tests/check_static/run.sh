#!/bin/bash
# (c) 2015-2016 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#

which ldd
if [ "$?" == "1" ]; then
	# No ldd -- infer that this is Travis OS X
	exit 0
fi

RET=0
for T in  $(find $(dirname $(which parameterize))/../lib/*/*/htmd/parameterization/share/bin/ -type f); do
	ldd $T | grep -q "not a dynamic executable"
	if [ "$?" == "1" ]; then
		RET=1
	fi
done

exit $RET
