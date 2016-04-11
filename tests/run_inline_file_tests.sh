#!/bin/bash
RET=0
TF=$(mktemp)
for T in $(find htmd -name "*.py" -type f); do 
	python $T 2>&1
    echo "$?"
	if [ "$?" == "1" ]; then
		echo "$T   Failed"
		cat $T
		RET=1
#	else
#		echo "$T   Passed"
	fi
done
rm $TF
if [ "$RET" == "0" ]; then
	echo "All inline file tests passed"
fi

exit $RET
