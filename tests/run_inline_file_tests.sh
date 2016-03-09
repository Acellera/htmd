#!/bin/sh
RET=0
TF=$(mktemp)
for T in $(find htmd -name "*.py" -type f); do 
	python $T  > $T 2>&1
	if [ "$?" != "0" ]; then
		echo "$T   Failed"
		cat $T
		RET=1
#	else
#		echo "$T   Passed"
	fi
done
rm $TF
if [ "$RET" == "0" ]; then
	echo "All inline file tests passwd"
fi

exit $RET
