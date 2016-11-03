#!/bin/sh

PACKAGES="
ambermini
bhmm
funcsigs
mdtraj
msmtools
openbabel
pint
progress_reporter
pyemma
thermotools
"

DIR=$(mktemp -d)
cd $DIR
echo "Using TMPDIR $DIR"

mkdir win-64 win-32 osx-32 osx-64 linux-32 linux-64 noarch

for T in $PACKAGES; do
	echo "Downloading [$T]"
	anaconda download omnia/$T &
done	
wait

for F in noarch/* */*py35*; do
	echo "Uploading [$F]"
<<<<<<< HEAD
	anaconda upload $F -u acellera  &
=======
	anaconda upload $F -u acellera 
>>>>>>> d702fe0d1faf3f424d82174ee1284c8ad2e42a58
done
wait

cd /
rm -rf $DIR

exit 0 

