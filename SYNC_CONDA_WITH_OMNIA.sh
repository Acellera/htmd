#!/bin/sh

PACKAGES="
mdtraj
openbabel
progress_reporter
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
	anaconda upload $F -u acellera 
done

cd /
rm -rf $DIR

exit 0 

