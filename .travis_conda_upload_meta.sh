#!/usr/bin/env bash

export ANACONDA_TOKEN=$ANACONDA_TOKEN_BASIC

export MAJOR=$(echo $TRAVIS_BRANCH | sed 's/[\.-]/ /g'  | awk '{print $1}')
export MINOR=$(echo $TRAVIS_BRANCH | sed 's/[\.-]/ /g'  | awk '{print $2}')
export BUGFIX=$(echo $TRAVIS_BRANCH | sed 's/[\.-]/ /g'  | awk '{print $3}')

echo "$MINOR" | egrep -qe "[02468]$"
if [ "$?" == "0" ]; then
 export META_PACKAGE_NAME=htmd-stable
else
 export META_PACKAGE_NAME=htmd-devel
fi

export CHANNEL=acellera
echo "Uploading to channel: $CHANNEL ; META_PACKAGE: $META_PACKAGE_NAME (based on $PACKAGE_NAME version $MAJOR.$MINOR.$BUGFIX)"

if [ "$CROSS_COMPILE" == "1" ]; then
    conda convert -f -p win-64 $HOME/miniconda/conda-bld/linux-64/$PACKAGE_NAME-[0-9]*.tar.bz2
    export pkgtar=$(ls win-64/$META_PACKAGE_NAME-[0-9]*.tar.bz2)
    cp -rp $pkgtar ${pkgtar/$PACKAGE_NAME/$META_PACKAGE_NAME}
    anaconda -t $ANACONDA_TOKEN upload ${pkgtar/$PACKAGE_NAME/$META_PACKAGE_NAME} -u $CHANNEL
else
    export pkgtar=$(ls $HOME/miniconda/conda-bld/*-64/$META_PACKAGE_NAME-[0-9]*.tar.bz2)
    cp -rp $pkgtar ${pkgtar/$PACKAGE_NAME/$META_PACKAGE_NAME}
	anaconda -t $ANACONDA_TOKEN upload ${pkgtar/$PACKAGE_NAME/$META_PACKAGE_NAME} -u $CHANNEL
fi