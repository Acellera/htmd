#!/bin/bash
mkdir -p nlopt/build
mkdir -p nlopt/lib
cd nlopt/build
tar -zxf ../nlopt-2.4.2.tar.gz
cd nlopt*
echo "CXX=$CC"
echo "CC=$CC"
echo "FC=$FC"
AR=ar
LD=ld
RANLIB=ranlib
NM=nm
HOST=""
if [ "$CROSS_COMPILE" == "1" ]; then
	HOST="--host=mingw32"
  AR=x86_64-w64-mingw32-ar
  LD=x86_64-w64-mingw32-ld
  RANLIB=x86_64-w64-mingw32-ranlib
  NM=x86_64-w64-mingw32-nm
fi
echo "HOST=$HOST"

rm -rf test
mkdir test
echo "all:" > test/Makefile.in 
echo "" >> test/Makefile.in 
echo "install:" > test/Makefile.in 
./configure --prefix="$PWD/../.." $HOST CXX=$CC CC=$CC FC=$FC CFLAGS="-g" AR=$AR LD=$LD RANLIB=$RANLIB NM=$NM LDFLAGS="" || cat config.log

make && make install

