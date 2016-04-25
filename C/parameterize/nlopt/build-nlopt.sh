#!/bin/sh
mkdir -p nlopt/build
mkdir -p nlopt/lib
cd nlopt/build
tar -zxf ../nlopt-2.4.2.tar.gz
cd nlopt*
echo "CXX=$CC"
echo "CC=$CC"
echo "FC=$FC"
HOST=""
if [ "$CROSS_COMPILE" == "1" ]; then
	HOST="--host=mingw32"
fi
echo "HOST=$HOST"

./configure --prefix="$PWD/../.." $HOST CXX=$CC CC=$CC FC=$FC CFLAGS="-g" LDFLAGS="" || cat config.log

make && make install

