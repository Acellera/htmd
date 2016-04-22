#!/bin/sh
mkdir -p nlopt/build
cd nlopt/build
tar -zxf ../nlopt-2.4.2.tar.gz
cd nlopt*
echo "CXX=$CC"
echo "CC=$CC"
echo "FC=$FC"
./configure --prefix="$PWD/../.." CXX=$CC CC=$CC FC=$FC --without-threads  CFLAGS="-g" || cat config.log

make && make install

