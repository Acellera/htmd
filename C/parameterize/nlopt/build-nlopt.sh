#!/bin/sh
mkdir -p nlopt/build
cd nlopt/build
tar -zxvf ../nlopt-2.4.2.tar.gz
cd nlopt*
./configure --prefix="$PWD/../.."
make && make install

