#!/usr/bin/env bash
# Run this
make clean

for TYPE in pro basic; do
#for TRAVIS_OS_NAME in "linux" "windows"; do #  "windows" "osx"; do
	if [ "$1" == "linux"   ]; then
		PLATFORM=Linux; 
		CC=gcc
		CXX=g++	
	elif [ "$1" == "windows" ]; then 
		PLATFORM=Windows;
		CC=x86_64-w64-mingw32-gcc
		CXX=x86_64-w64-mingw32-g++
	elif [ "$1" == "osx"     ]; then 
		PLATFORM=Darwin; 
		CC=gcc
		CXX=g++
	else
		echo "Unknown OS [$1]"
		exit 1
	fi
	make all CC=$CC CXX=$CXX PLATFORM=$PLATFORM  TYPE=$TYPE
#done
done
