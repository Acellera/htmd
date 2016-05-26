#!/usr/bin/env bash
if [ "$TRAVIS_OS_NAME" == "linux" ]; then
	parameterize --delete ethanol
	parameterize --input input
fi

