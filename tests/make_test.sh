#!/bin/sh
name="$1"
let idx=0

if [ "$name" == "" ]; then
    echo "Syntax: $0 test-name"
    exit 1
fi

mkdir -p tests

fn=$(printf "%05d" $idx)
    while [ -e tests/$fn-* ]; do
        let idx=idx+1
        fn=$(printf "%05d" $idx)
    done

name=$(echo "$name" |  sed 's/[ \/]/_/g')
fn="$fn-$name"

mkdir  -p tests/$fn
mkdir  -p tests/$fn/reference

echo "#!/usr/bin/env python
import os
from htmd import *
from htmd.molecule import *

DATADIR=os.environ['DATADIR']

" >> tests/$fn/test.py

chmod +x tests/$fn/test.py

echo "$fn"

git add tests/$fn/test.py
git add tests/$fn/reference


        
