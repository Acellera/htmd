#!/bin/bash
#PBS -lselect=1:ncpous=1:mem=2gb
#PBS -lwalltime=72:0:0
#PBS -J 1-119
if [ "$PBS_O_WORKDIR" != "" ]; then
 cd "$PBS_O_WORKDIR"
fi

D=$(head -$PBS_ARRAY_INDEX list | tail -1)
cd $D

rm -rf minimize esp dihed* param*

export PYTOHNPATH="$PWD/../../../../"
python ../../../../htmd/newparameterization/cli.py input.mol2 > log.txt 2>&1

