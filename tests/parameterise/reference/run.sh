#!/bin/sh
#PBS -J1-119
#PBS -lselect=1:ncpus=1:mem=2gb
#PBS -lwalltime=72:0:0


if "$PBS_O_WORKDIR" != ""; then
 cd $PBS_O_WORKDIR
fi

export PYTHONPATH="$PWD/../../.."

newparameterize="python $PWD/../../../htmd/newparameterization/cli.py"

if "$PBS_ARRAY_INDEX" != ""; then
 DD=$(head -$PBS_ARRAY_INDEX list | tail -1)
 cd $DD
 $newparameterize input.mol2 > log.txt 2>&1
else

for T in m*; do
  echo "PARAMETERIZE $T"
  cd $T
  $newparameterize input.mol2 > log.txt 2>&1
  cd -
done

fi
