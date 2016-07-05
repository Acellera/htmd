#!/bin/sh

for T in m*; do
  cd $T
  newparameterize input.mol2
  cd -
done
