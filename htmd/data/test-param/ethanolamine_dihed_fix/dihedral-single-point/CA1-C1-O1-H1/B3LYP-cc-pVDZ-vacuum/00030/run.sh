#!/bin/sh

export HTMD_PSI4_WORKDIR=$(pwd)
psi4 -i psi4.in -o psi4.out &> psi4.log
