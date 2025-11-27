#!/bin/bash


trap "touch /home/sdoerr/Work/htmd/tutorials/adaptivemd/input/e3s1_e2s1p0f6/jobqueues.done" EXIT SIGTERM

export CUDA_VISIBLE_DEVICES=0

cd /home/sdoerr/Work/htmd/tutorials/adaptivemd/input/e3s1_e2s1p0f6
/home/sdoerr/Work/htmd/tutorials/adaptivemd/input/e3s1_e2s1p0f6/run.sh
mv *.xtc /home/sdoerr/Work/htmd/tutorials/adaptivemd/data/e3s1_e2s1p0f6