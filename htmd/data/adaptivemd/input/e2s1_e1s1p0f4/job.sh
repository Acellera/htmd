#!/bin/bash


trap "touch /home/sdoerr/Work/htmd/tutorials/adaptivemd/input/e2s1_e1s1p0f4/jobqueues.done" EXIT SIGTERM

export CUDA_VISIBLE_DEVICES=0

cd /home/sdoerr/Work/htmd/tutorials/adaptivemd/input/e2s1_e1s1p0f4
/home/sdoerr/Work/htmd/tutorials/adaptivemd/input/e2s1_e1s1p0f4/run.sh
mv *.xtc /home/sdoerr/Work/htmd/tutorials/adaptivemd/data/e2s1_e1s1p0f4