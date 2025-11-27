#!/bin/bash


trap "touch /home/sdoerr/Work/htmd/tutorials/adaptivemd/input/e3s2_e1s1p0f1/jobqueues.done" EXIT SIGTERM

export CUDA_VISIBLE_DEVICES=0

cd /home/sdoerr/Work/htmd/tutorials/adaptivemd/input/e3s2_e1s1p0f1
/home/sdoerr/Work/htmd/tutorials/adaptivemd/input/e3s2_e1s1p0f1/run.sh
mv *.xtc /home/sdoerr/Work/htmd/tutorials/adaptivemd/data/e3s2_e1s1p0f1