#!/bin/bash


trap "touch /home/sdoerr/Work/htmd/tutorials/adaptivemd/input/e2s2_e1s2p0f3/jobqueues.done" EXIT SIGTERM

export CUDA_VISIBLE_DEVICES=0

cd /home/sdoerr/Work/htmd/tutorials/adaptivemd/input/e2s2_e1s2p0f3
/home/sdoerr/Work/htmd/tutorials/adaptivemd/input/e2s2_e1s2p0f3/run.sh
mv *.xtc /home/sdoerr/Work/htmd/tutorials/adaptivemd/data/e2s2_e1s2p0f3