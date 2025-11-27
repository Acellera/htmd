#!/bin/bash


trap "touch /home/sdoerr/Work/htmd/tutorials/adaptivemd/input/e1s1_ntl9_1ns_0/jobqueues.done" EXIT SIGTERM

export CUDA_VISIBLE_DEVICES=0

cd /home/sdoerr/Work/htmd/tutorials/adaptivemd/input/e1s1_ntl9_1ns_0
/home/sdoerr/Work/htmd/tutorials/adaptivemd/input/e1s1_ntl9_1ns_0/run.sh
mv *.xtc /home/sdoerr/Work/htmd/tutorials/adaptivemd/data/e1s1_ntl9_1ns_0