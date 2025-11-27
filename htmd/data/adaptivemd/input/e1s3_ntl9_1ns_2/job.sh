#!/bin/bash


trap "touch /home/sdoerr/Work/htmd/tutorials/adaptivemd/input/e1s3_ntl9_1ns_2/jobqueues.done" EXIT SIGTERM

export CUDA_VISIBLE_DEVICES=0

cd /home/sdoerr/Work/htmd/tutorials/adaptivemd/input/e1s3_ntl9_1ns_2
/home/sdoerr/Work/htmd/tutorials/adaptivemd/input/e1s3_ntl9_1ns_2/run.sh
mv *.xtc /home/sdoerr/Work/htmd/tutorials/adaptivemd/data/e1s3_ntl9_1ns_2