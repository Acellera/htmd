#!/bin/bash
#
#SBATCH --job-name=prerun_04960
#SBATCH --partition=normalGPU
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4000
#SBATCH --priority=None
#SBATCH --workdir=/shared/joao/prerun
#SBATCH --output=slurm.%N.%j.out
#SBATCH --error=slurm.%N.%j.err
#SBATCH --export=ACEMD_HOME,HTMD_LICENSE_FILE

trap "touch /shared/joao/prerun/htmd.queues.done" EXIT SIGTERM

cd /shared/joao/prerun
/shared/joao/prerun/run.sh