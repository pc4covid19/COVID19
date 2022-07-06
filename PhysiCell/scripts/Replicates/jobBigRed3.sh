#!/bin/bash

#SBATCH --mail-user=getzm@iu.edu
#SBATCH --job-name=COVID19
#SBATCH -p general
#SBATCH -o COVID19_%j.txt
#SBATCH -e COVID19_%j.err
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=3
#SBATCH --cpus-per-task=1
#SBATCH --time=1-0:00:00
#SBATCH --mail-type=FAIL,BEGIN,END

module load python/3.6.11
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

srun --cpu-bind=sockets python RunModelReplicate.py
