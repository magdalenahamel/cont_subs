#!/bin/sh
#SBATCH -N 1      # nodes requested
#SBATCH -n 1      # tasks requested
#SBATCH -c 4      # cores requested
#SBATCH --mem=10  # memory in Mb
#SBATCH -o outfile  # send stdout to outfile
#SBATCH -e errfile  # send stderr to errfile
#SBATCH -t 0:01:00  # time requested in hour:minute:second

module load python
python hello.py