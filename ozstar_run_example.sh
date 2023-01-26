#!/bin/sh
#SBATCH -N 1      # nodes requested
#SBATCH -n 1      # tasks requested
#SBATCH -c 4      # cores requested
#SBATCH --mem=5G  # memory in Mb
#SBATCH -t 0:01:00  # time requested in hour:minute:second
#SBATCH --mail-type=ALL
#SBATCH --mail-user=magdalena.hamel@gmail.com

#SBATCH --account=oz088

module load anaconda3/5.0.1

python ppxf/ppxf_kinematics_example_sdss.py