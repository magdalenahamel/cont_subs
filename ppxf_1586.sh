#!/bin/bash
#
#SBATCH --job-name=ppxf_single_spec
#SBATCH --output=ppxf_fitting_results/output_sigle_spec.txt
#
#SBATCH --ntasks=1
#SBATCH --time=5:00
#SBATCH --mem=40G

#SBATCH --mail-type=ALL
#SBATCH --mail-user=magdalena.hamel@gmail.com

#SBATCH --account=oz088

module load anaconda3/5.0.1

srun python python_scripts/apply_ppxf.py