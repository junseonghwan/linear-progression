#!/bin/bash -l
 
#SBATCH -A snic2018-8-281
#SBATCH -p core
#SBATCH -n 16
#SBATCH -t 24:00:00
#SBATCH -J lpm_model_selection

python run_model_selection.py $1 $2 $3
