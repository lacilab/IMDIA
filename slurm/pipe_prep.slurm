#!/bin/bash

#SBATCH --job-name=prep_5FKW
#SBATCH --account=lp_lacilab
#SBATCH --partition=batch
#SBATCH --clusters=wice
#SBATCH --nodes=1 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=18
#SBATCH --time=01:00:00
#SBATCH --output=logs/mdrun_5FKW_prep_%j.out
#SBATCH --error=logs/mdrun_5FKW_prep_%j.err

export PATH=/data/leuven/359/vsc35916/miniconda3/envs/pipemd_test2/bin/:$PATH

cd /data/leuven/359/vsc35916/thesis_main
python -m scripts.main --params params/params_5FKW.yml --prep