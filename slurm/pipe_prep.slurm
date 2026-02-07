#!/bin/bash

#SBATCH --job-name=prep_1YVL_WT
#SBATCH --account=lp_lacilab
#SBATCH --partition=batch
#SBATCH --clusters=wice
#SBATCH --nodes=1 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=02:00:00
#SBATCH --output=logs/mdrun_1YVL_WT_prep_%j.out
#SBATCH --error=logs/mdrun_1YVL_WT_prep_%j.err

export PATH=/data/leuven/359/vsc35916/miniconda3/envs/pipemd_test2/bin/:$PATH

cd /staging/leuven/stg_00187/users/Irati/thesis_main/
python -m scripts.main --params params/params_1YVL_WT.yml --prep