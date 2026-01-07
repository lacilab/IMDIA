#!/bin/bash

#SBATCH --job-name=prep_1BF5_N460H
#SBATCH --account=lp_lacilab
#SBATCH --partition=batch
#SBATCH --clusters=wice
#SBATCH --nodes=1 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=18
#SBATCH --time=04:00:00
#SBATCH --output=logs/mdrun_1BF5_N460H_prep_%j.out
#SBATCH --error=logs/mdrun_1BF5_N460H_prep_%j.err

export PATH=/lustre1/scratch/359/vsc35916/conda_envs/pipemd/bin/:$PATH

cd /data/leuven/359/vsc35916/thesis_main
python -m scripts.main --params params/params_1BF5_N460H.yml --prep
