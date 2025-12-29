#!/bin/bash

#SBATCH --job-name=prep_1BF5_E411A
#SBATCH --account=lp_lacilab
#SBATCH --partition=batch
#SBATCH --clusters=wice
#SBATCH --nodes=1 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=18
#SBATCH --time=02:00:00
#SBATCH --output=logs/mdrun_1BF5_E411A_prep_%j.out
#SBATCH --error=logs/mdrun_1BF5_E411A_prep_%j.err

export PATH=/lustre1/scratch/359/vsc35916/conda_envs/pipemd/bin/:$PATH

cd /data/leuven/359/vsc35916/thesis_main
python -m scripts.main --params scripts/params_1BF5.yml --prep
