#!/bin/bash

#SBATCH --job-name=post_1BF5_E411A_rep1
#SBATCH --account=lp_lacilab
#SBATCH --cluster=wice
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=00:30:00
#SBATCH --output=./logs/mdrun_1BF5_E411A_rep1_post_%j.out
#SBATCH --error=./logs/mdrun_1BF5_E411A_rep1_post_%j.err
module purge
module load cluster/wice/batch
module load GROMACS/2024.4-foss-2023b-CUDA-12.4.0
module load PyYAML/6.0.1-GCCcore-13.2.0

cd /staging/leuven/stg_00187/users/Irati/thesis_main
python -m scripts.main --params params/params_1BF5_E411A.yml --post
