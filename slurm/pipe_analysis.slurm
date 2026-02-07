#!/bin/bash

#SBATCH --job-name=an_N460H_gbsa_rep1
#SBATCH --account=lp_lacilab
#SBATCH --cluster=wice
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=36
#SBATCH --time=16:00:00
#SBATCH --output=logs/mdrun_an_N460H_gbsa_rep1_%j.out
#SBATCH --error=logs/mdrun_an_N460H_gbsa_rep1_%j.err

module purge
module load cluster/wice/batch
module load GROMACS/2024.4-foss-2023b-CUDA-12.4.0
export PATH=/data/leuven/359/vsc35916/miniconda3/envs/pipemd_test2/bin/:$PATH

cd /staging/leuven/stg_00187/users/Irati/thesis_main/
python -m scripts.main --params params/params_1BF5_N460H.yml --analysis --rules gbsa