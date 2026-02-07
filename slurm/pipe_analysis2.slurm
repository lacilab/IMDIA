#!/bin/bash

#SBATCH --job-name=E411A_rep1_qc
#SBATCH --account=lp_lacilab
#SBATCH --cluster=wice
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=14
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --output=logs/mdrun_E411A_rep1_qc_%j.out
#SBATCH --error=logs/mdrun_E411A_rep1_qc_%j.err

module purge
module load cluster/wice/batch
module load GROMACS/2024.4-foss-2023b-CUDA-12.4.0
export PATH=/data/leuven/359/vsc35916/miniconda3/envs/pipemd_test2/bin/:$PATH

cd /staging/leuven/stg_00187/users/Irati/thesis_main/
python -m scripts.main --params params/params_1BF5_E411A.yml --analysis --rules qc
