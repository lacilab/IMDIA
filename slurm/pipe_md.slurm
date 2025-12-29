#!/bin/bash

#SBATCH --job-name=1BF5_WT_rerun
#SBATCH --account=lp_lacilab
#SBATCH --cluster=wice
#SBATCH --partition=gpu_h100
#SBATCH --nodes=1
#SBATCH --gpus-per-node=2
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=16
#SBATCH --mem=12G
#SBATCH --time=00:05:00
#SBATCH --output=logs/mdrun_1BF5_WT_rerun_%j.out
#SBATCH --error=logs/mdrun_1BF5_WT_rerun_%j.err


module purge
module load cluster/wice/gpu_h100
module load GROMACS/2024.4-foss-2023b-CUDA-12.4.0
module load PyYAML

cd /lustre1/scratch/359/vsc35916/thesis_main/
python -m scripts.main --md
