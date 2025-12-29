#!/bin/bash

#SBATCH --job-name=post_WT_rep2
#SBATCH --account=lp_lacilab
#SBATCH --cluster=wice
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=1:00:00
#SBATCH --output=logs/mdrun_WT_rep2_post_%j.out
#SBATCH --error=logs/mdrun_WT_rep2_post_%j.err


module purge
module load cluster/wice/gpu_h100
module load GROMACS/2024.4-foss-2023b-CUDA-12.4.0
module load PyYAML

cd /lustre1/scratch/359/vsc35916/thesis_main/
python -m scripts.main --post
