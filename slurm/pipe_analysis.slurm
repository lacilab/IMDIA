#!/bin/bash

#SBATCH --job-name=an_all_hbonds
#SBATCH --account=lp_lacilab
#SBATCH --cluster=wice
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=9
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00
#SBATCH --output=logs/mdrun_an_all_hbonds%j.out
#SBATCH --error=logs/mdrun_an_all_hbonds%j.err

module purge
module load cluster/wice/batch
module load GROMACS/2024.4-foss-2023b-CUDA-12.4.0
export PATH=/lustre1/scratch/359/vsc35916/conda_envs/pipeanalysis/bin/:$PATH

cd /lustre1/scratch/359/vsc35916/thesis_main/
python -m scripts.main --params scripts/params.yml --analysis --rules hbonds
