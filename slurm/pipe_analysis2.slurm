#!/bin/bash

#SBATCH --job-name=an_T385M_wb_contacts
#SBATCH --account=lp_lacilab
#SBATCH --cluster=wice
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=3
#SBATCH --cpus-per-task=12
#SBATCH --time=24:00:00
#SBATCH --output=logs/mdrun_T385M_an_wb_contacts%j.out
#SBATCH --error=logs/mdrun_T385M_an_wb_contacts%j.err

module purge
module load cluster/wice/batch
module load GROMACS/2024.4-foss-2023b-CUDA-12.4.0
export PATH=/lustre1/scratch/359/vsc35916/conda_envs/pipeanalysis/bin/:$PATH

cd /lustre1/scratch/359/vsc35916/thesis_main/
python -m scripts.main --params scripts/params_T385M.yml --analysis --rules wb contacts

