#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=bloomfilter
#SBATCH --gres-flags=enforce-binding
#SBATCH --mem 50G
#SBATCH -t 8:00:00
#SBATCH --cpus-per-task=20


source /work/FAC/FBM/DBC/cdessim2/default/dmoi/miniconda3/etc/profile.d/conda.sh
#conda init bash
conda activate ML2

python bloomSTRING.py