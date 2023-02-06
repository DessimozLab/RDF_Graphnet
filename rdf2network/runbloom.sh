#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=foldseek
#SBATCH --gres-flags=enforce-binding
#SBATCH --mem 50G
#SBATCH -t 8:00:00
#SBATCH --cpus-per-task=20


module load  gcc/10.4.0 cuda/11.6.2  cudnn/8.4.0.27-11.6 miniconda3/4.10.3
source /work/FAC/FBM/DBC/cdessim2/default/dmoi/condaenvs/etc/profile.d/conda.sh
#conda init bash
conda activate ML2

python bloomSTRING.py