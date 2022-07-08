#!/usr/bin/env bash
#SBATCH -J extractBCs
#SBATCH -p hns,dpetrov,normal,owners
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -t 0-01:00
#SBATCH --array=1-13
#SBATCH --mem-per-cpu=2G
#SBATCH --requeue
#SBATCH -o SlurmFiles/slurm-%A_%a_%x.out
#SBATCH --mail-user=grantkinsler@gmail.com
#SBATCH --mail-type=END

files=$(sed -n "$SLURM_ARRAY_TASK_ID"p yeast_alignment_samples.inp)

module load python/3.6.1
module load py-numpy/1.18.1_py36
module load py-biopython/1.79_py36

python3 extract_barcodes.py $files