#!/usr/bin/env bash
#SBATCH -J EvolutionTracking
#SBATCH -p owners,normal,hns,dpetrov
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --mem=10000
#SBATCH -t 5-00:00
#SBATCH --array=1-5
#SBATCH --requeue
#SBATCH -o SlurmFiles/slurm-%A_%a_%x.out
#SBATCH --mail-user=grantkinsler@gmail.com
#SBATCH --mail-type=END

parameters=$(sed -n "$SLURM_ARRAY_TASK_ID"p evolutiontracking_lane2.inp)

module load python/2.7.1

python bartender_BC1_BC2_bfa.py $parameters

