#!/usr/bin/env bash
#SBATCH -J Smoove1
#SBATCH -p hns,dpetrov,normal,owners
#SBATCH -n 16
#SBATCH -N 1
#SBATCH -t 2-00:00
#SBATCH --array=1-457
#SBATCH --mem-per-cpu=2G
#SBATCH --requeue
#SBATCH -o SlurmFiles/slurm-%A_%a_%x.out
#SBATCH --mail-user=grantkinsler@gmail.com
#SBATCH --mail-type=END

bn=$(sed -n "$SLURM_ARRAY_TASK_ID"p EvolvingFront_WGS_samples.inp)

ref_genome=$"/home/groups/dpetrov/grant/S288C_reference_genome_R64-1-1_20110203/S288C_reference_sequence_R64-1-1_20110203.fasta"

bam_dir=$"Output/"
out_dir=$"smoove_out/"

singularity exec /home/groups/dpetrov/SOFTWARE/simg/smoove_latest.sif smoove call --outdir "$out_dir" --name "$bn" --fasta "$ref_genome" -p 1 --genotype "$bam_dir$bn".sorted.rg.md.bam



