#!/usr/bin/env bash
#SBATCH -J WGS_Alignment
#SBATCH -p hns,dpetrov,normal,owners
#SBATCH -n 16
#SBATCH -N 1
#SBATCH -t 2-00:00
#SBATCH --array=1-1
#SBATCH --mem-per-cpu=2G
#SBATCH --requeue
#SBATCH -o SlurmFiles/slurm-%A_%a_%x.out
#SBATCH --mail-user=grantkinsler@gmail.com
#SBATCH --mail-type=END

bn=$(sed -n "$SLURM_ARRAY_TASK_ID"p EvolvingFront_WGS_samples.inp)

ml system
ml ncurses/6.0
ml biology
ml bedtools
ml java


ref_genome=$"/home/groups/dpetrov/grant/S288C_reference_genome_R64-1-1_20110203/S288C_reference_sequence_R64-1-1_20110203.fasta"

bam_dir=$"Output/"
out_dir=$"smoove_out/"

smoove genotype -d -x -p 1 --name "$bn"-joint --outdir results-genotped/ --fasta "$ref_genome" --vcf merged.sites.vcf.gz "$bam_dir$bn".sorted.rg.md.bam
