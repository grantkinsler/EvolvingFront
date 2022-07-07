#!/usr/bin/env bash
#SBATCH -J WGS_Alignment
#SBATCH -p hns,dpetrov,normal,owners
#SBATCH -n 16
#SBATCH -N 1
#SBATCH -t 2-00:00
#SBATCH --array=170-361
#SBATCH --mem-per-cpu=2G
#SBATCH --requeue
#SBATCH -o SlurmFiles/slurm-%A_%a_%x.out
#SBATCH --mail-user=grantkinsler@gmail.com
#SBATCH --mail-type=END

bn=$(sed -n "$SLURM_ARRAY_TASK_ID"p EvolvingFront_WGS_samples.inp)

ml system
ml ncurses/6.0
ml biology
ml bwa/0.7.17
ml java

ref_genome=$"/home/groups/dpetrov/grant/S288C_reference_genome_R64-1-1_20110203/S288C_reference_sequence_R64-1-1_20110203.fasta"
dir=$"AllRawData/"
out_dir=$"Output/"



# bwa index "$ref_genome"
# /home/groups/dpetrov/SOFTWARE/gatk-4.2.0.0/gatk CreateSequenceDictionary -R "$ref_genome"
# samtools faidx "$ref_genome"

bwa mem -t 16 "$ref_genome" "$dir$bn"_R1.fastq.gz "$dir$bn"_R2.fastq.gz | samtools view -bS - | samtools sort > "$out_dir$bn".sorted.bam

### fix read groups
java -jar /home/groups/dpetrov/SOFTWARE/picard-tools/picard.jar AddOrReplaceReadGroups I="$out_dir$bn".sorted.bam O="$out_dir$bn".sorted.rg.bam RGID="$bn" RGLB=lib"$bn" RGPL=illumina RGPU=unit"$bn" RGSM="$bn"
samtools index "$out_dir$bn".sorted.rg.bam

### mark duplicates
java -jar /home/groups/dpetrov/SOFTWARE/picard-tools/picard.jar MarkDuplicates I="$out_dir$bn".sorted.rg.bam O="$out_dir$bn".sorted.rg.md.bam M="$out_dir$bn".md.metrics.txt
samtools index "$out_dir$bn".sorted.rg.md.bam

### recalibrate base calling [we don't have a known VCF!!!]
# /home/groups/dpetrov/SOFTWARE/gatk-4.2.0.0/gatk BaseRecalibrator -R "$ref_genome" -I "$out_dir$bn".sorted.rg.md.bam -O "$out_dir$bn".recal_data.table
# /home/groups/dpetrov/SOFTWARE/gatk-4.2.0.0/gatk ApplyBQSR -R "$ref_genome" -I "$out_dir$bn".sorted.rg.md.bam --bqsr-recal-file "$out_dir$bn".recal_data.table -O "$out_dir$bn".sorted.rg.md.br.bam

/home/groups/dpetrov/SOFTWARE/gatk-4.2.0.0/gatk --java-options "-Xmx4g" HaplotypeCaller -R "$ref_genome" -I "$out_dir$bn".sorted.rg.md.bam -ERC GVCF -O "$out_dir$bn".g.vcf.gz
