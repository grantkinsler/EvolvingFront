#!/usr/bin/env bash
#SBATCH -J WGS_Alignment
#SBATCH -p hns,dpetrov,normal,owners
#SBATCH -n 16
#SBATCH -N 1
#SBATCH -t 2-00:00
#SBATCH --array=1-96
#SBATCH --mem-per-cpu=2G
#SBATCH --requeue
#SBATCH -o SlurmFiles/slurm-%A_%a_%x.out
#SBATCH --mail-user=grantkinsler@gmail.com
#SBATCH --mail-type=END

bn=$(sed -n "$SLURM_ARRAY_TASK_ID"p EvolvingFront_WGS_Plate5_samples.inp)

in1="$(cut -d',' -f1 <<<"$bn")"
in2="$(cut -d',' -f2 <<<"$bn")"
in3="$(cut -d',' -f3 <<<"$bn")"
in4="$(cut -d',' -f4 <<<"$bn")"
out="$(cut -d',' -f5 <<<"$bn")"

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

# bwa mem -t 16 "$ref_genome" "$dir$in1"_R1.fastq.gz "$dir$in1"_R2.fastq.gz | samtools view -bS - | samtools sort > "$out_dir$in1".sorted.bam
# bwa mem -t 16 "$ref_genome" "$dir$in2"_R1.fastq.gz "$dir$in2"_R2.fastq.gz | samtools view -bS - | samtools sort > "$out_dir$in2".sorted.bam
# bwa mem -t 16 "$ref_genome" "$dir$in3"_R1.fastq.gz "$dir$in3"_R2.fastq.gz | samtools view -bS - | samtools sort > "$out_dir$in3".sorted.bam
# bwa mem -t 16 "$ref_genome" "$dir$in4"_R1.fastq.gz "$dir$in4"_R2.fastq.gz | samtools view -bS - | samtools sort > "$out_dir$in4".sorted.bam

# ### fix read groups
java -jar /home/groups/dpetrov/SOFTWARE/picard-tools/picard.jar AddOrReplaceReadGroups I="$out_dir$in1".sorted.bam O="$out_dir$in1".sorted.rg.bam RGID="$in1" RGLB=lib"$out" RGPL=illumina RGPU=unit"$in1" RGSM="$out"
samtools index "$out_dir$in1".sorted.rg.bam
java -jar /home/groups/dpetrov/SOFTWARE/picard-tools/picard.jar AddOrReplaceReadGroups I="$out_dir$in2".sorted.bam O="$out_dir$in2".sorted.rg.bam RGID="$in2" RGLB=lib"$out" RGPL=illumina RGPU=unit"$in2" RGSM="$out"
samtools index "$out_dir$in2".sorted.rg.bam
java -jar /home/groups/dpetrov/SOFTWARE/picard-tools/picard.jar AddOrReplaceReadGroups I="$out_dir$in3".sorted.bam O="$out_dir$in3".sorted.rg.bam RGID="$in3" RGLB=lib"$out" RGPL=illumina RGPU=unit"$in3" RGSM="$out"
samtools index "$out_dir$in3".sorted.rg.bam
java -jar /home/groups/dpetrov/SOFTWARE/picard-tools/picard.jar AddOrReplaceReadGroups I="$out_dir$in4".sorted.bam O="$out_dir$in4".sorted.rg.bam RGID="$in4" RGLB=lib"$out" RGPL=illumina RGPU=unit"$in4" RGSM="$out"
samtools index "$out_dir$in4".sorted.rg.bam

# samtools merge "$out_dir$out".sorted.rg.bam "$out_dir$in1".sorted.rg.bam "$out_dir$in2".sorted.rg.bam "$out_dir$in3".sorted.rg.bam "$out_dir$in4".sorted.rg.bam
### mark duplicates
java -jar /home/groups/dpetrov/SOFTWARE/picard-tools/picard.jar MarkDuplicates I="$out_dir$in1".sorted.bam I="$out_dir$in2".sorted.bam I="$out_dir$in3".sorted.bam I="$out_dir$in4".sorted.bam O="$out_dir$out".sorted.md.bam M="$out_dir$out".md.metrics.txt
samtools index "$out_dir$out".sorted.md.bam

java -jar /home/groups/dpetrov/SOFTWARE/picard-tools/picard.jar AddOrReplaceReadGroups I="$out_dir$out".sorted.md.bam O="$out_dir$out".sorted.rg.md.bam RGID="$out" RGLB=lib"$out" RGPL=illumina RGPU=unit"$out" RGSM="$out"
samtools index "$out_dir$out".sorted.rg.md.bam

### recalibrate base calling [we don't have a known VCF!!!]
# /home/groups/dpetrov/SOFTWARE/gatk-4.2.0.0/gatk BaseRecalibrator -R "$ref_genome" -I "$out_dir$bn".sorted.rg.md.bam -O "$out_dir$bn".recal_data.table
# /home/groups/dpetrov/SOFTWARE/gatk-4.2.0.0/gatk ApplyBQSR -R "$ref_genome" -I "$out_dir$bn".sorted.rg.md.bam --bqsr-recal-file "$out_dir$bn".recal_data.table -O "$out_dir$bn".sorted.rg.md.br.bam

/home/groups/dpetrov/SOFTWARE/gatk-4.2.0.0/gatk --java-options "-Xmx4g" HaplotypeCaller -R "$ref_genome" -I "$out_dir$out".sorted.rg.md.bam -ERC GVCF -O "$out_dir$out".g.vcf.gz
