#!/bin/bash
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --mem=48G
#SBATCH -t 12:00:00
#SBATCH -o logs/slurm/tomato_%A_%a.log
#SBATCH -e logs/slurm/tomato_%A_%a.err
#SBATCH --account=barobinson4742
#SBATCH --partition=silver
#SBATCH --array=1-15%10
#
# Job array: runs 81 tasks (one per sample), max 10 concurrent.
# Each task reads its sample ID from samples.txt using SLURM_ARRAY_TASK_ID.
# Adjust %10 to run more/fewer simultaneously (max depends on cluster load).

set -euo pipefail

# ---- Modules ----
module purge
module load biological/samtools_1.23
module load biological/java
module load biological/perl_5.40
module load conda/miniconda3

set +u
source activate bio   # provides bcftools, bwa, fastqc, fastp
set -u

# ---- Paths ----
PICARD=/export/share/software/biological/picard/picard.jar
export PROJ_DIR=/export/home/bio_class/barobinson4742/Bioinformatics
cd $PROJ_DIR

THREADS=16
REF=genome/SL4.0.genome.fasta

# ---- Get sample ID for this array task ----
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" samples_subset.txt)
if [ -z "$SAMPLE" ]; then
    echo "ERROR: No sample found for array task $SLURM_ARRAY_TASK_ID"
    exit 1
fi
echo "===== Processing $SAMPLE (task $SLURM_ARRAY_TASK_ID of $SLURM_ARRAY_TASK_COUNT) ====="

# ---- Verify input files exist ----
if [ ! -f sra/${SAMPLE}_1.fastq.gz ] || [ ! -f sra/${SAMPLE}_2.fastq.gz ]; then
    echo "ERROR: FASTQ files not found for $SAMPLE in sra/"
    exit 1
fi

# ---- Create per-sample output dirs ----
mkdir -p trimmed alignment variants logs

# ---- Reference prep (one-time, idempotent) ----
# Multiple array tasks may hit this simultaneously on first run.
# File existence checks prevent redundant work; flock prevents corruption.
if [ ! -f ${REF}.bwt ]; then
    (
        flock -n 200 || { echo "BWA index being built by another task, waiting..."; flock 200; }
        [ ! -f ${REF}.bwt ] && bwa index $REF
    ) 200>/tmp/bwa_index.lock
fi
if [ ! -f genome/SL4.0.genome.dict ]; then
    (
        flock -n 201 || { echo "Dict being built by another task, waiting..."; flock 201; }
        [ ! -f genome/SL4.0.genome.dict ] && java -jar $PICARD CreateSequenceDictionary R=$REF O=genome/SL4.0.genome.dict
    ) 201>/tmp/picard_dict.lock
fi
[ ! -f ${REF}.fai ] && samtools faidx $REF

# ---- 1. FastQC on raw reads ----
if [ ! -f logs/${SAMPLE}_1_fastqc.html ]; then
    fastqc -t $THREADS -o logs/ sra/${SAMPLE}_1.fastq.gz sra/${SAMPLE}_2.fastq.gz
fi

# ---- 2. Trim adapters and filter reads with >10% Ns (per Kim & Lee 2024) ----
if [ ! -f trimmed/${SAMPLE}_1.fastq.gz ]; then
    fastp -i sra/${SAMPLE}_1.fastq.gz -I sra/${SAMPLE}_2.fastq.gz \
        -o trimmed/${SAMPLE}_1.fastq.gz -O trimmed/${SAMPLE}_2.fastq.gz \
        -n 10 \
        --thread $THREADS \
        -j logs/${SAMPLE}.fastp.json \
        -h logs/${SAMPLE}.fastp.html
fi

# ---- 3. Align with BWA-MEM (per Kim & Lee 2024: BWA 0.7.17) ----
if [ ! -f alignment/${SAMPLE}.sorted.bam ]; then
    bwa mem -t $THREADS \
        -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tLB:${SAMPLE}\tPL:ILLUMINA" \
        $REF \
        trimmed/${SAMPLE}_1.fastq.gz trimmed/${SAMPLE}_2.fastq.gz \
        2> logs/${SAMPLE}.bwa.log \
      | samtools sort -@ $THREADS -m 2G -o alignment/${SAMPLE}.sorted.bam -
    samtools index alignment/${SAMPLE}.sorted.bam
fi

# ---- 4. Mark duplicates ----
if [ ! -f alignment/${SAMPLE}.dedup.bam ]; then
    java -Xmx16G -jar $PICARD MarkDuplicates \
        I=alignment/${SAMPLE}.sorted.bam \
        O=alignment/${SAMPLE}.dedup.bam \
        M=logs/${SAMPLE}.dup_metrics.txt
    samtools index alignment/${SAMPLE}.dedup.bam
fi

# ---- 5. Per-sample variant calling ----
if [ ! -f variants/${SAMPLE}.vcf.gz ]; then
    bcftools mpileup -f $REF --threads $THREADS alignment/${SAMPLE}.dedup.bam \
      | bcftools call -mv -Oz --threads $THREADS -o variants/${SAMPLE}.vcf.gz
    bcftools index variants/${SAMPLE}.vcf.gz
fi

# ---- 6. QC ----
samtools flagstat alignment/${SAMPLE}.dedup.bam > logs/${SAMPLE}.flagstat.txt
samtools depth -a alignment/${SAMPLE}.dedup.bam \
    | awk '{sum+=$3; n++} END {print "Mean coverage:", sum/n}' \
    > logs/${SAMPLE}.coverage.txt
bcftools stats variants/${SAMPLE}.vcf.gz > logs/${SAMPLE}.vcf_stats.txt

# ---- 7. Cleanup sorted BAM (redundant after dedup) ----
if [ -f alignment/${SAMPLE}.dedup.bam ] && [ -f alignment/${SAMPLE}.sorted.bam ]; then
    rm alignment/${SAMPLE}.sorted.bam alignment/${SAMPLE}.sorted.bam.bai
    echo "Cleaned up sorted BAM to save disk space"
fi

echo "===== Done: $SAMPLE ====="
