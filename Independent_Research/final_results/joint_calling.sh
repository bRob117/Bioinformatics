#!/bin/bash
#SBATCH -N 1
#SBATCH -n 32
#SBATCH --mem=96G
#SBATCH -t 48:00:00
#SBATCH -o logs/slurm/joint_calling_%j.log
#SBATCH -e logs/slurm/joint_calling_%j.err
#SBATCH --account=barobinson4742
#SBATCH --partition=silver
#
# Joint variant calling across all 81 samples.
# Run AFTER all array tasks from tomato_array_slurm.sh have completed.
# Submit with: sbatch joint_calling.sh

set -euo pipefail

module purge
module load biological/samtools_1.23
module load biological/java
module load biological/perl_5.40
module load conda/miniconda3

set +u
source activate bio
set -u

export PROJ_DIR=/export/home/bio_class/barobinson4742/Bioinformatics
cd $PROJ_DIR

THREADS=32
REF=genome/SL4.0.genome.fasta

# ---- Verify all samples completed ----
EXPECTED=$(wc -l < samples_subset.txt)
ACTUAL=$(ls alignment/*.dedup.bam 2>/dev/null | wc -l)
echo "Expected $EXPECTED samples, found $ACTUAL dedup BAMs"

if [ "$ACTUAL" -lt "$EXPECTED" ]; then
    echo "WARNING: Not all samples have completed alignment."
    echo "Missing samples:"
    while read SAMPLE; do
        [ ! -f "alignment/${SAMPLE}.dedup.bam" ] && echo "  $SAMPLE"
    done < samples_subset.txt
    echo "Proceeding with available samples..."
fi

# ---- Joint variant calling with bcftools ----
# mpileup all BAMs together, then call variants jointly.
# This is the correct approach for population-level variant analysis.
echo "Starting joint variant calling across $ACTUAL samples..."

if [ ! -f variants/all_samples.vcf.gz ]; then
    bcftools mpileup -f $REF \
        --threads $THREADS \
        --annotate FORMAT/AD,FORMAT/DP \
        alignment/*.dedup.bam \
      | bcftools call -mv \
        --threads $THREADS \
        -Oz -o variants/all_samples.vcf.gz

    bcftools index variants/all_samples.vcf.gz
fi

# ---- Filter joint VCF ----
if [ ! -f variants/all_samples.filtered.vcf.gz ]; then
    bcftools filter \
        -e 'QUAL<30 || INFO/DP<10' \
        -Oz -o variants/all_samples.filtered.vcf.gz \
        variants/all_samples.vcf.gz

    bcftools index variants/all_samples.filtered.vcf.gz
fi

# ---- Joint VCF stats ----
bcftools stats variants/all_samples.filtered.vcf.gz > logs/joint_vcf_stats.txt

# ---- Summary ----
VARIANTS=$(zcat variants/all_samples.filtered.vcf.gz | grep -v "^#" | wc -l)
SNPS=$(bcftools stats variants/all_samples.filtered.vcf.gz | grep "^SN" | grep "number of SNPs" | awk '{print $NF}')
INDELS=$(bcftools stats variants/all_samples.filtered.vcf.gz | grep "^SN" | grep "number of indels" | awk '{print $NF}')
TSTV=$(bcftools stats variants/all_samples.filtered.vcf.gz | grep "^TSTV" | awk '{print $5}')

echo ""
echo "========== JOINT CALLING COMPLETE =========="
echo "Samples:  $ACTUAL"
echo "Variants: $VARIANTS"
echo "SNPs:     $SNPS"
echo "Indels:   $INDELS"
echo "Ts/Tv:    $TSTV"
echo "Output:   variants/all_samples.filtered.vcf.gz"
echo "============================================"
