#!/bin/bash
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -o lab9.log
#SBATCH --account=barobinson4742
#SBATCH --partition=silver

#make directories to sort outputs
mkdir -p genome fastq alignment variants

#purge modules to make a clear workspace
module purge

#load samtools and java(for picard)
module load biological/samtools_1.23
module load biological/java
module load biological/perl_5.40

#make my account directory a variable as PROJ_DIR
export PROJ_DIR=/export/home/bio_class/barobinson4742/Bioinformatics/
cd $PROJ_DIR
export SRR=SRR5324768

#copy reference genome in genome/ directory
cp $PROJ_DIR/GCF_900604845.1_TTHNAR1_genomic.fna genome/

#use picard to make reference genome .dict file
java -jar /export/share/software/biological/picard/picard.jar CreateSequenceDictionary \
R=genome/GCF_900604845.1_TTHNAR1_genomic.fna \
O=genome/GCF_900604845.1_TTHNAR1_genomic.dict

#run bowtie2 build
/export/share/software/biological/bowtie2-2.4.2-sra-linux-x86_64/bowtie2-build \
genome/GCF_900604845.1_TTHNAR1_genomic.fna \
genome/GCF_900604845.1_TTHNAR1_genomic

#run bowtie2
/export/share/software/biological/bowtie2-2.4.2-sra-linux-x86_64/bowtie2 -x \
genome/GCF_900604845.1_TTHNAR1_genomic \
-1 fastq/${SRR}_pass_1.fastq.gz \
-2 fastq/${SRR}_pass_2.fastq.gz --sensitive-local \
--rg-id ${SRR} --rg SM:${SRR} --rg PL:ILLUMINA \
> alignment/${SRR}.sam

#run samtools
samtools view -hb alignment/${SRR}.sam | samtools sort -l 5 -o alignment/${SRR}.bam

#download consensus
samtools consensus -f fasta -o ${SRR}_consensus.fasta alignment/${SRR}.bam
