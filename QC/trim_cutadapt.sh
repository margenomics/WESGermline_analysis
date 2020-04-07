#!/bin/bash
#SBATCH -p short            # Partition to submit to
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu 24Gb     # Memory in MB
#SBATCH -J Cutadapt           # job name
#SBATCH -o logs/Cutadapt.%j.out    # File to which standard out will be written
#SBATCH -e logs/Cutadapt.%j.err    # File to which standard err will be written

#-------------------------------

module purge  ## Why? Clear out .bashrc /.bash_profile settings that might interfere
module load Python/3.5.2-foss-2016b

#-------------------------------
# grab filename base and create output directory
name=$1
indir=$2
outdir=$3

#-------------------------------


cutadapt -j $SLURM_CPUS_PER_TASK -m 20 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT  -o ${outdir}/${name}_read1.fastq.gz -p ${outdir}/${name}_read2.fastq.gz ${indir}/${name}_read1.fastq.gz ${indir}/${name}_read2.fastq.gz 


# Adapters (R1 and R2) recommended by cutadapt and Illumina in https://cutadapt.readthedocs.io/en/stable/guide.html#illumina-truseq and https://support.illumina.com/downloads/illumina-adapter-sequences-document-1000000002694.html
# If not gone, we can still try with Universal adapter: 5' AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT

module load FastQC/0.11.5-Java-1.7.0_80

fastqc --outdir $outdir --threads $SLURM_CPUS_PER_TASK ${outdir}/${name}_read1.fastq.gz 

fastqc --outdir $outdir --threads $SLURM_CPUS_PER_TASK ${outdir}/${name}_read2.fastq.gz 

