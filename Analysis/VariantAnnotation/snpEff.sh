#!/bin/bash
#SBATCH -p short,normal           # Partition to submit to (12h)
#SBATCH --cpus-per-task=1 		  # Number of cpus (cores) per task
#SBATCH --mem-per-cpu 30Gb           # Memory in MB
#SBATCH -J snpEff          	# job name
#SBATCH -o logs/snpEff.%j.out    # File to which standard out will be written
#SBATCH -e logs/snpEff.%j.err    # File to which standard err will be written



module purge  ## Why? Clear out .bashrc /.bash_profile settings that might interfere
module load snpEff/20180918


# grab filename base and create output directory
SAMPLE=$1
VCF=$2
OUTDIR=$3


snpEff -dataDir /bicoh/MARGenomics/Ref_Genomes_fa/snpEff -csvStats ${OUTDIR}/${SAMPLE}.snpEff -v hg38 $VCF > ${OUTDIR}/${SAMPLE}-PASS.snpEFF.vcf


