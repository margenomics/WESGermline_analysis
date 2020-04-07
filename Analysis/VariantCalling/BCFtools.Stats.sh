#!/bin/bash
#SBATCH -p normal         # Partition to submit to
#SBATCH --cpus-per-task=1 		  # Number of cpus (cores) per task
#SBATCH --mem-per-cpu 30Gb           # Memory in MB
#SBATCH -J BCFtools.Stats          	# job name
#SBATCH -o logs/BCFtools.Stats.%j.out    # File to which standard out will be written
#SBATCH -e logs/BCFtools.Stats.%j.err    # File to which standard err will be written



module purge  ## Why? Clear out .bashrc /.bash_profile settings that might interfere
module load BCFtools/1.8

SAMPLE=$1
VCF=$2
OUTDIR=$3



bcftools stats --apply-filters "PASS,." $VCF > ${OUTDIR}/${SAMPLE}_recalibrated.PASS.bcftools



