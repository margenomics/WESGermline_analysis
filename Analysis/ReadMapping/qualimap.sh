#!/bin/bash
#SBATCH -p normal,long           # Partition to submit to (12h)
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu 30Gb     # Memory in MB
#SBATCH -J qualimap           	# job name
#SBATCH -o logs/qualimap.%j.out    # File to which standard out will be written
#SBATCH -e logs/qualimap.%j.err    # File to which standard err will be written



# grab filename base and create output directory

SAMPLE=$1
INFILE=$2
OUDIR=$3


TARGETS=/bicoh/MARGenomics/annotationData/ExomeTargetRegions/SureSelect.Human.AllExon.V6.60Mb_hs_hg38/S07604514_Covered_6fields.bed


cd $OUTDIR

/bicoh/MARGenomics/Software/qualimap_v2.2.1/qualimap bamqc -bam $INFILE -gff $TARGETS --java-mem-size=30G -sd


# -sd : skip duplicates marked with flag
