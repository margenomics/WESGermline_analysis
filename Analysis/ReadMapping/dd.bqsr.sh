#!/bin/bash
#SBATCH -p normal,long           	     # Partition to submit to (12h)
#SBATCH --cpus-per-task=1            # Number of cores (chosen from scale test)
#SBATCH --mem-per-cpu 91Gb           # Memory in MB
#SBATCH -J dd.bqsr           	# job name
#SBATCH -o logs/dd.bqsr.%j.out    # File to which standard out will be written
#SBATCH -e logs/dd.bqsr.%j.err    # File to which standard err will be written



module purge  ## Why? Clear out .bashrc /.bash_profile settings that might interfere
module load picard/2.2.4-Java-1.8.0_92
module load GATK/4.0.8.1



# Same REF used in the alignment
REF='/bicoh/MARGenomics/Ref_Genomes_fa/GATK_bundle/hg38/ref/Homo_sapiens_assembly38.fasta'
SNP1='/bicoh/MARGenomics/Ref_Genomes_fa/GATK_bundle/hg38/dbsnp_146.hg38.vcf'
SNP2='/bicoh/MARGenomics/Ref_Genomes_fa/GATK_bundle/hg38/1000G_phase1.snps.high_confidence.hg38.vcf'
INDELS='/bicoh/MARGenomics/Ref_Genomes_fa/GATK_bundle/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf'

SAMPLE=$1
INFILE=$2
DIR=$3

# Remove duplicates
java -jar $EBROOTPICARD/picard.jar MarkDuplicates \
	INPUT=$INFILE \
	OUTPUT=${DIR}/${SAMPLE}_merged_dedup.bam \
	M=${DIR}/${SAMPLE}.dedup.txt \
	REMOVE_DUPLICATES=FALSE \
	CREATE_INDEX=True


# Build base recalibration model:
gatk BaseRecalibrator \
	-R $REF \
	-I ${DIR}/${SAMPLE}_merged_dedup.bam \
	--known-sites $SNP1 \
	--known-sites $SNP2\
	--known-sites $INDELS \
	-O ${DIR}/${SAMPLE}_recal.table



# Recalibrate scores based on recal.table:
gatk ApplyBQSR \
	-R $REF \
	-I ${DIR}/${SAMPLE}_merged_dedup.bam \
	-bqsr ${DIR}/${SAMPLE}_recal.table \
	-O ${DIR}/${SAMPLE}_merged_dedup_bqsr.bam

