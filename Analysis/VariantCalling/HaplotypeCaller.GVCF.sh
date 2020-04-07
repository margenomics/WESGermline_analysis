#!/bin/bash
#SBATCH -p normal,long           	     # Partition to submit to (12h)
#SBATCH --cpus-per-task=1            # Number of cores (chosen from scale test)
#SBATCH --mem-per-cpu 91Gb           # Memory in MB
#SBATCH -J HaplotypeCaller.GVCF           	# job name
#SBATCH -o logs/HaplotypeCaller.GVCF.%j.out    # File to which standard out will be written
#SBATCH -e logs/HaplotypeCaller.GVCF.%j.err    # File to which standard err will be written



module purge  ## Why? Clear out .bashrc /.bash_profile settings that might interfere
module load GATK/4.1.4.0



# Same REF used in the alignment
REF='/bicoh/MARGenomics/Ref_Genomes_fa/GATK_bundle/hg38/ref/Homo_sapiens_assembly38.fasta'
INTERVALS='/bicoh/MARGenomics/annotationData/ExomeTargetRegions/SureSelect.Human.AllExon.V6.60Mb_hs_hg38/SureSelect.Human.AllExon.V6.60Mb_hs_hg38.interval_list'

SAMPLE=$1
INFILE=$2
DIR=$3

gatk --java-options '-Xmx60g' HaplotypeCaller \
	-ERC GVCF \
	-R $REF \
	-L $INTERVALS \
	-I $INFILE \
	-O ${DIR}/${SAMPLE}.g.vcf \
	--bam-output ${DIR}/${SAMPLE}.bam 





