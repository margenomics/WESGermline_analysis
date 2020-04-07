#!/bin/bash
#SBATCH -p normal,long           	     # Partition to submit to (12h)
#SBATCH --cpus-per-task=1            # Number of cores (chosen from scale test)
#SBATCH --mem-per-cpu 91Gb           # Memory in MB
#SBATCH -J JointCaller.GVCF           	# job name
#SBATCH -o logs/JointCaller.GVCF.%j.out    # File to which standard out will be written
#SBATCH -e logs/JointCaller.GVCF.%j.err    # File to which standard err will be written



module purge  ## Why? Clear out .bashrc /.bash_profile settings that might interfere
module load GATK/4.1.4.0



# Same REF used in the alignment
REF='/bicoh/MARGenomics/Ref_Genomes_fa/GATK_bundle/hg38/ref/Homo_sapiens_assembly38.fasta'
INTERVALS='/bicoh/MARGenomics/annotationData/ExomeTargetRegions/SureSelect.Human.AllExon.V6.60Mb_hs_hg38/SureSelect.Human.AllExon.V6.60Mb_hs_hg38.interval_list'

DIR=$1

gatk --java-options '-Xmx90g' CombineGVCFs \
	-R $REF \
	-L $INTERVALS \
	-V ${DIR}/VHIJA_01064AAB_TCCTCT.g.vcf \
	-V ${DIR}/VMADRE_01065AAB_GTTCAG.g.vcf \
	-V ${DIR}/VTIO_01062AAB_ATTCTC.g.vcf \
	-V ${DIR}/VHIJO_01063AAB_CAAGGA.g.vcf \
	-V ${DIR}/VPADRE_01061AAB_TAGTCG.g.vcf \
	-O ${DIR}/Combined.g.vcf
	
gatk --java-options '-Xmx90g' GenotypeGVCFs \
	-R $REF \
	-L $INTERVALS \
	-V ${DIR}/Combined.g.vcf \
	-O ${DIR}/Jointcalls.vcf





