#!/bin/bash
#SBATCH -p normal,long           	     # Partition to submit to (12h)
#SBATCH --cpus-per-task=1            # Number of cores (chosen from scale test)
#SBATCH --mem-per-cpu 91Gb           # Memory in MB
#SBATCH -J VQSR.filtering           	# job name
#SBATCH -o logs/VQSR.filtering.%j.out    # File to which standard out will be written
#SBATCH -e logs/VQSR.filtering.%j.err    # File to which standard err will be written



module purge  ## Why? Clear out .bashrc /.bash_profile settings that might interfere
module load GATK/4.1.4.0
module load R


# Same REF used in the alignment
REF='/bicoh/MARGenomics/Ref_Genomes_fa/GATK_bundle/hg38/ref/Homo_sapiens_assembly38.fasta'
HAPMAP='/bicoh/MARGenomics/Ref_Genomes_fa/GATK_bundle/hg38/hapmap_3.3.hg38.vcf.gz'
OMNI='/bicoh/MARGenomics/Ref_Genomes_fa/GATK_bundle/hg38/1000G_omni2.5.hg38.vcf.gz'
K1000G='/bicoh/MARGenomics/Ref_Genomes_fa/GATK_bundle/hg38/1000G_phase1.snps.high_confidence.hg38.vcf'
DBSNP='/bicoh/MARGenomics/Ref_Genomes_fa/GATK_bundle/hg38/dbsnp_146.hg38.vcf'
MILLS='/bicoh/MARGenomics/Ref_Genomes_fa/GATK_bundle/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf'



DIR=$1


gatk --java-options '-Xmx90g' VariantRecalibrator \
	-R $REF \
	-V ${DIR}/Jointcalls.vcf \
	--resource:hapmap,known=false,training=true,truth=true,prior=15.0 $HAPMAP \
	--resource:omni,known=false,training=true,truth=false,prior=12.0 $OMNI \
	--resource:1000G,known=false,training=true,truth=false,prior=10.0 $K1000G \
	--resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $DBSNP \
	-an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
	-mode SNP \
	-tranche 100.0 \
	-tranche 99.9 \
	-tranche 99.0 \
	-tranche 90.0 \
	-O ${DIR}/recalibrate_SNP.recal \
	--tranches-file ${DIR}/recalibrate_SNP.tranches \
	--rscript-file ${DIR}/recalibrate_SNP.plots.R
	

gatk --java-options '-Xmx90g' ApplyVQSR \
	-R $REF \
	-V ${DIR}/Jointcalls.vcf \
	--mode SNP \
	--truth-sensitivity-filter-level 99.0 \
	--recal-file ${DIR}/recalibrate_SNP.recal \
	--tranches-file ${DIR}/recalibrate_SNP.tranches \
	-O ${DIR}/Jointcalls_recalibrated_snps_raw_indels.vcf



gatk --java-options '-Xmx90g' VariantRecalibrator \
	-R $REF \
	-V ${DIR}/Jointcalls_recalibrated_snps_raw_indels.vcf \
	--resource:mills,known=false,training=true,truth=true,prior=15.0 $MILLS \
	--resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $DBSNP \
	-an QD -an FS -an SOR -an MQRankSum -an ReadPosRankSum \
	-mode INDEL \
	-tranche 100.0 \
	-tranche 99.9 \
	-tranche 99.0 \
	-tranche 90.0 \
	--max-gaussians 4 \
	--tranches-file ${DIR}/recalibrate_INDEL.tranches \
	--rscript-file ${DIR}/recalibrate_INDEL.plots.R \
	-O ${DIR}/recalibrate_INDEL.recal 



gatk --java-options '-Xmx90g' ApplyVQSR \
	-R $REF \
	-V ${DIR}/Jointcalls_recalibrated_snps_raw_indels.vcf \
	--mode INDEL \
	--truth-sensitivity-filter-level 99.0 \
	--recal-file ${DIR}/recalibrate_INDEL.recal \
	--tranches-file ${DIR}/recalibrate_INDEL.tranches \
	-O ${DIR}/Jointcalls_recalibrated.vcf


#### PASS Joint
gatk --java-options '-Xmx90g' SelectVariants \
	-R $REF \
	-V ${DIR}/Jointcalls_recalibrated.vcf \
	-O ${DIR}/Jointcalls_recalibrated.PASS.vcf \
	--exclude-filtered \
	--exclude-non-variants \
	--remove-unused-alternates 
	
#### PASS sample-wise
for sample in $(echo VHIJA VHIJO VMADRE VPADRE VTIO) 
	do
	echo $sample
	gatk --java-options '-Xmx90g' SelectVariants \
		-R $REF \
		-V ${DIR}/Jointcalls_recalibrated.vcf \
		-O ${DIR}/${sample}_recalibrated.PASS.vcf \
		--exclude-filtered \
		--exclude-non-variants \
		--remove-unused-alternates \
		--sample-name $sample
done


#### ALL (not filtered but with PASS flag) sample-wise
for sample in $(echo VHIJA VHIJO VMADRE VPADRE VTIO) 
	do
	echo $sample
	gatk --java-options '-Xmx90g' SelectVariants \
		-R $REF \
		-V ${DIR}/Jointcalls_recalibrated.vcf \
		-O ${DIR}/${sample}_recalibrated.vcf \
		--exclude-non-variants \
		--remove-unused-alternates \
		--sample-name $sample
done

