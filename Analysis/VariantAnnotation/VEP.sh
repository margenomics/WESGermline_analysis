#!/bin/bash
#SBATCH -p short,normal           # Partition to submit to (12h)
#SBATCH --cpus-per-task=1 		  # Number of cpus (cores) per task
#SBATCH --mem-per-cpu 30Gb           # Memory in MB
#SBATCH -J VEP         	# job name
#SBATCH -o logs/VEP.%j.out    # File to which standard out will be written
#SBATCH -e logs/VEP.%j.err    # File to which standard err will be written


module pure
module load VEP/98
module load vcfmaf/1.6.16



# grab filename base and create output directory
SAMPLE=$1
VCF=$2
OUTDIR=$3


vcf2maf.pl --input-vcf $VCF --output-maf ${OUTDIR}/${SAMPLE}-PASS.VEP.maf --tumor-id $SAMPLE --ref-fasta /bicoh/MARGenomics/Ref_Genomes_fa/GATK_bundle/hg38/ref/Homo_sapiens_assembly38.fasta --ncbi-build GRCh38 --vep-path /soft/EB_repo/bio/sequence/programs/foss/2016b/VEP/98 --vep-data /bicoh/MARGenomics/Ref_Genomes_fa/vep/

#Source			Version (GRCh38)
#Ensembl DB version	98
#Genome assembly	GRCh38.p12
#GENCODE		31
#RefSeq			"2019-06-28(GCF_000001405.39_GRCh38.p13_genomic.gff)"
#Regulatory build	1
#PolyPhen		2.2.2
#SIFT			5.2.2
#dbSNP			152
#COSMIC			89
#HGMD-PUBLIC		2018.4
#ClinVar		2019-04
#1000 Genomes		Phase 3 (remapped)
#NHLBI-ESP		V2-SSA137 (remapped)
#gnomAD			r2.1, exomes only (remapped)
