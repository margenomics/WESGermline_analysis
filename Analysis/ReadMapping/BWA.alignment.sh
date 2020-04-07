#!/bin/bash
#SBATCH -p short             # Partition to submit to
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu 11Gb     # Memory in MB
#SBATCH -J BWA           	# job name
#SBATCH -o logs/BWA.%j.out    # File to which standard out will be written
#SBATCH -e logs/BWA.%j.err    # File to which standard err will be written



module purge  ## Why? Clear out .bashrc /.bash_profile settings that might interfere
module load BWA/0.7.17
module load SAMtools/1.8-foss-2016b

# grab filename base and create output directory
name=$1
indir=$2
outdir=$3

# For read group
SAMPLE=$(echo $1 | cut -d"_" -f1)
RGLB=$(echo $1 | cut -d"_" -f2)
RGPU=$(echo $1 | cut -d"_" -f2,3)

REF='/bicoh/MARGenomics/Ref_Genomes_fa/GATK_bundle/hg38/ref/Homo_sapiens_assembly38.fasta'

# Change prefix
R1=_read1.fastq.gz
R2=_read2.fastq.gz

bwa mem -t $SLURM_CPUS_PER_TASK -R $(echo "@RG\tID:${SAMPLE}_${RGPU}\tLB:$RGLB\tSM:$SAMPLE\tPL:ILLUMINA\tPU:$RGPU") $REF ${indir}/${name}${R1} ${indir}/${name}${R2} | samtools sort -o ${outdir}/${name}-pe_sorted.bam

# Afegeixo aquí el RG per heredar-lo si hem d'ajuntar bams.
# Add RG -R @RG\tID:${1}\tLB:Exome1\tSM:Exome1\t:PL:ILLUMINA
# Create index file (we do it later in MarkDuplicates)
# samtools index 04_MAP/${1}-pe_sorted.bam


# Validate BAM (in stdout):
module load picard/2.2.4-Java-1.8.0_92


java -jar $EBROOTPICARD/picard.jar  ValidateSamFile \
	I=${outdir}/${name}-pe_sorted.bam \
	MODE=SUMMARY
