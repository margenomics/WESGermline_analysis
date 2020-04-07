#!/bin/bash
#SBATCH -p normal         # Partition to submit to
#SBATCH --cpus-per-task=1 		  # Number of cpus (cores) per task
#SBATCH --mem-per-cpu 30Gb           # Memory in MB
#SBATCH -J annovar          	# job name
#SBATCH -o logs/annovar.%j.out    # File to which standard out will be written
#SBATCH -e logs/annovar.%j.err    # File to which standard err will be written



module purge  ## Why? Clear out .bashrc /.bash_profile settings that might interfere


SAMPLE=$1
VCF=$2
OUTDIR=$3
built=$4



/bicoh/MARGenomics/Software/annovar/convert2annovar.pl -format vcf4 $VCF > ${OUTDIR}/${SAMPLE}-PASS.avinput


/bicoh/MARGenomics/Software/annovar/table_annovar.pl ${OUTDIR}/${SAMPLE}-PASS.avinput /bicoh/MARGenomics/Software/annovar/humandb -buildver $built -out ${OUTDIR}/${SAMPLE}-PASS.annovar -remove -protocol  refGene,knownGene,avsnp150,dbnsfp33a,clinvar_20190305,cosmic70,gnomad211_exome -operation g,g,f,f,f,f,f -nastring . -csvout 

rm ${OUTDIR}/${SAMPLE}-PASS.avinput
