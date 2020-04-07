#!/bin/bash
#SBATCH -p long            # Partition to submit to
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu 4Gb     # Memory in MB
#SBATCH -J FastQC           # job name
#SBATCH -o logs/FastQC.%j.out    # File to which standard out will be written
#SBATCH -e logs/FastQC.%j.err    # File to which standard err will be written



module purge 
module load FastQC/0.11.5-Java-1.7.0_80        

#------------------------
INFILE=$1
OUTDIR=$2


#------------------------

fastqc --outdir $OUTDIR --threads $SLURM_CPUS_PER_TASK $INFILE

