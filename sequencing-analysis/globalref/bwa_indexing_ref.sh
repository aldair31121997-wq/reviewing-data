#!/bin/bash
#SBATCH -p workq
#SBATCH -J ref_indexing
#SBATCH -e error_ref_indexing
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mem=10G

#modules
module purge
module load bioinfo/bwa/0.7.17
module load bioinfo/samtools/1.19



for i in *.fasta
do
# commands
 bwa index ${i}

 samtools faidx ${i}


 java -jar -Xmx4g /usr/local/bioinfo/src/picard-tools/picard-2.20.7/picard.jar CreateSequenceDictionary REFERENCE=${i}

done

