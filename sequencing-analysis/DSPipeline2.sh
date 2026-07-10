#!/bin/bash
### 20221130 Gabriel Dupré - Aldair Martinez
### script for NGS SSCSreads analysis
#SBATCH -p workq
#SBATCH -J SSCSxxx
#SBATCH -e error_SSCSxxxx
#SBATCH --mem=40G

#load modules, replace the bioinfo/ section for the location of your packages in your work machine
module purge
module load bioinfo/bwa/0.7.17
module load bioinfo/samtools/1.19
module load bioinfo/Cutadapt/4.3
module load devel/python/Python-2.7.18

# Setup variables for run:
clear

# Set up error checking
# Stop on any error
set -e
# Stop on an error inside a pipeline
set -o pipefail
# Throw an error on calling an unassigned variable
set -u

#SET_PARAMETERS, replace the reference name, sample, read1in and read2in by the samples you want to work with.
alignRef="$1"
Sample="$2"
read1in="$3"
read2in="$4"
minMem=3
maxMem=1000
cutOff=0.7
NcutOff=0.3
barcodeLength=8
spacerLength=0
repFilt=9

#LOG_FILE_NAME
logFile=${Sample}.log.txt

#Export all variables
export alignRef
export Sample
export read1in
export read2in
export minMem
export maxMem
export cutOff
export NcutOff
export barcodeLength
export spacerLength
export repFilt

# Print out options used to log file
touch $logFile
echo "Sample: " $Sample | tee -a ${logFile}
echo "Reference genome: " $alignRef | tee -a ${logFile}
echo "Barcode length: " $barcodeLength | tee -a ${logFile}
echo "Spacer length: " $spacerLength | tee -a ${logFile}
echo "Minimum family size: " $minMem | tee -a ${logFile}
echo "Maximum family size: " $maxMem | tee -a ${logFile}
echo "Consensus cutoff: " $cutOff | tee -a ${logFile}
echo "Consensus N cutoff: " $NcutOff | tee -a ${logFile}
echo "Rep filt" $repFilt | tee -a ${logFile}


## 1.Create unaligned bam file from the two fqreads, replace /usr/local/bioinfo/src/picard-tools/picard-2.20.7 by the location of piccard tools in your work machine
java -jar -Xmx4g /usr/local/bioinfo/src/picard-tools/picard-2.20.7/picard.jar FastqToSam FASTQ=$read1in FASTQ2=$read2in O=${Sample}.unaligned.pe.bam SM=$Sample

## 2. Run UnifiedConsensusMaker
echo "Starting Run" | tee -a ${logFile}
echo "UnifiedConsensusMaker runs"  | tee -a ${logFile}
date | tee -a ${logFile}
echo "" | tee -a ${logFile}
python ../Duplex-Sequencing/UnifiedConsensusMaker.py --input ${Sample}.unaligned.pe.bam --prefix $Sample --taglen $barcodeLength --spacerlen $spacerLength --minmem $minMem --maxmem $maxMem --tagstats --cutoff $cutOff --Ncutoff $NcutOff --write-sscs --rep_filt $repFilt

### 3. avant enlever les infos supp au niveau de la 3ieme ligne "+chiffre" qui bloque cutadapt, renommer fichier *sscs.f.fq
echo "Modifying SSCS fq file for cutadapt compatibility" | tee -a ${logFile}
gunzip *fq.gz
cat ${Sample}_read1_sscs.fq | sed -e "s/+[^a-zA-Z]*/+/g" > ${Sample}_read1_sscs.f.fq

##run cutadapt

cutadapt -a AGATCGGAAGAGC -o ${Sample}_read1_sscs.cuta.fq  -m 36 --trim-n -q 20,20 --max-n 0 ${Sample}_read1_sscs.f.fq > cutadapt.log 

### 4. alignement bwa mem des reads dcs sur séquence de référence, puis créaction d'un fichier .bam, avec bwamem_trimmedsscs1_sort_filt_UniConsMakDCS.sh
echo "Run bwa mem" | tee -a ${logFile}
#alignment
bwa mem $alignRef ${Sample}_read1_sscs.cuta.fq > ${Sample}_trimmedsscs1.sam
#convert to bam and sort
echo "Convert sam to bam and sorting" | tee -a ${logFile}
samtools view -Sbu ${Sample}_trimmedsscs1.sam | samtools sort -o ${Sample}_trimmedsscs1.sort.bam
#index sorted bam file
echo "Indexing sorted bam file" | tee -a ${logFile}
samtools index ${Sample}_trimmedsscs1.sort.bam
#filter out unmapped reads from the final dcs.pe.sort.bam file
echo "Filter out unmapped reads" | tee -a ${logFile}
samtools view -F 4 -b ${Sample}_trimmedsscs1.sort.bam > ${Sample}.trimmedsscs1.sort.filt.bam

### 5. add read group field to the header of the final DCS.bam file to allow compatibility with GATK using picard tools + indexing, replace /usr/local/bioinfo/src/picard-tools/picard-2.20.7 by the location of piccard tools in your work machine
echo "Add read group field to the header of the final bam file" | tee -a ${logFile}
java -jar -Xmx4g /usr/local/bioinfo/src/picard-tools/picard-2.20.7/picard.jar AddOrReplaceReadGroups INPUT=${Sample}.trimmedsscs1.sort.filt.bam OUTPUT=${Sample}.cleanedsscs1.sort.filt.readgroups.bam RGLB=UW RGPL=Illumina RGPU=ATATAT RGSM=default
echo "Indexing" | tee -a ${logFile}
samtools index ${Sample}.cleanedsscs1.sort.filt.readgroups.bam

### 6. local realignment: 2 steps
# before create .fai index file of the reference with samtools faidx function
# before create .dict index file of the reference with picard 2.14.1 CreateSequenceDictionary function, replace /usr/local/bioinfo/src/GATK/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/ by the location of GATK in your work machine
echo "local realignment" | tee -a ${logFile}
java -jar /usr/local/bioinfo/src/GATK/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar -T RealignerTargetCreator -R $alignRef -I ${Sample}.cleanedsscs1.sort.filt.readgroups.bam -o ${Sample}.cleanedsscs1.sort.filt.readgroups.intervals
java -jar /usr/local/bioinfo/src/GATK/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar -T IndelRealigner -R $alignRef -I ${Sample}.cleanedsscs1.sort.filt.readgroups.bam -targetIntervals ${Sample}.cleanedsscs1.sort.filt.readgroups.intervals -o ${Sample}.cleanedsscs1.sort.filt.readgroups.realign.bam

### 7. end-trimming, remove 5 bases at both 5' and 3' ends
## no endtrimming was done for the moment beacause bwa mem already performed soft clipping in its algorithm and we have already 3'end trimmed our sscs1 reads before
## java -jar /usr/local/bioinfo/src/GATK/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar -T ClipReads -I UniConsMak.dcs.pe.sort.filt.readgroups.realign.bam -o DCS.UniConsMak_final.bam -R /home/gdupre/work/HTS/50-778968403/ref/HA.full.H5N8_mut1_ref.fasta --cyclesToTrim "1-5,138-142" --clipRepresentation SOFTCLIP_BASES

### 8. make pileup file from the final DCS reads:
echo "mpileup construction" | tee -a ${logFile}
samtools mpileup -B -A -d 500000 -f $alignRef ${Sample}.cleanedsscs1.sort.filt.readgroups.realign.bam > ${Sample}.sscs1_final.bam.pileup

### 9. Count the number of unique mutations present in the final DCS sequences and calculate their frequencies
echo "Count the number of unique mutations" | tee -a ${logFile}
cat ${Sample}.sscs1_final.bam.pileup | python ../Duplex-Sequencing/Nat_Protocols_Version/CountMuts.py -d 100 -c 0 -C 0.99 -u > ${Sample}.sscs1_final.bam.pileup.countmuts

### 10. locate the genomic position of each mutation
echo "genomic location of each mutation" | tee -a ${logFile}
python ../Duplex-Sequencing/Nat_Protocols_Version/mut_pos_v2.py -i ${Sample}.sscs1_final.bam.pileup -o ${Sample}.sscs1_final.bam.pileup.mutpos -d 100 -C 0.99

# Step 11: Finishing
echo "Finishing with run.. " $Sample | tee -a ${logFile}
echo " Run finished.." | tee -a ${logFile}
