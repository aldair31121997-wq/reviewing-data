# Pipeline for processing SSCS libraries

## Dependencies

The SSCS analysis pipeline is written in Bash and was developed for Linux systems, specifically Ubuntu.

The following software is required:

* [bwa](https://github.com/lh3/bwa), for Illumina read alignment (version 0.7.17)
* [samtools](https://github.com/samtools/samtools), for manipulation of SAM/BAM files (version 1.19)
* [cutadapt](https://github.com/marcelm/cutadapt), for adapter trimming (version 5.1.0)
* [Python](https://www.python.org/downloads) 3.7

The `DSPipeline2` tools and `Consensusmaker.py` originate from the Duplex-Seq-Pipeline developed by the Kennedy Lab:

https://github.com/Kennedy-Lab-UW/Duplex-Seq-Pipeline

In this repository, we re-uploaded the full package because key components of the original pipeline were modified to generate the insertion lists used in this study. We retained the original license notice from the Kennedy Lab pipeline.

## SSCS library analysis

Assuming that the sequencing libraries were prepared as described in the Materials and Methods section of the manuscript, the starting point of this analysis should be paired-end FASTQ files generated from the experimental infections.

To run the analysis, you will need:

* a reference genome sequence in FASTA format corresponding to the parental virus used in the clonal experiment;
* the `Duplex-sequencing/` folder;
* all software listed under **Dependencies** installed on your computational platform, for example on a computing cluster or another high-performance machine.

## Repository structure

After downloading this repository, the `sequencing-analysis/` directory should contain the following folders:

```
Duplex-sequencing/
references/
rawdata/
mutposreportsA/

```

2) you will put your fastq paired end files in your rawdata/ directory, it is important to DO NOT clean the sequences in this step, the pipeline has a cleaning step integrated and you need your full length reads for SSCS construction, to get the paper data
  get into the rawdata/ folder and run the following command:

*IMPORTANT* paper data has not been released yet

```
# BioProject accession
PRJ=PRJNA1481509

# Install tools if needed
sudo apt update
sudo apt install -y sra-toolkit ncbi-entrez-direct

# Get all SRA runs associated with the BioProject
esearch -db sra -query "${PRJ}[BioProject]" | \
  efetch -format runinfo > SraRunInfo_${PRJ}.csv

# Extract SRR accessions
cut -d',' -f1 SraRunInfo_${PRJ}.csv | grep '^SRR' > SraAccList_${PRJ}.txt

# Download SRA files and convert to paired-end FASTQ
mkdir -p sra fastq tmp

prefetch \
  --option-file SraAccList_${PRJ}.txt \
  --max-size u \
  -O sra

while read SRR; do
  fasterq-dump "$SRR" \
    --split-files \
    --threads 8 \
    --outdir fastq \
    --temp tmp
done < SraAccList_${PRJ}.txt

# Compress FASTQ files
gzip -9 fastq/*.fastq

```



4) your reference sequences in fasta format go in references, the example folder provides all neded fasta files to replicate the paper data analysis, get into the references dir and run in place the file called bwa_indexing_ref.sh, this command will index all the references for the subsequent analysis.

```

cd references

sbatch bwa_indexing_ref.sh

```
 

8) now modify the file automaticrun.tsv, this is a table with 4 columns separated by tab, the first one corresponds to the reference of
your sample, the second one is the sample name you want to attribute, I suggest you to format in the following manner:VIRUS-CELLTYPE-REPLICATE, the third collumn corresponds to the read1.fastq file name and location and finally the fourth column is the read2.fastq file name and location.




alignRef                              | Sample           | read1in           | read2in
--------------------------------- | -------------------- | ------------------- | ------------------
 references/HA.full.ShimH5_24a2b.fasta | Shim23a2b-MDCK-1 | rawdata/Shim24a2b-MDCKR1_R1.fastq.gz | rawdata/Shim24a2b-MDCKR1_R2.fastq.gz





6) finally run the script automaticrun.sh, this script takes the information given in the automaticrun.tab table and launches the direct pipeline script, for each sample it will create all intermediate files in `Sample/`
and it will create a table with all the data regarding insertions that will be stocked in Sample/*.mutpos

```

sbatch automaticrun.sh


```

8) next step once all your samples are processed, is that you run the script treatment.sh, it will merge all the tables *mutpos and create a new column with the sample name that you provided, the result will be outputed in the compilatedclean.tab file, while a copy of all the tables *mutpos will be present at mutposreportsA/ 
```

sbatch treatment.sh


```


8) take your table compilatedclean.tab and run insertioncounter.py;

```

python3 insertioncounter.py

```


9) this last step will give you a insertions.xlsx file that contains the following columns

Sample  Cell  Rep

ready for mapping and plotting, the mapping step refers to asign detected insertions to the real positions they emerge acording to our backtrack-model, currently this step is manual and you require the 
backtrack-prediction output to have the insertions and the real positions they correspond, I give you an example of how mapping works:

imagine the following sequence:

CAAAGAAAAAAAAGAGG 

if it aquires an insertion of "GA" product of backtrack that emerges in position 1016, this will create the following sequence:

CAAAGAGAAAAAAAAGAGG

however this sequence can be seen as a sequence with an insertion of AG in front of position 1014 or a sequence with an insertion of GA exactly in position 1016; for this reason you need the plot provided 
by backtrack-prediction, this plot will tell you the position of insertions, you must verify manually if the sequence they create corresponds to the sequence reported by the insertions obseved, if thats the
case you can re-map the position of the insertion. This step is intensive and probably will be automated in a future version of this repository, but for the moment you must do it manually.









