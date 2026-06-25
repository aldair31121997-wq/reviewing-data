# Pipeline for processing of SSCS libraries.


## Dependencies
The SSCS analysis pipeline is written in bash, on linux OS (ubuntu):
* [bwa](https://github.com/lh3/bwa) illumina sequence alignment (version 0.7.17)
* [samtools](https://github.com/samtools/samtools) manipulation of sam files (version 1.19)
* [cutadapt](https://github.com/marcelm/cutadapt) removal of adapter sequences (version 5.15.0)
* [python](https://www.python.org/downloads) (version 2.7.18)


the DSPipeline2 tools and Consensusmaker.py come from https://github.com/Kennedy-Lab-UW/Duplex-Seq-Pipeline here we decided to re-upload the whole package as we modified key components of the original pipeline to get the insertion list, we keept the original licence advise coming from kennedylab.



# SSCS library Analysis.
Asuming you implemented your library construction as described in material and methods, you should have now a paired end R1 and R2 couple of fastq files that come from your experimental infections, for this analysis you require 3 things, a reference sequence (corresponding to the parent virus of your clonal experiment), the DuplexSequencing folder and the #dependencies1 installed in your work machine (a cluster or any other high resource machine)

1) first step is to download this repository, inside you should have 4 total dir, Duplex-sequencing/, Unifiedworkflow/, references/, rawdata/.


2) you will put your fastq paired end files in your rawdata/ directory, it is important to DO NOT clean the sequences in this step, the pipeline has a cleaning step integrated and you need your full length reads for SSCS construction, to get the example data
  get into the rawdata/ folder and run the following command:





4) your reference sequences in fasta format go in references, the example folder provides all neded fasta files and indexed files required for the analysis of the example and paper data, however, IF YOU ADD NEW DATA WITH DIFFERENT VIRUSES THAN THE EXAMPLE AND THE PAPER you must add the fasta reference in this folder an run in place the file inside called index.sh

5) now you modify the file automaticrun.tsv, this is a table with 4 columns separated by tab, the first one corresponds to the reference of
your sample, the second one is the sample name you want to attribute, I suggest you to format in the following manner:VIRUS-CELLTYPE-REPLICATE, the third collumn corresponds to the read1.fastq file name and location and finally the fourth column is the read2.fastq file name and location.




alignRef                              | Sample           | read1in           | read2in
--------------------------------- | -------------------- | ------------------- | ------------------
 references/HA.full.ShimH5_24a2b.fasta | Shim23a2b-MDCK-1 | rawdata/Shim24a2b-MDCKR1_R1.fastq.gz | rawdata/Shim24a2b-MDCKR1_R2.fastq.gz





6) finally run the script automaticrun.sh, this script takes the information given in the automaticrun.tab table and launches the direct pipeline script, for each sample it will create all intermediate files in Sample/
and finally it will create a table with all the data regarding insertions that will be stocked in Sample/*.mutpos

7) next step once all your samples are processed, is that you run the script treatment.sh, it will merge all the tables *mutpos and create a new column with the sample name that you provided

8) take your table with all data, completedata.tab and run insertioncounter.py, this last step will give you a table that contains the following columns

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


