# Pipeline for processing SSCS libraries

## Dependencies

The SSCS analysis pipeline is written in Bash and was developed for Linux systems, specifically Ubuntu.

The following software is required:

* [bwa](https://github.com/lh3/bwa), for Illumina read alignment (version 0.7.17)
* [samtools](https://github.com/samtools/samtools), for manipulation of SAM/BAM files (version 1.19)
* [cutadapt](https://github.com/marcelm/cutadapt), for adapter trimming (version 5.1.0)
* [Python](https://www.python.org/downloads) version 2.7 is used in `automaticrun.sh`, while Python >=3.7 can be used for the rest of the pipeline
* [Pysam](https://pysam.readthedocs.io/en/latest/api.html), directly available in Python 2.7 and required as a module in Python >=3.7

The `DSPipeline2` tools and `Consensusmaker.py` originate from the Duplex-Seq-Pipeline developed by the Kennedy Lab:

https://github.com/Kennedy-Lab-UW/Duplex-Seq-Pipeline

In this repository, we re-uploaded the full package because key components of the original pipeline were modified to generate the insertion lists used in this study. We retained the original license notice from the Kennedy Lab pipeline. All Python packages derived from this pipeline are written in Python 2.7.

## SSCS library analysis

Assuming that the sequencing libraries were prepared as described in the Materials and Methods section of the manuscript, the starting point of this analysis should be paired-end FASTQ files generated from the experimental infections.

To run the analysis, you will need:

* a reference genome sequence in FASTA format corresponding to the parental virus used in the clonal experiment;
* the `Duplex-sequencing/` folder;
* all software listed under **Dependencies** installed on your computational platform, for example, on a computing cluster or another high-performance machine.

## Repository structure

After downloading this repository, the `sequencing-analysis/` directory should contain the following folders:

```text
Duplex-sequencing/
references/
rawdata/
mutposreportsA/
```

## Analysis workflow

### 1. Add the paired-end FASTQ files

Place your paired-end FASTQ files in the `rawdata/` directory.

It is important not to clean the sequences at this step. The pipeline has an integrated cleaning step, and untrimmed raw reads are required for SSCS construction.

To obtain the paper data, we provide the accession number of the currently unreleased SRA BioProject. A secondary FileSender link is also provided, where a dataset derived from a PA-HA-Shim24a2b infection in MDCK cells is available as example data.

**IMPORTANT: The paper data have not yet been released.**

```bash
# BioProject accession — NOT RELEASED YET
PRJ=PRJNA1481509

# Example FASTQ files available via FileSender
https://filesender.renater.fr/?s=download&token=ac6b07ab-9ad8-448c-a888-9089670daf42
```

### 2. Index the reference sequences

Place your reference sequences in FASTA format in the `references/` directory.

The example folder provides all the FASTA files needed to replicate the paper data analysis. Navigate to the `references/` directory and run the file called `bwa_indexing_ref.sh`. This command will index all the references for the subsequent analysis.

```bash
cd references

sbatch bwa_indexing_ref.sh
```

### 3. Configure `automaticrun.tsv`

Modify the `automaticrun.tsv` file.

This is a table containing four tab-separated columns. The first column corresponds to the reference sequence for your sample. The second column corresponds to the sample name that you want to assign. We suggest formatting it in the following manner: `VIRUS-CELLTYPE-REPLICATE`. The third column corresponds to the name and location of the Read 1 FASTQ file. Finally, the fourth column corresponds to the name and location of the Read 2 FASTQ file.

| alignRef                              | Sample           | read1in                              | read2in                              |
| ------------------------------------- | ---------------- | ------------------------------------ | ------------------------------------ |
| references/HA.full.ShimH5_24a2b.fasta | Shim23a2b-MDCK-1 | rawdata/Shim24a2b-MDCKR1_R1.fastq.gz | rawdata/Shim24a2b-MDCKR1_R2.fastq.gz |

### 4. Run `automaticrun.sh`

Run the `automaticrun.sh` script.

This script takes the information provided in the `automaticrun.tsv` table and launches the pipeline directly. For each sample, it will create all intermediate files in `Sample/`, for example:

```text
Shim24a2b-MDCK-1/
```

It will also create a table containing all the insertion data, which will be stored as a `.mutpos` file in the sample directory, for example:

```text
Shim24a2b-MDCK-1/*.mutpos
```

Before launching the analysis, adapt the following line to the number of samples being analysed, corresponding to the number of lines in your table. For example, for four samples:

```bash
#SBATCH --array=1-4
```

Launch the analysis with:

```bash
sbatch automaticrun.sh
```

### 5. Run `treatment.sh`

Once all your samples have been processed, run the `treatment.sh` script.

This script will merge all the `.mutpos` tables and create a new column containing the sample name that you provided. The result will be written to the `compilatedclean.tab` file, while a copy of all the `.mutpos` tables will be placed in the `mutposreportsA/` directory.

```bash
sbatch treatment.sh
```

### 6. Run `insertioncounter.py`

Use the `compilatedclean.tab` table as input for `insertioncounter.py`:

```bash
python3 insertioncounter.py compilatedclean.tab
```

### 7. Interpret and map the output

This final step will generate an `insertions.xlsx` file containing the following columns:

| sample | Pos | maxrate | insertion | FreqInser | Template | Depths | rep | cell |
| ------ | --- | ------- | --------- | --------- | -------- | ------ | --- | ---- |

The `sample` column contains the sample name. `Pos` contains a given position in the HA sequence that is being analysed. `maxrate` contains the total number of insertions detected at a given position. `insertion` contains the size and sequence of the observed insertion.

`FreqInser` contains the total frequency of the observed insertion, obtained by dividing the total number of insertions detected at a given position, `maxrate`, by the total sequencing depth at that position, `Depths`. The `rep` column contains the replicate number of the experiment, and the `cell` column contains the cell type in which the experiment was performed.

This dataset is ready for mapping and plotting.

The mapping step consists of assigning detected insertions to the actual positions from which they emerged according to our backtracking model. Currently, this step is performed manually and requires the backtracking-prediction output to identify the insertion sequences and the actual positions to which they correspond.

The following example illustrates how the mapping works.

Consider the following sequence:

```text
CAAAGAAAAAAAAGAGG
```

If it acquires a `GA` insertion resulting from a backtracking event emerging at position 1016, it will produce the following sequence:

```text
CAAAGAGAAAAAAAAGAGG
```

However, this sequence can be interpreted either as containing an `AG` insertion before position 1014 or as containing a `GA` insertion exactly at position 1016.

For this reason, you need the plot provided by the backtracking prediction. This plot indicates the predicted insertion positions. You must manually verify whether the sequence generated by a predicted insertion corresponds to the observed insertion sequence. If this is the case, you can remap the position of the insertion accordingly.

This step is labor-intensive and will probably be automated in a future version of this repository. For the moment, it must be performed manually.
