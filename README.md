# HA insertion

This GitHub repository contains the bioinformatics resources required to reproduce the analyses presented in the preprint *“Insertion-driven emergence of highly pathogenic avian influenza viruses is governed by product–template duplex thermodynamics.”*

In our study, polymerase slippage is proposed as the main mechanism driving recurrent insertions at the HA cleavage site of H5 and H7 influenza viruses. Slippage is determined by sequence-dependent RNA–RNA product–template interactions within the catalytic site of the viral RNA-dependent RNA polymerase (RdRp). This model is evaluated using `HPAIVpredict`, our slippage-prediction tool, and validated experimentally through the analysis of sequencing data generated using PA-HA and minigenome systems.

## Reproducibility scope

This repository contains three complementary components.

### `HPAIVpredict/`

`HPAIVpredict/` contains our slippage-based prediction tool. It is designed to predict the occurrence and identity of insertions in an RNA sequence based on sequence properties that modulate product–template duplex interactions within the catalytic site of the viral RdRp.

This directory constitutes the main reproducible computational component of the repository. The workflow can be executed using the provided example dataset to regenerate the prediction tables and plots presented in the study.

### `sequencing-analysis/`

`sequencing-analysis/` documents the SSCS sequencing-analysis pipeline used to analyse data generated from PA-HA and minigenome products.

This workflow depends on local computing infrastructure, reference-sequence indexing, raw FASTQ availability, and sample-specific configuration files. It is therefore provided for transparency and reuse, but is not intended to function as a fully machine-independent, one-click workflow.

### `structure-analysis/`

`structure-analysis/` provides the sequences and instructions used to perform the RNA structure analyses included in the study.

Some steps rely on external software or web-server execution and are therefore documented rather than fully encapsulated within the repository.
