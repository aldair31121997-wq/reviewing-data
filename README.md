# HA-insertion

This github repository contains all the bioinformatic resources coming from the paper "Insertion-driven emergence of high pathogenicity avian influenza viruses is governed by product–template duplex thermodynamics".

## Reproducibility scope

This repository contains three complementary components.

`HPAIVpredict/` is the main reproducible computational component and can be executed on the provided example dataset to regenerate the prediction tables and plots.

`sequencing-analysis/` documents the SSCS sequencing-analysis pipeline used in the study. This workflow depends on local computing infrastructure, reference indexing, raw FASTQ availability, and sample-specific configuration files. It is therefore provided for transparency and reuse, but is not intended as a fully machine-independent one-click workflow.

`stationary and transient structure analysis/` provides the sequences and instructions used for RNA structure analyses. Some steps rely on external tools or web-server execution and are therefore documented rather than fully encapsulated.
