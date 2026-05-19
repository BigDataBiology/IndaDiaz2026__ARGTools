# How Antibiotic Resistance Gene Detection Pipelines Shape Our View of the Resistome

This repository contains files and scripts to generate analysis and figures in the manuscript _How Antibiotic Resistance Gene Detection Pipelines Shape Our View of the Resistome_:

> Juan S. Inda-Díaz, Faith Adegoke, Ulrike Löber, Víctor Hugo Jarquín-Díaz, Yiqian Duan, Johan Bengtsson-Palme, Svetlana Ugarcina Perovic, Luis Pedro Coelho

See also the **interactive** app at the [ARG pipelines](https://arg-pipelines.big-data-biology.org/).

## Structure

The folders `dna`, `protein`, and `resfinder_dna` contain Snakemake workflows used to run the ARG detection pipelines on the GMGC data. The Snakemake workflows and `The pipeline.md` document describe how each tool was executed within this study, including command structure, input/output organization, and downstream integration into the analysis workflow.

The folder `code_R_analysis` contains the script `retrieve_aros_abundances_diversity.R`, which compiles and harmonizes the results from the detection pipelines, and `plots_2.R`, which generates the figures and analyses presented in the manuscript from precomputed files.

We direct users to the [Global Microbial Gene Catalog v1.0](https://gmgc.embl.de/) for access to the source unigene sequences, metagenomic sample metadata, and abundance estimates.

Installation instructions, software dependencies, databases, and usage information for each ARG detection tool are maintained by the respective software developers. Therefore, we refer users to the official repositories and documentation pages of each tool listed below for installation and configuration instructions. 

The detection tools were ran on Linux (RHEL 9-compatible distribution; kernel 5.14.x; x86_64 architecture). 

## ARG detection

### Tools

The tools used to detect ARGs are listed below.

| **Tool** | **Availability** | 
| :---: | :---: | 
| fARGene (v0.1) | https://github.com/fannyhb/fargene | 
| DeepARG (v2) | https://github.com/gaarangoa/deeparg | 
| AMRFinderPlus (v4.0.15), database 2024-12-18.1 | https://github.com/ncbi/amr |
| RGI (v6.0.3), database CARD (v4.0.0) | https://github.com/arpcard/rgi | 
| ResFinder (v2.4.0) | https://github.com/cadms/resfinder | 
| ABRicate v1.0.1, databases: ARG-ANNOT, CARD, MEGARes v2.0, ResFinder, and NCBI (all updated 2025-01-14) | https://github.com/tseemann/ABRICATE |  

### Normalization

The outputs of DeepARG, AMRFinderPlus, ABRicate, and ResFinder were processed with [argNorm v1.0.0](https://github.com/BigDataBiology/argNorm).

### Workflow

[Snakemake v8.27.1](https://snakemake.github.io/)

### Gene clustering

[VSEARCH v2.30.0](https://github.com/torognes/vsearch)

### Preprocessed data

ARG pipeline and normalization outputs as well as the estimated abundance, and richness have been deposited at Zenodo (link) and are available in this repository under `code_R_analysis/data_to_Zenodo`.

### Compilation and analysis of results 

Analyses were performed using R version 4.5.2 (2025-10-31).

| Package        | Version   |
|----------------|----------|
| ggbreak        | v0.1.7   |
| scales         | v1.4.0   |
| cowplot        | v1.2.0   |
| Cairo          | v1.7-0   |
| ggpattern      | v1.3.1   |
| RColorBrewer    | v1.1-3   |
| lubridate      | v1.9.5   |
| forcats        | v1.0.1   |
| stringr        | v1.6.0   |
| purrr          | v1.2.2   |
| readr          | v2.1.6   |
| tidyr          | v1.3.2   |
| tibble         | v3.3.1   |
| tidyverse      | v2.0.0   |
| gridExtra      | v2.3     |
| ggplot2        | v4.0.2   |
| dplyr          | v1.2.1   |


### Installation time
Typical installation time of the ARG detection tools on a standard desktop computer is approximately 1–3 hours, depending on internet bandwidth and system configuration. This includes installation of all software dependencies and setup of required databases. No compilation of core software is required for most tools, as precompiled binaries or conda packages are used.

No special hardware needed — standard compute is sufficient.

### Expected runtime for ARG detection tools

Approximate runtimes for the demo dataset (~1000 gene sequences in FASTA format) on a standard desktop computer (8–16 CPU cores, 16–32 GB RAM) are shown below.

| Tool | Version | Expected runtime |
|---|---|---|
| RGI (DIAMOND mode) | ~1–10 minutes |
| DeepARG  | ~1–10 minutes |
| ResFinder  | ~1–10 minutes |
| fARGene (all models) | ~5–20 minutes |
| ABRicate | ~1-10 minutes per database |
| AMRFinderPlus | ~1–10 minutes |
| argNorm | ~1–10 minutes |
| VSEARCH | ~1–10 minutes |
`
Full-scale analyses on the GMGC dataset were performed on HPC infrastructure and are not intended for execution on standard desktop machines.

### Outputs

The ouputs of each tool are different. Examples of input fasta files and ouputs can be found in `test/`.


