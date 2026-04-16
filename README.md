# How Antibiotic Resistance Gene Detection Pipelines Shape Our View of the Resistome

This repository contains files and scripts to generate analysis and figures in the manuscript _How Antibiotic Resistance Gene Detection Pipelines Shape Our View of the Resistome_:

> Juan S. Inda-Díaz, Faith Adegoke, Ulrike Löber, Víctor Hugo Jarquín-Díaz, Yiqian Duan, Johan Bengtsson-Palme, Svetlana Ugarcina Perovic, Luis Pedro Coelho

See also the **interactive** app at the [ARG pipelines](https://arg-pipelines/big-data-biology.org).

## Structure

The folders `dna`,`protein`,`resfinder_dna` contain Snakemake scripts used to run the detection pipelines on the GMGC data. The file `The pipeline.md` goes through each step of the detection.

The folder `code_R_analysis` contains the script `retrieve_aros_abundances_diversity.R` that compiles the results from the detection pipelines, and `plots_2.R` that plots the results pre-computed files and scripts to run the analysis and generate figures included in the manuscript. We direct the users to the [Global Microbial Gene Catalog v1.0](https://gmgc.embl.de/) for access to the source unigene sequences and abundances. 

## ARG detection

### Tools

The tools used to detect ARGs are listed below.

| **Tool** | **Availability** | **Database** |
| :---: | :---: | :---: |
| fARGene (v0.1) | https://github.com/fannyhb/fargene | |
| DeepARG (v2) | https://github.com/gaarangoa/deeparg | |
| AMRFinderPlus (v4.0.15) | https://github.com/ncbi/amr | 2024-12-18.1 |
| RGI (v6.0.3) | https://github.com/arpcard/rgi | CARD (v4.0.0) |
| ResFinder (v2.4.0) | https://github.com/cadms/resfinder | |
| ABRicate v1.0.1 | https://github.com/tseemann/ABRICATE | ARG-ANNOT, CARD, MEGARes v2.0, ResFinder, NCBI (all databases updated 2025-01-14) |

### Normalization

The outputs of DeepARG, AMRFinderPlus, ABRicate, and ResFinder were processed with [argNorm v1.0.0](https://github.com/BigDataBiology/argNorm).


### Preprocessed data

ARG pipeline and normalization outputs as well as the estimated abundance, and richness have been deposited at Zenodo (link) and are available in this repository under `code_R_analysis/data_to_Zenodo`.
