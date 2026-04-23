

# Antibiotic resistance gene detection on the Global Microbial Gene Catalog

## Overview

This dataset contains antibiotic resistance gene (ARG) detection results from unigenes in the Global Microbial Gene Catalog v1.0 ([GMGC](https://gmgc.embl.de/)) using multiple widely used detection pipelines (including DeepARG, RGI, fARGene, ResFinder, AMRFinderPlus, and ABRicate with multiple databases).

Outputs were harmonized using argNorm and manually curated, estimations of ARG abundance and richness as well as pan- and core-resistome were done across 11,519 metagenomic samples spanning 13 habitats (human-associated, animal-associated, and environmental), at the ARG class level. 

This resource supports comparative analyses of how ARG detection pipelines influence resistome interpretation.

---

## Tools and databases 

ARG detection was performed using:

| **Tool** | **Availability** | 
| :---: | :---: | 
| fARGene (v0.1) | https://github.com/fannyhb/fargene | 
| DeepARG (v2) | https://github.com/gaarangoa/deeparg | 
| AMRFinderPlus (v4.0.15), database 2024-12-18.1 | https://github.com/ncbi/amr |
| RGI (v6.0.3), database CARD (v4.0.0) | https://github.com/arpcard/rgi | 
| ResFinder (v2.4.0) | https://github.com/cadms/resfinder | 
| ABRicate v1.0.1, databases: ARG-ANNOT, CARD, MEGARes v2.0, ResFinder, and NCBI (all updated 2025-01-14) | https://github.com/tseemann/ABRICATE |  


## Normalization

Outputs from DeepARG, AMRFinderPlus, ABRicate, and ResFinder were standardized using [argNorm v1.0.0](https://github.com/BigDataBiology/argNorm).

---

## Datasets contents

### 1. Pipeline outputs 

`pipelines_output/pipeline`

Contains ARG predictions per pipeline after normalization and manual curation.

Each record includes:

- ARO: identifiers for the reported unigene (after normalization)
- ARO\_name: (ARO identifier annotation)
- pipeline: pipeline identifier
- id: standardize name for the identity level column

As well as resistance ontology mapping:

- parent: first level of ARO aggregation
- parent\_description: description of the first level of aggregation
- geneclass\_argcompare: gene class - final level of aggregation, should perhaps rename it, in the abundance file it’s called gene class


**Pipeline identifiers.** Pipeline names encode tool, database, and sequence type (nucleotide or amino acid). Additional filtered variants (70/80/90% identity thresholds) were used for analyses.


### 2. Unified ARG unigene list

`reported_unigenes_as_ARG_per_habitat.csv`

List of all unigenes detected as ARGs by any pipeline, stratified by habitat.


### 3. Ontology mapping

`conversion_aro_geneclass.csv`

Maps ARO terms to hierarchical ontology levels and curated ARG classes.

The columns include: 

- Term\_ID: ARO identifier
- Term\_Label: description of the ARO identifier
- Parent\_ID: ARO of the first higher aggregation
- Parent\_Label: description of the Parent\_ID
- geneclass: manually curated gene class


### 4. Abundance and richness

`abundance_richness.csv.gz`

ARG abundance and richness across 11,519 metagenomic samples from [GMGC](https://gmgc.embl.de/). The metagenomic samples span 13 habitats, including human-associated (gut, skin, nose, oral, and vaginal), non-human and host-associated (pig, mouse, dog, and cat), and external habitats (wastewater, marine, freshwater, and soil). The abundance and richness is aggregated at the ARG class level.

Columns:

- sample: [GMGC](https://gmgc.embl.de/) metagenome ID
- geneclass: curated ARG class
- pipeline: detection pipeline
- abundance: normalized read counts (per million, gene-length corrected) and aggregated at the ARG class level
- richness: rarefied ARG diversity, number of unique unigenes detected as ARG each metagenomic sample
- richness\_no\_rarified: aw ARG diversity


### 5. Unigene-level abundance

`args_abundances.tsv.gz`

Abundance of individual ARG-associated unigenes ([GMGC](https://gmgc.embl.de/) format), scaled per 10 million reads per metagenome.


### 6. Metadata

`metagenomes_metadata.csv`

Columns:
- sample\_id: metagenome identifier
- habitat
- insertsHQ: high-quality reads


