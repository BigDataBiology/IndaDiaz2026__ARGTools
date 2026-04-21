# Running the R Scripts

---

## Overview

This repository contains the R scripts used to generate figures and analysis results. It covers:

---

### R Analysis File Structure

```
IndaDiaz2026__ARGTools/
├── pixi.toml                                  # Pixi environment 
├── code_R_analysis/
│   ├── helper.R                               # Shared functions used across scripts
│   ├── plots_2.R                              # Main figure generation script
│   ├── retrieve_aros_abundances_diversity.R   # ARG abundance & diversity retrieval
│   ├── generate_data_to_zenodo.R              # Prepares data for Zenodo deposit
│   ├── data_for_shinny.R                      # Prepares data for Shiny app
│   ├── output_abundance_diversity_resistome/  # Intermediate .rds data files
│   └── output_plots/                          # Generated figures (SVG)
└── data/
    └── metadata_GMGC10.sample.meta.tsv        # Sample metadata
```

---

## Installation

This pipeline uses [Pixi](https://pixi.sh) to manage all R dependencies reproducibly.

### 1. Install Pixi

```bash
curl -fsSL https://pixi.sh/install.sh | bash
```

### 2. Clone the repository

```bash
git clone https://github.com/indajuan/IndaDiaz2026__ARGTools.git
cd IndaDiaz2026__ARGTools
```

### 3. Install the environment

```bash
pixi install
```

This will automatically install R and all required packages as defined in `pixi.toml`.

---

## Usage

All commands should be run from the **project root** (`IndaDiaz2026__ARGTools/`), not from inside `code_R_analysis/`.


### Running the Scripts

```bash
# 1. Generate all figures
pixi run Rscript code_R_analysis/plots_2.R

# 2. Generate data files for Zenodo
pixi run Rscript code_R_analysis/generate_data_to_zenodo.R

# 3. Prepare Shiny app data
pixi run Rscript code_R_analysis/data_for_shinny.R
```

Output figures are saved as SVG files in `code_R_analysis/output_plots/`.

