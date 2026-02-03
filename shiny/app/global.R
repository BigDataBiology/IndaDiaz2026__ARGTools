library(shiny)
library(bslib)
library(tidyverse)
library(shinyWidgets)
#library(plotly)
library(dplyr)
library(ggplot2)
#library(ggVennDiagram)
library(gridExtra)
library(RColorBrewer)
library(ggpattern)
library(grid)
#library(eulerr)
library(reactable)
library(Cairo)
library(ggalluvial)
library(cowplot)
library(scales)

options(dplyr.summarise.inform = FALSE)
options(shiny.usecairo = TRUE)
options(shiny.userrender.type = "svg")

# Path Configuration
DATA_DIR <- "../../code_R_analysis/output_abundance_diversity_resistome"
METADATA_PATH <- "../../data/metadata_GMGC10.sample.meta.tsv"

# Constants
general_size <- 10
pal_7 <- brewer.pal(8, "BrBG")[c(4, 5, 3, 6, 2, 7, 1)]
pal_10_q <- pal_7[c(1, 2, 3, 4, 5, 5, 6, 6, 7, 7)]
pal_10_complete <- brewer.pal(10, "BrBG")
pal_10_complete <- pal_10_complete[c(-1,-10)]

# pattern for plots

pattern_density <- 0.001 
pattern_spacing <- 0.025
pattern_fill <- "white"
pattern_size <- 0.12

# all HABITATS
EN <- c("human gut", "human oral",  "human skin", "human nose", "human vagina", 
        "dog gut", "cat gut", "mouse gut", "pig gut", "wastewater", "marine", 
        "freshwater", "soil" , "amplicon", "isolate",  "built-environment" )

# SOURCE FOR EACH HABITAT

SO <- c(rep("humans", 5), rep("mammals", 4),  
        "wastewater", "marine", "freshwater", 
        "soil", rep("other", 2), "built-environment")

names(SO) <- EN

not_env <- c("amplicon", "isolate", "built-environment" )

EN2  <- EN[!EN %in% not_env]
h2 <- c("humans","mammals","wastewater", "freshwater", "soil", "marine")

# Tool Definitions
tools_levels <- c("DeepARG", "fARGene", "ABRicate-ARGANNOT", "ABRicate-MEGARes", 
                  "RGI-DIAMOND", "ABRicate-CARD", "AMRFinderPlus", 
                  "ABRicate-NCBI", "ResFinder", "ABRicate-ResFinder")

tools_labels <- c("DeepARG", "fARGene", "ABRicate-ARGANNOT", "ABRicate-MEGARes", 
                  "RGI-CARD", "ABRicate-CARD", "AMRFinderPlus-NCBI", 
                  "ABRicate-NCBI", "ResFinder", "ABRicate-ResFinder")

tools_texture <- c("ABRicate-ARGANNOT", "ABRicate-MEGARes", "ABRicate-CARD", 
                   "ABRicate-NCBI", "ABRicate-ResFinder")


tool_choices <- c("DeepARG" = "DeepARG", "fARGene" = "fARGene", 
                  "ABRicate-ARGANNOT" = "ABRicate-ARGANNOT" , 
                  "ABRicate-MEGARes" = "ABRicate-MEGARes", 
                  "RGI-CARD" = "RGI-DIAMOND", "ABRicate-CARD" = "ABRicate-CARD", 
                  "AMRFinderPlus-NCBI" = "AMRFinderPlus", 
                  "ABRicate-NCBI" = "ABRicate-NCBI", "ResFinder" = "ResFinder", 
                  "ABRicate-ResFinder" = "ABRicate-ResFinder")


# Source External Scripts
source("../../code_R_analysis/helper.R")
source("load_data.R")

# Load Data
data_list <- list()

data_list <- c(data_list, 
               load_results_tools(
                 DATA_DIR = DATA_DIR, 
                 file = "results_tools.rds"))

data_list <- c(data_list, 
               load_metadata(
                 DATA_DIR = DATA_DIR,  
                 file = "metadata_GMGC10.sample.meta.tsv"))


data_list <- c(data_list, 
               load_abundances(
                 DATA_DIR = DATA_DIR, 
                 file1 = "abundance_diversity_part1.rds", 
                 file2 = "abundance_diversity_part2.rds",
                 file3 = "abundance_diversity_part3.rds",
                 data_list$metadata))
  

data_list <- c(data_list, 
               load_pan_core(
                 DATA_DIR = DATA_DIR,
                 core_file = "core_resistome.rds",
                 pan_file = "pan_resistome.rds"))

top20 <- c("van", "efflux pump", "cell wall charge", "rpoB", "tet RPG", "class A beta-lactamase", 
           "aph", "MFS efflux pump", "erm", "target-modifying enzyme", "tet enzyme")

abundance_threshold <- load_abundances_thresholds(
  DATA_DIR = "../../code_R_analysis/output_abundance_diversity_resistome", 
  file = "abundance_diversity_60.rds",
  data_list$metadata)

data_list <- c(data_list, 
               list(abundance60 = abundance_threshold$abundance,
               abundance_tool_sample60 = abundance_threshold$abundance_tool_sample,
               abundance_class60 = abundance_threshold$abundance_class))

abundance_threshold <- load_abundances_thresholds(
  DATA_DIR = "../../code_R_analysis/output_abundance_diversity_resistome", 
  file = "abundance_diversity_70.rds",
  data_list$metadata)

data_list <- c(data_list, 
               list(abundance70 = abundance_threshold$abundance,
               abundance_tool_sample70 = abundance_threshold$abundance_tool_sample,
               abundance_class70 = abundance_threshold$abundance_class))

abundance_threshold <- load_abundances_thresholds(
  DATA_DIR = "../../code_R_analysis/output_abundance_diversity_resistome", 
  file = "abundance_diversity_80.rds",
  data_list$metadata)

data_list <- c(data_list, 
               list(abundance80 = abundance_threshold$abundance,
               abundance_tool_sample80 = abundance_threshold$abundance_tool_sample,
               abundance_class80 = abundance_threshold$abundance_class))
