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
library(tidyr)

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




# Precomputing everything needed to speed up load time for abundance and diversity overview

abundance_tools_excl <- c("DeepARG", "RGI-DIAMOND")

# --- Ensure factor levels are consistent across all abundance tables ---
# Sample as factor once (master)
data_list$abundance <- data_list$abundance %>%
  dplyr::mutate(sample = factor(sample))

# Make sure tool is a factor too (and keep the full tool universe from the master table)
# If it's already a factor, this is harmless.
data_list$abundance <- data_list$abundance %>%
  dplyr::mutate(tool = factor(tool))

all_samples <- levels(data_list$abundance$sample)
all_tools   <- levels(data_list$abundance$tool)  # <- this is your full "10 tools" universe

# Align threshold tables to the master factor levels
data_list$abundance60 <- data_list$abundance60 %>%
  dplyr::mutate(
    sample = factor(sample, levels = all_samples),
    tool   = factor(tool,   levels = all_tools)
  )

data_list$abundance70 <- data_list$abundance70 %>%
  dplyr::mutate(
    sample = factor(sample, levels = all_samples),
    tool   = factor(tool,   levels = all_tools)
  )

data_list$abundance80 <- data_list$abundance80 %>%
  dplyr::mutate(
    sample = factor(sample, levels = all_samples),
    tool   = factor(tool,   levels = all_tools)
  )

# --- Build a fixed sample x tool grid (forces all tools to exist, even if absent in data) ---
sample_tool_grid <- tidyr::expand_grid(
  sample = all_samples,
  tool   = all_tools
)

# --- Precompute summaries (matches original behavior, but faster + stable tool counts) ---
prep_abundance <- function(df) {
  
  # Map habitat/habitat2 from sample (like your original left_join distinct(sample, habitat...))
  habitat_map <- df %>%
    dplyr::group_by(sample) %>%
    dplyr::summarise(
      habitat  = dplyr::first(habitat[!is.na(habitat)]),
      habitat2 = dplyr::first(habitat2[!is.na(habitat2)]),
      .groups = "drop"
    )
  
  # Aggregate abundance + diversity per tool x sample
  summ <- df %>%
    dplyr::group_by(tool, sample) %>%
    dplyr::summarise(
      normed10m = sum(normed10m, na.rm = TRUE),
      unigenes  = sum(unigenes,  na.rm = TRUE),
      .groups = "drop"
    )
  
  sample_tool_grid %>%
    dplyr::left_join(summ, by = c("sample", "tool")) %>%
    dplyr::left_join(habitat_map, by = "sample") %>%
    dplyr::mutate(
      normed10m = tidyr::replace_na(normed10m, 0),
      unigenes  = tidyr::replace_na(unigenes, 0)
    ) %>%
    dplyr::arrange(tool, sample)
}

# --- Build prepped datasets ---
# Default branch in your original code did NOT exclude DeepARG/RGI-DIAMOND
abundance_base_all <- data_list$abundance

# Threshold branches in your original code excluded DeepARG/RGI-DIAMOND before binding extras
abundance_base_excl <- data_list$abundance %>%
  dplyr::filter(!tool %in% abundance_tools_excl)

abundance_prepped <- list(
  "default" = prep_abundance(abundance_base_all),
  "60"      = prep_abundance(dplyr::bind_rows(abundance_base_excl, data_list$abundance60)),
  "70"      = prep_abundance(dplyr::bind_rows(abundance_base_excl, data_list$abundance70)),
  "80"      = prep_abundance(dplyr::bind_rows(abundance_base_excl, data_list$abundance80))
)

dplyr::count(abundance_prepped[["default"]], tool) %>% print(n = Inf)
dplyr::count(abundance_prepped[["60"]], tool) %>% print(n = Inf)


#For median abundance and diversity plot
pal_for_tools <- function(selected_tools, tools_levels, pal_10_q) {
  sel <- intersect(tools_levels, selected_tools)
  cols <- pal_10_q[match(sel, tools_levels)]
  stats::setNames(cols, sel)   # named vector is key
}

