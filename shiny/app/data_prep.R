library(tidyverse)

DATA_DIR <- "code_R_analysis/output_abundance_diversity_resistome"
PREP_DIR <- "data" 

#All Habitats
EN <- c("human gut", "human oral", "human skin", "human nose", "human vagina", 
        "dog gut", "cat gut", "mouse gut", "pig gut", "wastewater", "marine", 
        "freshwater", "soil", "amplicon", "isolate", "built-environment")

#Source for each habitats
SO <- c(rep("humans", 5), rep("mammals", 4), "wastewater", "marine", "freshwater", 
        "soil", rep("other", 2), "built-environment")

names(SO) <- EN
not_env <- c("amplicon", "isolate", "built-environment")
EN2 <- EN[!EN %in% not_env]
h2 <- c("humans","mammals","wastewater", "freshwater", "soil", "marine")
tools_levels <- c("DeepARG", "fARGene", "ABRicate-ARGANNOT", "ABRicate-MEGARes", 
                  "RGI-DIAMOND", "ABRicate-CARD", "AMRFinderPlus", 
                  "ABRicate-NCBI", "ResFinder", "ABRicate-ResFinder")

# Sourcing helper.R for create_class_overlaps
source("code_R_analysis/helper.R")

data_list <- list()

# Load Metadata
metadata <- read.delim("data/metadata_GMGC10.sample.meta.tsv")
data_list$metadata <- metadata

# Load and Prep Results Tools (Unigenes & Recall/FNR)
lst <- readRDS(file.path(DATA_DIR, "results_tools.rds"))
unigenes <- as_tibble(do.call(rbind, lapply(lst, function(x) x[, c("query","tool","ARO","parent","parent_description","new_level","id")]))) %>%
  filter(tool %in% tools_levels) %>%
  mutate(tool = factor(tool, levels = tools_levels))

data_list$levels_unigenes <- unigenes %>%
  group_by(query) %>% slice_head(n = 1) %>% ungroup() %>%
  group_by(new_level) %>% summarise(n = n(), .groups = "drop") %>%
  arrange(desc(n)) %>% pull(new_level)

data_list$unigenes_base <- unigenes
data_list$recall_fnr    <- create_class_overlaps(unigenes)
data_list$recall_fnr60  <- create_class_overlaps(unigenes %>% filter(!(tool %in% c("DeepARG","RGI-DIAMOND") & id < 60)))
data_list$recall_fnr70  <- create_class_overlaps(unigenes %>% filter(!(tool %in% c("DeepARG","RGI-DIAMOND") & id < 70)))
data_list$recall_fnr80  <- create_class_overlaps(unigenes %>% filter(!(tool %in% c("DeepARG","RGI-DIAMOND") & id < 80)))

# Load and Prep Abundances
process_abundance_file <- function(file_paths) {
  df <- bind_rows(lapply(file_paths, function(x) readRDS(file.path(DATA_DIR, x)))) %>%
    mutate(unigenes = distinct_unigenes_rarefied,
           habitat = factor(habitat, levels = EN),
           habitat2 = factor(SO[habitat], levels = h2),
           tool = factor(tool, levels = tools_levels)) %>%
    filter(tool %in% tools_levels & !habitat %in% not_env, aggregation == "new_level")
  return(df)
}

ab_base <- process_abundance_file(c("abundance_diversity_part1.rds", "abundance_diversity_part2.rds", "abundance_diversity_part3.rds"))
ab_60   <- process_abundance_file("abundance_diversity_60.rds")
ab_70   <- process_abundance_file("abundance_diversity_70.rds")
ab_80   <- process_abundance_file("abundance_diversity_80.rds")

# Pre-bind the excluded tool logic for thresholds offline
ab_base_excl <- ab_base %>% filter(!tool %in% c("DeepARG", "RGI-DIAMOND"))

# Function to build the sample_tool grids 
prep_abundance_grid <- function(df) {
  df <- df %>% mutate(sample = factor(sample), tool = factor(tool, levels = tools_levels))
  all_samples <- levels(df$sample)
  all_tools   <- levels(df$tool)
  grid <- tidyr::expand_grid(sample = all_samples, tool = all_tools)
  
  habitat_map <- df %>% group_by(sample) %>% summarise(habitat = first(habitat[!is.na(habitat)]), .groups = "drop")
  
  summ <- df %>% group_by(tool, sample) %>% summarise(normed10m = sum(normed10m, na.rm = TRUE), unigenes = sum(unigenes, na.rm = TRUE), .groups = "drop")
  
  grid %>%
    left_join(summ, by = c("sample", "tool")) %>%
    left_join(habitat_map, by = "sample") %>%
    mutate(normed10m = replace_na(normed10m, 0), unigenes = replace_na(unigenes, 0)) %>%
    arrange(tool, sample)
}

data_list$abundance_prepped <- list(
  "default" = prep_abundance_grid(ab_base),
  "60"      = prep_abundance_grid(bind_rows(ab_base_excl, ab_60)),
  "70"      = prep_abundance_grid(bind_rows(ab_base_excl, ab_70)),
  "80"      = prep_abundance_grid(bind_rows(ab_base_excl, ab_80))
)

# Similar pre-binding for abundance_class
process_ab_class <- function(df) {
  df %>% ungroup() %>%
    mutate(gene = factor(gene), sample = factor(sample)) %>%
    select(sample, gene, tool, normed10m, unigenes) %>%
    complete(sample, gene, tool, fill = list(normed10m = 0, unigenes = 0)) %>%
    mutate(
      aggregation = "new_level",
      gene = as.character(gene), # <--- This prevents the 'Other' genes from being dropped!
      habitat = metadata$habitat[match(sample, metadata$sample_id)],
      habitat2 = metadata$habitat2[match(sample, metadata$sample_id)],
      location = metadata$location[match(sample, metadata$sample_id)]
    )
}

class_base <- process_ab_class(ab_base)
class_excl <- class_base %>% filter(!tool %in% c("DeepARG", "RGI-DIAMOND"))

data_list$abundance_class_prepped <- list(
  "default" = class_base,
  "60"      = bind_rows(class_excl, process_ab_class(ab_60)),
  "70"      = bind_rows(class_excl, process_ab_class(ab_70)),
  "80"      = bind_rows(class_excl, process_ab_class(ab_80))
)


# Load and Prep Pan/Core Resistome

# Helper function to process core files
process_core_file <- function(file_name) {
  readRDS(file.path(DATA_DIR, file_name)) %>%
    rename(new_level = new_level_centroid,
           X = centroid) %>%
    filter(tool %in% tools_levels,
           !habitat %in% not_env) %>%
    mutate(habitat = factor(habitat, levels = EN2),
           tool = factor(tool, levels = tools_levels))
}

# Helper function to process pan files
process_pan_file <- function(file_name) {
  readRDS(file.path(DATA_DIR, file_name)) %>%
    filter(tool %in% tools_levels,
           !habitat %in% not_env,
           aggregation %in% "new_level_centroid") %>%
    mutate(habitat = factor(habitat, levels = EN2),
           tool = factor(tool, levels = tools_levels))
}

# Helper function to calculate sumpan2 from pan data
calculate_sumpan2 <- function(pan_df) {
  pan_df %>% 
    ungroup() %>%
    group_by(tool, habitat, aggregation, epoch) %>%
    summarise(s = sum(unigenes), .groups = "drop") %>%
    group_by(tool, habitat, aggregation) %>%
    summarise(md = median(s), mn = mean(s), sd = sd(s), .groups = "drop")
}

# Load Base Files
core_base <- process_core_file("core_resistome.rds")
pan_base  <- process_pan_file("pan_resistome.rds")
sumpan2_base <- calculate_sumpan2(pan_base)

# Load Threshold 60 
core_60 <- process_core_file("core_resistome60.rds")
pan_60  <- process_pan_file("pan_resistome60.rds")
sumpan2_60 <- calculate_sumpan2(pan_60)

# Load Threshold 70
core_70 <- process_core_file("core_resistome70.rds")
pan_70  <- process_pan_file("pan_resistome70.rds")
sumpan2_70 <- calculate_sumpan2(pan_70)

# Load Threshold 80 
core_80 <- process_core_file("core_resistome80.rds")
pan_80  <- process_pan_file("pan_resistome80.rds")
sumpan2_80 <- calculate_sumpan2(pan_80)

# Pre-bind Excluded Tool Logic
tools_excl <- c("DeepARG", "RGI-DIAMOND")

core_base_excl <- core_base %>% filter(!tool %in% tools_excl)
pan_base_excl  <- pan_base %>% filter(!tool %in% tools_excl)
sumpan2_base_excl <- sumpan2_base %>% filter(!tool %in% tools_excl)

# Save Pre-bound Lists
data_list$core_prepped <- list(
  "default" = core_base,
  "60"      = bind_rows(core_base_excl, core_60),
  "70"      = bind_rows(core_base_excl, core_70),
  "80"      = bind_rows(core_base_excl, core_80)
)

data_list$pan_prepped <- list(
  "default" = pan_base,
  "60"      = bind_rows(pan_base_excl, pan_60),
  "70"      = bind_rows(pan_base_excl, pan_70),
  "80"      = bind_rows(pan_base_excl, pan_80)
)

data_list$sumpan2_prepped <- list(
  "default" = sumpan2_base,
  "60"      = bind_rows(sumpan2_base_excl, sumpan2_60),
  "70"      = bind_rows(sumpan2_base_excl, sumpan2_70),
  "80"      = bind_rows(sumpan2_base_excl, sumpan2_80)
)

# Saving using rds
# saveRDS(data_list, file.path(PREP_DIR, "precomputed_data.rds"))
qs::qsave(data_list, file.path(PREP_DIR, "precomputed_data.qs"))