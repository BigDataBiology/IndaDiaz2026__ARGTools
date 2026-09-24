library(dplyr)
library(ggplot2)
library(gridExtra)
library(tidyverse)
library(RColorBrewer)
library(ggpattern)
library(grid)
library(Cairo)
library(cowplot)
library(scales)
library(ggbreak)

setwd("~/Documents/GitHub/arg_compare/")
options(dplyr.summarise.inform = FALSE)
source("code_R_analysis/helper.R")

# Sourced gene classes 

general_size <- 6
lab_fn <- function(x) {
  x <- gsub("-", "-\n", x)
  x <- gsub(" ", "\n", x)
  x <- gsub("/", "/\n", x)
  x
}

# FORMAT PLOTS 
pal_7 <- brewer.pal(8, "Dark2")
pal_7 <- pal_7[-7]
pal_7 <- pal_7[c(1,2,3,4,6,5,7)]


pal_10_q <- pal_7[c(1,2,3,4,5,5,6,6,7,7)]
pal_10_complete <- brewer.pal(7, "Dark2")
pal_10_complete <- pal_7

# pattern for plots

pattern_density <- 0.001 
pattern_spacing <- 0.025
pattern_fill <- "black"
pattern_size <- 0.12

# all HABITATS
EN <- c("human gut", "human oral",  "human skin", 
        "human nose", "human vagina", 
        "dog gut", "cat gut", "mouse gut", 
        "pig gut", "wastewater", "marine", 
        "freshwater", "soil" )

# SOURCE FOR EACH HABITAT

SO <- c(rep("humans", 5), rep("mammals", 4),  
        "wastewater", "marine", "freshwater", "soil")

names(SO) <- EN

# tools compared in the paper 

basic_tools <- c(
  "DeepARG", "fARGene","ABRicate-ARGANNOT", "ABRicate-MEGARes",
  "RGI-DIAMOND", "ABRicate-CARD","AMRFinderPlus", "ABRicate-NCBI",
  "ResFinder", "ABRicate-ResFinder")

# all tool levels
tools_levels <- c(
  "DeepARG", "fARGene","ABRicate-ARGANNOT", "ABRicate-MEGARes",
  "RGI-DIAMOND", "ABRicate-CARD","AMRFinderPlus", "ABRicate-NCBI",
  "ResFinder", "ABRicate-ResFinder","DeepARG70","DeepARG80","DeepARG90",
  "RGI-DIAMOND70","RGI-DIAMOND80","RGI-DIAMOND90",
  "DeepARG-aa", "RGI-BLAST", "RGI-DIAMOND-aa", "fARGene-aa", "AMRFinderPlus-nt")

# tag reference genes 

tag_to_tool <- c(
  "resfinder"          = "ResFinder",
  "abricate-resfinder" = "ABRicate-ResFinder",
  "abricate-argannot"  = "ABRicate-ARGANNOT",
  "abricate-megares"   = "ABRicate-MEGARes",
  "abricate-card"      = "ABRicate-CARD",
  "abricate-ncbi"      = "ABRicate-NCBI",
  "deeparg"            = "DeepARG",
  "rgi-card"           = "RGI-DIAMOND",
  "amrfinderplus"      = "AMRFinderPlus"
)


# add the name of each tool to the color palet 
names(pal_10_q) <- basic_tools

# repeat the color for the different thresholds in deeparg and rgi

da <- rep(pal_10_q["DeepARG"],3)
names(da) <- c("DeepARG70","DeepARG80","DeepARG90")
rgi <- rep(pal_10_q["RGI-DIAMOND"],3)
names(rgi) <- c("RGI-DIAMOND70","RGI-DIAMOND80","RGI-DIAMOND90")
other <- c(pal_10_q["DeepARG"], pal_10_q["RGI-DIAMOND"],
           pal_10_q["RGI-DIAMOND"],pal_10_q["fARGene"],pal_10_q["AMRFinderPlus"])
names(other) <- c("DeepARG-aa", "RGI-BLAST", 
                  "RGI-DIAMOND-aa", "fARGene-aa","AMRFinderPlus-nt")

pal_10_q <- c(pal_10_q, da, rgi, other)
rm(da, rgi, other)

# The name of eaach pipeline in the plots
tools_labels <- c(
  "DeepARG", "fARGene", "ABRicate-\nARGANNOT", "ABRicate-\nMEGARes",
  "RGI", "ABRicate-\nCARD", "AMRFinder-\nPlus", "ABRicate-\nNCBI",
  "ResFinder", "ABRicate-\nResFinder",
  "DeepARG-70%","DeepARG-80%","DeepARG-90%","RGI-70%","RGI-80%","RGI-90%",
  "DeepARG-aa", "RGI/nBLAST", "RGI-aa", "fARGene-aa", "AMRFinder-\nPlus-nt")

names(tools_labels) <- tools_levels

# add factor for ordering the tools within plots 

tools_labels_factor <- c(
  "DeepARG", "fARGene", "ABRicate-\nARGANNOT", "ABRicate-\nMEGARes", 
  "RGI", "ABRicate-\nCARD", "AMRFinder-\nPlus", "ABRicate-\nNCBI",
  "ResFinder", "ABRicate-\nResFinder", "DeepARG-70%","DeepARG-80%","DeepARG-90%",
  "RGI-70%","RGI-80%","RGI-90%",
  "DeepARG-aa", "RGI/nBLAST", "RGI-aa", "fARGene-aa", "AMRFinder-\nPlus-nt")

# one space for DeepARG
# two spaces for fARGene

tools_db <- c(" ", "  ", "   ", "    ", "CARD", "CARD","NCBI","NCBI", 
              "ResFinder","ResFinder"," "," "," ","CARD","CARD","CARD",
              " ", "CARD", "CARD", "  ", "NCBI")

# 
# two spaces for fARGene, one space for DeepARG
tools_db_factor <- c(" ", "  ", "   ", "    ",
                     "CARD", "NCBI", "ResFinder")

tools_texture <- c("ABRicate-CARD", "ABRicate-NCBI", "ABRicate-ResFinder")

ARO <- read.csv("code_R_analysis/output_abundance_diversity_resistome/conversion_ARO_parent_new_level.csv")

unigenes <- tibble(readRDS(file = "code_R_analysis/output_abundance_diversity_resistome/unigenes_per_tool.rds")) %>% 
  # convert MFS to efflux pump
  mutate(gene_class = ifelse(new_level == "MFS efflux pump", "efflux pump", new_level)) %>%  
  mutate(tool = factor(tool, levels = tools_levels)) %>%
  mutate(tools_labels = factor(tools_labels[tool], levels = tools_labels_factor),
         texture = ifelse(tool %in% tools_texture, "yes", "no"),
         tools_db = factor(tools_db[tool], levels = tools_db_factor))


aros_per_unigene <- unigenes %>% filter(tool %in% basic_tools) %>% 
  group_by(query) %>% 
  summarise(n = n_distinct(ARO)) %>% 
  ungroup() %>% 
  group_by(n) %>% 
  summarise(aro_per_unigene = n()) %>% 
  mutate(p = round(100*aro_per_unigene/sum(aro_per_unigene),2))

library(widyr)

aro_pairs <- unigenes %>% 
  filter(tool %in% basic_tools) %>% 
  group_by(query) %>% 
  filter(n_distinct(ARO) > 1) %>%
  ungroup() %>%
  pairwise_count(ARO, query, sort = TRUE, upper = FALSE)

aro_pairs <- aro_pairs %>% mutate(description_1 = ARO$Term_Label[match(item1, ARO$Term_ID)],
                     description_2 = ARO$Term_Label[match(item2, ARO$Term_ID)],
                     parent_description_1 = ARO$Parent_Label[match(item1, ARO$Term_ID)],
                     parent_description_2 = ARO$Parent_Label[match(item2, ARO$Term_ID)])

aro_pairs_mismatch <- aro_pairs %>% 
  mutate(term1 = pmin(description_1, description_2),
         term2 = pmax(description_1, description_2)) %>% 
  group_by(term1, term2) %>%
  summarise(ARO1 = item1[1], ARO2 = item2[1], total_n = sum(n), .groups = "drop") %>%
  arrange(desc(total_n)) %>% 
  rename(n_unigenes = total_n) %>%
  select(ARO1, term1, ARO2, term2, n_unigenes) %>% 
  mutate(gene_class1 = ARO$new_level[match(ARO1, ARO$Term_ID)],
         gene_class2 = ARO$new_level[match(ARO2, ARO$Term_ID)])

gene_class_pairs_mismatch <- aro_pairs_mismatch %>% 
  filter(gene_class1 != gene_class2)

unigenes %>% 
       filter(tool %in% basic_tools) %>% 
       group_by(query) %>% 
       filter(n_distinct(ARO) > 1) %>% ungroup() %>% summarise(n_distinct(query))
# A tibble: 1 × 1
#`n_distinct(query)`
#<int>
#  1                9542

unigenes %>% 
  filter(tool %in% basic_tools) %>% 
  group_by(query) %>% 
  summarise(n_aro = n_distinct(ARO) ) %>% ungroup() %>% 
  group_by(n_aro) %>% summarise(n = n())

# # A tibble: 6 × 2
# n_aro      n
# <int>  <int>
#   1     1 169053
# 2     2   8937
# 3     3    515
# 4     4     66
# 5     5     18
# 6     6      6

write.csv(aro_pairs_mismatch, "arg_norm_correction/number_times_aro_mismatch_per_ARO_pair.csv")
sum(aro_pairs_mismatch$n_unigenes)
sum(gene_class_pairs_mismatch$n_unigenes)



risk_category <- table(paste("aro",unigenes$rank_aro[unigenes$tool %in% basic_tools]), 
                       paste("bit",unigenes$rank_highest_bit_80[unigenes$tool %in% basic_tools]))

overlap_aro_bitscore_risk <- unigenes %>% group_by(rank_aro) %>% 
  filter(tool %in% basic_tools, rank_aro %in% c("I","II","III","IV")) %>%
  mutate(same_as_bitscore_blast = (rank_aro == rank_highest_bit_80)) %>% 
  summarise(n_by_aro = n(), overlap = sum(same_as_bitscore_blast)) %>% bind_cols(

unigenes %>% group_by(rank_highest_bit_80) %>% 
  filter(tool %in% basic_tools, rank_highest_bit_80 %in% c("I","II","III","IV")) %>%
  mutate(same_as_bitscore_blast = (rank_aro == rank_highest_bit_80)) %>% 
  summarise(n_by_alignment = n(), same_as_aro = sum(same_as_bitscore_blast)) %>%
  select(n_by_alignment)) %>% 
  rename(rank = rank_aro) %>% select(rank, n_by_aro, n_by_alignment, overlap)
overlap_aro_bitscore_risk


core_levels   <- c("I", "II", "III", "IV")
margin_levels <- c("Not found", "notassessed")
total_level   <- "Total assessed"
all_levels    <- c(core_levels, total_level, margin_levels)

df <- as.data.frame(as.table(risk_category)) %>%
  rename(aro = Var1, bit = Var2, n = Freq) %>%
  mutate(
    aro = gsub("^aro ", "", aro),
    bit = gsub("^bit ", "", bit)
  )

# "Total assessed" row: sum across aro I–IV, per existing bit column
row_totals <- df %>%
  filter(aro %in% core_levels) %>%
  group_by(bit) %>%
  summarise(n = sum(n), .groups = "drop") %>%
  mutate(aro = total_level)

# "Total assessed" column: sum across bit I–IV, per existing aro row
col_totals <- df %>%
  filter(bit %in% core_levels) %>%
  group_by(aro) %>%
  summarise(n = sum(n), .groups = "drop") %>%
  mutate(bit = total_level)

# corner cell: grand total of the I–IV core block
corner <- df %>%
  filter(aro %in% core_levels, bit %in% core_levels) %>%
  summarise(n = sum(n)) %>%
  mutate(aro = total_level, bit = total_level)

df <- bind_rows(df, row_totals, col_totals, corner) %>%
  mutate(
    aro     = factor(aro, levels = all_levels),
    bit     = factor(bit, levels = all_levels),
    is_core = aro %in% core_levels & bit %in% core_levels
  )

core_max <- df %>% filter(is_core) %>% summarise(m = max(n)) %>% pull(m)

df <- df %>%
  mutate(
    text_color = case_when(
      !is_core            ~ "grey20",
      n > core_max * 0.55 ~ "white",
      TRUE                ~ "grey20"
    )
  )

risk_ambiguity <- ggplot(df, aes(x = bit, y = aro)) +
  geom_tile(data = filter(df, !is_core),
            fill = "grey95", color = "white", linewidth = 0.6) +
  geom_tile(data = filter(df, is_core),
            aes(fill = n), color = "white", linewidth = 0.6) +
  geom_text(aes(label = comma(n), color = I(text_color)), size = 3.2) +
  scale_fill_viridis_c(option = "mako", direction = -1,
                       limits = c(0, core_max), name = "Count\n(I–IV scale)") +
  scale_x_discrete(position = "top") +
  scale_y_discrete(limits = rev(all_levels)) +
  labs(
    x = "Highest bit-score alignment", y = "ARO matching risk",
    title = "Risk category cross-tabulation",
    subtitle = "Color intensity scaled to the I\u2013IV block only \u2014 totals and margins shown as counts"
  ) +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank(), axis.ticks = element_blank())


ambiguity_blast_unigenes_risk <- unigenes %>% 
  filter(tool %in% basic_tools) %>% 
  group_by(query) %>% 
  summarise(n_rank_aro = n_distinct(rank_aro),
            n_rank_highest_bit_80 = n_distinct(rank_highest_bit_80)) %>% 
  ungroup() %>% 
  group_by(n_rank_highest_bit_80) %>% 
  summarise(n = n()) %>%
  mutate(p = round(n / sum(n)*100, 2))
ambiguity_blast_unigenes_risk


ambiguity_aro_unigenes_risk <- unigenes %>% 
  filter(tool %in% basic_tools) %>% 
  group_by(query) %>% 
  summarise(n_rank_aro = n_distinct(rank_aro),
            n_rank_highest_bit_80 = n_distinct(rank_highest_bit_80)) %>% 
  ungroup() %>% 
  group_by(n_rank_aro) %>% 
  summarise(n = n()) %>%
  mutate(p = round(n / sum(n)*100, 2))
ambiguity_aro_unigenes_risk


risk_ambiguity
