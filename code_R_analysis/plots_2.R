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
library(ggnewscale)
library(MASS)
library(car)
library(patchwork)


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
pal_7 <- pal_7[-8]
pal_7 <- pal_7[c(1,2,3,4,6,5,7)]

pal_10_q <- pal_7[c(1,2,3,4,5,5,6,6,7,7)]
pal_10_q <- pal_7[c(5,5,6,6,7,7,3,4,2,1)]

pal_7 <- pal_7[c(5,6,7,3,4,2,1)]
#pal_10_complete <- brewer.pal(7, "Dark2")
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
  "RGI-DIAMOND", "ABRicate-CARD","AMRFinderPlus", "ABRicate-NCBI",
  "ResFinder", "ABRicate-ResFinder", 
  "ABRicate-ARGANNOT", "ABRicate-MEGARes", "fARGene", "DeepARG" )

basic_tools <- c(
  "RGI-DIAMOND", "ABRicate-CARD","AMRFinderPlus", "ABRicate-NCBI",
  "ResFinder", "ABRicate-ResFinder",
  "ABRicate-ARGANNOT", "ABRicate-MEGARes", "fARGene", "DeepARG")

# all tool levels
tools_levels <- c(
  "RGI-DIAMOND", "ABRicate-CARD","AMRFinderPlus", "ABRicate-NCBI",
  "ResFinder", "ABRicate-ResFinder", "ABRicate-ARGANNOT", "ABRicate-MEGARes", "fARGene", "DeepARG",
  "DeepARG70","DeepARG80","DeepARG90",
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
  "amrfinderplus"      = "AMRFinderPlus",
  "fargene" = "fARGene"
)


tool_map <- c(
  "abricate-argannot" = "ABRicate-ARGANNOT" ,
  "abricate-card"      = "ABRicate-CARD",
  "abricate-megares"   = "ABRicate-MEGARes",
  "abricate-ncbi"      =  "ABRicate-NCBI",
  "abricate-resfinder" = "ABRicate-ResFinder",
  "amrfinderplus"      = "AMRFinderPlus",
  "deeparg"            = "DeepARG",
  "resfinder"          = "ResFinder",
  "rgi-card"        = "RGI-DIAMOND",
  "fargene" = "fARGene"
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
  "RGI", "ABRicate-\nCARD", "AMRFinder-\nPlus", "ABRicate-\nNCBI",
  "ResFinder", "ABRicate-\nResFinder",
  "ABRicate-\nARGANNOT", "ABRicate-\nMEGARes", "fARGene", "DeepARG", 
  "DeepARG-70%","DeepARG-80%","DeepARG-90%","RGI-70%","RGI-80%","RGI-90%",
  "DeepARG-aa", "RGI/nBLAST", "RGI-aa", "fARGene-aa", "AMRFinder-\nPlus-nt")

names(tools_labels) <- tools_levels

# add factor for ordering the tools within plots 

tools_labels_factor <- c(
  "RGI", "ABRicate-\nCARD", "AMRFinder-\nPlus", "ABRicate-\nNCBI",
  "ResFinder", "ABRicate-\nResFinder",
  "ABRicate-\nARGANNOT", "ABRicate-\nMEGARes", "fARGene" , "DeepARG", 
   "DeepARG-70%","DeepARG-80%","DeepARG-90%",
  "RGI-70%","RGI-80%","RGI-90%",
  "DeepARG-aa", "RGI/nBLAST", "RGI-aa", "fARGene-aa", "AMRFinder-\nPlus-nt")

# one space for DeepARG
# two spaces for fARGene

tools_db <- c(  "CARD", "CARD","NCBI","NCBI", 
              "ResFinder","ResFinder","   ", "    ", "  ", " ",
              " "," "," ", "CARD","CARD","CARD",
              " ", "CARD", "CARD", "  ", "NCBI")

# 
# two spaces for fARGene, one space for DeepARG
tools_db_factor <- c("CARD", "NCBI", "ResFinder", "   ", "    ",
                    "  ", " ")

tools_texture <- c("ABRicate-CARD", "ABRicate-NCBI", "ABRicate-ResFinder")


##


################################################################################################

# vector to encode the texture of the bars - 
# striped for abricate with card, ncbi, and resfinder

tools_texture_code <- rep("none", length(tools_labels))
tools_texture_code[tools_levels %in% tools_texture] <- "stripe"

# base format for the plots 

theme1 <- theme_minimal() +
  theme(
    legend.position = "none",
    text = element_text(size = general_size, color = "black"),
    title = element_text(size = general_size + 2, face = "bold"),
    axis.title = element_text(size = general_size , face = "bold"),
    strip.text = element_text(size = general_size, angle = 0),
    panel.background = element_blank(),
    axis.text.x = element_text(size = general_size, angle = 90,
                               vjust = 0.5, hjust = 1),
    axis.text.y = element_text(size = general_size),
    panel.border = element_rect(color =  "black"),
    plot.margin = margin(0, 0, 0, 0, unit = "pt"),
    legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
    legend.margin = margin(0, 0, 0, 0, unit = "pt"),
    panel.spacing = unit(0, "pt"),
    legend.text = element_text(size = general_size),
    panel.grid.minor.y = element_blank(),
    legend.spacing.y = unit(0, "cm"),
    legend.spacing.x = unit(0.1, "cm"),
    legend.key.height = unit(0.35, "cm"),
    legend.key.width  = unit(0.35, "cm"),
    legend.box.spacing = unit(0.1, "cm"))



################################################################################
################################################################################
## DATASETS 
unigenes_per_habitat <-  read.delim("data/abundances/args_by_habitat.tsv")

# gene classes 

ARO <- read.csv("code_R_analysis/output_abundance_diversity_resistome/conversion_ARO_parent_new_level.csv")

# results per individual tool 

lst <- readRDS("code_R_analysis/output_abundance_diversity_resistome/results_tools_all_GMGC.rds")
lst0 <- readRDS("code_R_analysis/output_abundance_diversity_resistome/results_tools.rds")
lst <- lst[c(1,4,6,9,10,11,12,13,14,15)]

# metagenomes' metadata

metadata <- read.delim("data/metadata_GMGC10.sample.meta.tsv")
metadata <- metadata %>% filter(!habitat %in% c("amplicon", "isolate", "built-environment"))

# metadata0 will be merged with abundance 
metadata0 <- metadata %>% select(sample_id, habitat) %>% rename(sample = sample_id)

abundance <- readRDS("code_R_analysis/output_abundance_diversity_resistome/abundance_diversity.rds")
abundance <- abundance %>% filter(tool %in% basic_tools)

# convert MFS to efflux pump
abundance <- abundance %>% mutate(gene = ifelse(gene == "MFS efflux pump", "efflux pump", gene))

metadata0 <- metadata0 %>% left_join(abundance, by = "sample")

# we add one class to avoid NA values and to be able to "complete" for all classes later
# we add one tool to avoid NA values and to be able to "complete" for all tools later

metadata0$gene[is.na(metadata0$gene)] <- "efflux pump"
metadata0$tool[is.na(metadata0$tool)] <- "DeepARG"
metadata0$abundance[is.na(metadata0$abundance)] <- 0
metadata0$richness[is.na(metadata0$richness)] <- 0
metadata0$richness_no_rarified[is.na(metadata0$richness_no_rarified)] <- 0


abundance <- metadata0 %>%
  mutate(habitat = factor(habitat, levels = EN),
         tool = factor(tool, levels = tools_levels)) %>%
  mutate(tools_labels = factor(tools_labels[tool], levels = tools_labels_factor),
         texture = ifelse(tool %in% tools_texture, "yes", "no"),
         tools_db = factor(tools_db[tool], levels = tools_db_factor))


# total abundance per sample

abundance_tool_sample <- abundance %>% 
  ungroup() %>%
  group_by(tool, sample, habitat) %>%  
  summarise(abundance = sum(abundance), richness = sum(richness)) %>% 
  ungroup() %>% 
  complete(sample, tool) %>% # complete all samples for all tools with NAs
  left_join(abundance %>% select(sample, habitat) %>% 
              distinct(), by = "sample") %>% 
  mutate(habitat  = coalesce(habitat.x, habitat.y)) %>%
  select(-habitat.x, -habitat.y) %>% 
  mutate(abundance = replace_na(abundance, 0)) %>%  
  mutate(richness = replace_na(richness, 0)) %>% 
  arrange(tool, sample) %>%
  mutate(tools_labels = factor(tools_labels[tool], levels = tools_labels_factor),
         texture = ifelse(tool %in% tools_texture, "yes", "no"),
         tools_db = factor(tools_db[tool], levels = tools_db_factor))

# abundance per sample and gene class

abundance_class <- abundance %>% 
  ungroup() %>% 
  select(-c(richness_no_rarified)) %>%
  mutate(gene = factor(gene), sample = factor(sample), 
         tool = factor(tool, levels = tools_levels)) %>%
  complete( sample, gene, tool, 
            fill = list(abundance = 0,  richness = 0)) %>% 
  mutate(gene = as.character(gene),
         habitat = metadata$habitat[match(sample, metadata$sample_id)])  %>%
  mutate(tools_labels = factor(tools_labels[tool], levels = tools_labels_factor),
         texture = ifelse(tool %in% tools_texture, "yes", "no"),
         tools_db = factor(tools_db[tool], levels = tools_db_factor))

## CORE RESISTOME 

core <- readRDS(file = "code_R_analysis/output_abundance_diversity_resistome/core_resistome.rds") %>% 
  mutate(new_level = ifelse(new_level == "MFS efflux pump", "efflux pump", new_level)) %>% 
  filter(tool %in% basic_tools) %>% 
  filter(habitat %in% EN) %>% 
  mutate(habitat = factor(habitat, levels = EN), 
         tool = factor(tool, levels =  tools_levels)) %>% 
  mutate(tool = factor(tool, levels = tools_levels)) %>%
  mutate(tools_labels = factor(tools_labels[tool], levels = tools_labels_factor),
         texture = ifelse(tool %in% tools_texture, "yes", "no"),
         tools_db = factor(tools_db[tool], levels = tools_db_factor))


pan <- readRDS(file = "code_R_analysis/output_abundance_diversity_resistome/pan_resistome.rds") %>% 
  mutate(new_level = ifelse(new_level == "MFS efflux pump", "efflux pump", new_level)) %>% 
  filter(tool %in% basic_tools) %>% 
  filter(habitat %in% EN) %>% 
  mutate(habitat = factor(habitat, levels = EN), 
         tool = factor(tool, levels =  tools_levels)) %>% 
  mutate(tool = factor(tool, levels = tools_levels)) %>%
  mutate(tools_labels = factor(tools_labels[tool], levels = tools_labels_factor),
         texture = ifelse(tool %in% tools_texture, "yes", "no"),
         tools_db = factor(tools_db[tool], levels = tools_db_factor))

## unigenes identified by tool

unigenes <- tibble(readRDS(file = "code_R_analysis/output_abundance_diversity_resistome/unigenes_per_tool.rds")) %>% 
  # convert MFS to efflux pump
  mutate(gene_class = ifelse(new_level == "MFS efflux pump", "efflux pump", new_level)) %>%  
  mutate(tool = factor(tool, levels = tools_levels)) %>%
  mutate(tools_labels = factor(tools_labels[tool], levels = tools_labels_factor),
         texture = ifelse(tool %in% tools_texture, "yes", "no"),
         tools_db = factor(tools_db[tool], levels = tools_db_factor)) %>% 
  filter(tool %in% basic_tools)

unigenes0 <- tibble(readRDS(file = "code_R_analysis/output_abundance_diversity_resistome/unigenes_per_tool.rds")) %>% 
  # convert MFS to efflux pump
  mutate(gene_class = ifelse(new_level == "MFS efflux pump", "efflux pump", new_level)) %>%  
  mutate(tool = factor(tool, levels = tools_levels)) %>%
  mutate(tools_labels = factor(tools_labels[tool], levels = tools_labels_factor),
         texture = ifelse(tool %in% tools_texture, "yes", "no"),
         tools_db = factor(tools_db[tool], levels = tools_db_factor))

sum((unigenes0 %>% filter(!tool %in% basic_tools) %>% pull(query)) %in% 
(unigenes0 %>% filter(tool %in% basic_tools) %>% pull(query)))

sum(!(unigenes0 %>% filter(!tool %in% basic_tools) %>% pull(query)) %in% 
      (unigenes0 %>% filter(tool %in% basic_tools) %>% pull(query)))

## database_clusters
db_cluster <- read.delim("db_cluster/nested_with_coverage_out/nested_cluster_membership.tsv") 


db_cluster <- db_cluster %>% mutate(tool = sapply(strsplit(protein_id, split = "@@@"), function(x) x[1]),
                                    gene = sapply(strsplit(protein_id, split = "@@@"), function(x) x[2])) %>% 
  mutate(tool = factor(tool_map[tool], levels = basic_tools)) 


db_cluster <- db_cluster %>% mutate(tool = sapply(strsplit(protein_id, split = "@@@"), function(x) x[1]),
                                    gene = sapply(strsplit(protein_id, split = "@@@"), function(x) x[2]))

norm_prime <- function(x) {
  x <- gsub("###", "'", x, fixed = TRUE)
  x <- gsub("\u2019", "'", x, fixed = TRUE)
  x <- gsub("\u2032", "'", x, fixed = TRUE)
  trimws(x)
}

db_cluster <- db_cluster %>%
  mutate(gene_revised = ifelse(tool %in% "rgi-card", 
                               as.numeric(sub(".*\\|ID:([0-9]+)\\|.*", "\\1", gene)),
                               NA))

db_cluster <- db_cluster %>%
  mutate(gene_revised = ifelse(tool %in% c("abricate-argannot", "abricate-card","abricate-ncbi", "abricate-megares","abricate-resfinder"), 
                               sapply(strsplit(gene, split = "~~~", fixed = TRUE), `[`, 2),gene_revised)) %>%
  mutate(gene_revised = case_when(
    tool == "abricate-argannot"  ~ sub("\\([^)]*\\)", "", gene_revised),
    tool == "abricate-resfinder" ~ sub("_[0-9]+$", "", gene_revised),
    TRUE ~ gene_revised
  ))

extract_resfinder_accession <- function(x) {
  refseq_style <- str_extract(x, "[A-Z]{1,4}_[0-9A-Za-z.]+$")
  ifelse(!is.na(refseq_style), refseq_style, sub(".*_", "", x))
}

db_cluster <- db_cluster %>%
  mutate(gene = trimws(gene)) %>%
  mutate(gene_revised = ifelse(tool %in% c("resfinder"), extract_resfinder_accession(gene), gene_revised)) %>%
  mutate(gene_revised = ifelse(tool == "resfinder" & grepl("^fosA", gene), sub("_.*", "", gene), gene_revised)) %>% 
  mutate(gene_revised = ifelse(tool %in% c("resfinder") & gene_revised == "Y15705", "1_Y15705", gene_revised)) %>%
  mutate(gene_revised = ifelse(tool %in% c("resfinder") & gene_revised == "Y15704", "1_Y15704", gene_revised))

db_cluster <- db_cluster %>%
  mutate(gene_revised = ifelse(tool %in% c("resfinder"), trimws(extract_resfinder_accession(gene)), gene_revised))

db_cluster$gene_revised[db_cluster$tool == "resfinder" & db_cluster$gene_revised == "Y15704"] <- "1_Y15704"
db_cluster$gene_revised[db_cluster$tool == "resfinder" & db_cluster$gene_revised == "Y15705"] <- "1_Y15705"

db_cluster <- db_cluster %>%
  mutate(gene_revised = ifelse(tool %in% "amrfinderplus", 
                               sapply(strsplit(gene, split = "|", fixed = TRUE), function(x) x[2]),
                               gene_revised))

db_cluster <- db_cluster %>%
  mutate(gene_revised = ifelse(tool %in% "deeparg", gene, gene_revised))

db_cluster <- db_cluster %>% mutate(gene_revised = ifelse(tool %in% "fargene", norm_prime(gene), gene_revised))

db_cluster <- db_cluster %>% bind_rows(db_cluster %>% filter(tool %in% "deeparg") %>% mutate(tool = "deeparg-aa"))
db_cluster <- db_cluster %>% bind_rows(db_cluster %>% filter(tool %in% "fargene") %>% mutate(tool = "fargene.prot"))
db_cluster <- db_cluster %>% bind_rows(db_cluster %>% filter(tool %in% "rgi-card") %>% mutate(tool = "rgi-blast"))
db_cluster <- db_cluster %>% bind_rows(db_cluster %>% filter(tool %in% "rgi-card") %>% mutate(tool = "rgi-diamond-prot-card"))
db_cluster <- db_cluster %>% bind_rows(db_cluster %>% filter(tool %in% "amrfinderplus") %>% mutate(tool = "amrfinderplus.nt"))

db_cluster0 <- db_cluster %>% 
  filter(tool %in% c("resfinder","abricate-resfinder","abricate-argannot",
                     "abricate-megares", "abricate-card", "abricate-ncbi",
                     "deeparg", "fargene", "amrfinderplus", "rgi-card"))



# Overlap between tools by gene class
# per class and tool

recall_fnr <- create_class_overlaps(unigenes %>% filter(tool %in% basic_tools))
#recall_fnr_db <- create_class_overlaps(unigenes %>% filter(tool %in% basic_tools))

# overlap by pipeline (without considering classes)

JI_all <- return_overlap_tools(unigenes %>% 
                                 filter(tool %in% basic_tools))

JI_all_no_rpob_van <- return_overlap_tools(unigenes %>% 
                                                    filter(tool %in% basic_tools, 
                                                           !new_level %in% c("rpoB","van")))

JI_all_no_rpob_van_efflux <- return_overlap_tools(unigenes %>% 
                                 filter(tool %in% basic_tools, 
                                        !new_level %in% c("rpoB","van","efflux pump")))

JI_all_no_rpob_van %>% 
  filter(as.numeric(tool_ref) < as.numeric(tool_comp)) %>% 
  summarise(mn = mean(jaccard))

JI_all_no_rpob_van_efflux %>% 
  filter(as.numeric(tool_ref) < as.numeric(tool_comp)) %>% 
  summarise(mn = mean(jaccard))


JI_db_99_90 <- return_overlap_tools(db_cluster %>% filter(tool %in% basic_tools) %>% 
                                          rename(query = cluster_99_cov90) %>% 
                                          distinct(tool, query))
JI_db_95_90  <- return_overlap_tools(db_cluster %>% filter(tool %in% basic_tools) %>% 
                                            rename(query = cluster_95_cov90) %>% 
                                            distinct(tool, query))
JI_db_90_90  <- return_overlap_tools(db_cluster %>% filter(tool %in% basic_tools) %>% 
                                            rename(query = cluster_90_cov90) %>% 
                                            distinct(tool, query))

JI_db_99_60 <- return_overlap_tools(db_cluster %>% filter(tool %in% basic_tools) %>% 
                                           rename(query = cluster_99_cov60) %>% 
                                           distinct(tool, query))
JI_db_95_60  <- return_overlap_tools(db_cluster %>% filter(tool %in% basic_tools) %>% 
                                            rename(query = cluster_95_cov60) %>% 
                                            distinct(tool, query))
JI_db_90_60  <- return_overlap_tools(db_cluster %>% filter(tool %in% basic_tools) %>% 
                                            rename(query = cluster_90_cov60) %>% 
                                            distinct(tool, query))


#JI_db_99_90
top_abundance_JI <- c("class A beta-lactamase",
                      "class B beta-lactamase",
                      "class C beta-lactamase",
                      "class D beta-lactamase",
                      "erm", "mph", "aac", "aph", "ant",
                      "tet RPG", "tet enzyme", "qnr")

JI_all_some_gene_classes <- 
  return_overlap_tools(unigenes %>% 
                         filter(tool %in% basic_tools,
                         new_level %in% top_abundance_JI)) %>% 
  mutate(rank = "Subset gene classes") 

JI_all_rank_I_aro <- 
  return_overlap_tools(unigenes %>% 
                         filter(tool %in% basic_tools,
                         rank_aro == "I")) %>% 
  mutate(rank = "Risk I") 

JI_all_rank_II_aro <- return_overlap_tools(unigenes %>% 
                                            filter(tool %in% basic_tools,
                                                   rank_aro == "II")) %>% 
  mutate(rank = "Risk II") 

JI_all_rank_III_aro <- return_overlap_tools(unigenes %>% 
                                            filter(tool %in% basic_tools,
                                                   rank_aro == "III")) %>% 
  mutate(rank = "Risk III") 

JI_all_rank_IV_aro <- 
  return_overlap_tools(unigenes %>% 
  filter(tool %in% basic_tools, rank_aro == "IV")) 



#

all_jaccard <- JI_all %>% 
  mutate(rank = "Reported ARGs") %>% 
  bind_rows(JI_all_rank_I_aro %>% 
              mutate(rank = "Risk I"), 
            JI_all_rank_II_aro%>% 
              mutate(rank = "Risk II"),  
            JI_all_rank_III_aro %>% 
              mutate(rank = "Risk III"),  
            JI_all_rank_IV_aro %>% 
              mutate(rank = "Risk IV")) %>% 
  bind_rows(JI_all_some_gene_classes %>% 
              mutate(rank = "Subset gene classes")) %>% 
  bind_rows(JI_db_99_90 %>% 
              mutate(rank = "Reference database (ID99%, COV90%)")) %>% 
  bind_rows(JI_db_95_90 %>% 
              mutate(rank = "Reference database (ID95%, COV90%)")) %>% 
  bind_rows(JI_db_90_90  %>% 
              mutate(rank = "Reference database (ID90%, COV90%)"))


all_jaccard %>% 
  group_by(rank) %>% 
  summarise(n = n(), mean = mean(jaccard), median = median(jaccard))

all_jaccard %>%
  mutate(rank = sapply(strsplit(as.character(rank), split = "sk"), function(x) x[1])) %>%
  group_by(rank) %>%
  summarise(n = n(), mean = mean(jaccard), median = median(jaccard))


all_jaccard <- all_jaccard %>% 
  mutate(rank_label = factor(rank)) %>%
  mutate(tool_ref = factor(tools_labels[tool_ref], levels = tools_labels_factor[1:10]),
         tool_comp = factor(tools_labels[tool_comp], levels = tools_labels_factor[1:10]))

all_jaccard2 <- all_jaccard %>%
  group_by(rank) %>%
  mutate(
    rank_label = sprintf(
      "%s (mean = %s, median = %s)",
      rank,
      percent(mean(jaccard), accuracy = 0.1),
      percent(median(jaccard), accuracy = 0.1))) %>%
  ungroup() 

all_jaccard2 <- all_jaccard2 %>% 
  bind_rows(all_jaccard2 %>% rename(tool_comp = tool_ref, tool_ref = tool_comp)) %>% 
  filter(as.numeric(tool_ref) > as.numeric(tool_comp)) 

  
all_jaccard_plot <- ggplot(all_jaccard2, aes(x = tool_comp, y = tool_ref, fill = jaccard)) +
  geom_tile(color = "white") +
  geom_text(
    data = all_jaccard2 %>% distinct(rank, tool_ref, tool_comp, .keep_all = TRUE) %>% filter(jaccard >= 0.0),  # only these get labels
    aes(label = scales::percent(jaccard, accuracy = 1),
        color = jaccard >= 0.65), 
    size = general_size / .pt, show.legend = F) +
  scale_color_manual(values = c( "black", "white")) +
  facet_wrap(~ rank) +
  scale_fill_gradientn(
    colors = brewer.pal(9, "YlOrBr"),
    labels = percent_format(accuracy = 1),
    breaks = c(0, 0.2, 0.4, 0.6, 0.8),
    limits = c(0, NA)) +
  theme1 +
  theme(panel.grid  = element_blank()) +
  labs(x = NULL, y = NULL, fill = "Jaccard\nindex") 

all_csc_plot <- 
  all_jaccard %>% filter(tool_ref != tool_comp) %>%
  ggplot( aes(x = tool_comp, y = tool_ref, fill = csc)) +
  geom_tile(color = "white") +
  geom_text(
      data = all_jaccard %>% 
      filter(tool_ref != tool_comp) %>% 
      filter(csc >= 0.05),  # only these get labels
    aes(label = scales::percent(csc, accuracy = 1),
        color = csc >= 0.65), 
    size = general_size / .pt, show.legend = F) +
  scale_color_manual(values = c( "black", "white")) +
  facet_wrap(~ rank_label) +
  scale_fill_gradientn(
    colors = brewer.pal(9, "YlOrBr"),
    labels = percent_format(accuracy = 1),
    breaks = c(0, 0.2, 0.4, 0.6, 0.8),
    limits = c(0, NA)) +
  theme1 +
  theme(panel.grid  = element_blank()) +
  labs(x = NULL, y = NULL, fill = "Jaccard\nindex") 


all_jaccard %>% filter(tool_ref != tool_comp) %>%
  filter(rank %in% c("Reported ARGs", "Risk I", "Subset gene classes",
                     "Reference database (ID99%, COV90%)")) %>% 
  ggplot(aes(x = tool_ref, y = csc, fill = rank)) +
    geom_boxplot()

all_jaccard %>% filter(tool_ref != tool_comp) %>%
  filter(rank %in% c("Reported ARGs", "Risk I", "Subset gene classes",
                     "Reference database (ID99%, COV90%)")) %>% 
  ggplot(aes(x = tool_comp, y = csc, fill = rank)) +
  geom_boxplot()

all_jaccard %>% 
  filter(tool_ref != tool_comp) %>% 
  select(tool_ref, tool_comp, csc, rank) %>% 
  rename(x = csc) %>% 
  left_join(all_jaccard %>% filter(tool_ref != tool_comp) %>% 
              select(tool_comp, csc, rank, tool_ref) %>% 
              rename(y = csc, tool_ref = tool_comp, tool_comp = tool_ref), by = c("tool_ref","tool_comp", "rank")) %>% 
  filter(rank %in% c("Reported ARGs")) %>% #, "Risk I", "Subset gene classes",
                     #"Reference database (ID99%, COV90%)")) %>% 
  ggplot(aes(x = y, y = x, color = tool_ref)) +
  geom_point() + geom_line()


all_jaccard %>% 
  filter(tool_ref != tool_comp) %>% 
  select(tool_ref, tool_comp, csc, rank) %>% 
  mutate(direction = "by") %>% 
  bind_rows(all_jaccard %>% filter(tool_ref != tool_comp) %>% 
              select(tool_comp, csc, rank, tool_ref) %>% 
              mutate(direction = "on") %>% 
              rename(tool_ref = tool_comp, tool_comp = tool_ref)) %>% 
  filter(rank %in% c("Reported ARGs", "Risk I", "Subset gene classes",
                     "Reference database (ID99%, COV90%)")) %>%
  ggplot(aes(x = tool_ref, y = csc, color = rank, fill = rank)) +
  geom_boxplot(whisker.linetype = NA, outliers = F) +
  facet_grid(. ~ direction) + theme1 +
  ylab("Inclusion ratio") 

all_jaccard %>% 
  filter(tool_ref != tool_comp) %>% 
  select(tool_ref, tool_comp, csc, rank) %>% 
  mutate(direction = "by") %>% 
  bind_rows(all_jaccard %>% filter(tool_ref != tool_comp) %>% 
              select(tool_comp, csc, rank, tool_ref) %>% 
              mutate(direction = "on") %>% 
              rename(tool_ref = tool_comp, tool_comp = tool_ref)) %>% 
  filter(rank %in% c("Reported ARGs", "Risk I", "Subset gene classes",
                     "Reference database (ID99%, COV90%)")) %>%
  ggplot(aes(x = tool_ref, y = csc, color = rank, fill = rank)) +
  geom_boxplot(whisker.linetype = NA, outliers = F) +
  facet_grid(. ~ direction) + theme1 +
  ylab("Inclusion ratio") +
  coord_flip()


delta_jaccard <- all_jaccard2 %>%
  group_by(tool_ref, tool_comp) %>%
  mutate(
    jaccard_reported = jaccard[match("Reported ARGs", rank)]
  ) %>%
  ungroup() 


bind_rows(delta_jaccard,
          delta_jaccard %>% 
            rename(tool_ref = tool_comp, tool_comp = tool_ref))



all_pairs_summary <- purrr::map_dfr(
  combn(as.character(basic_tools), 2, simplify = FALSE),
  function(p) {
    get_unigene_classification(unigenes, db_cluster = db_cluster,
                               tool_a = p[1], tool_b = p[2],
                               pipeline_a_label = p[1], 
                               pipeline_b_label = p[2],
                               both_label = "Both",
                               cluster_level = "cluster_99_cov90",
                               lst) %>%
      ungroup() %>%
      group_by(query) %>%
      slice_head(n = 1) %>%
      ungroup() %>%
      group_by(tool_a, tool_b, detected_by_both) %>%
      summarise(n = n(), .groups = "drop")
  }
)



all_pairs_summary_2 <- all_pairs_summary %>%  
  group_by(tool_a, tool_b, detected_by_both) %>% 
  summarise(n = sum(n)) %>% 
  mutate(N = sum(n)) %>% 
  mutate(p = n / N) 

  
non_overlap <- all_pairs_summary_2 %>% filter(detected_by_both != "Overlap\n(Jaccard Index)") %>% 
  ungroup() %>% 
  group_by(tool_a, tool_b) %>% 
  mutate(p = n/sum(n))

non_overlap <- non_overlap %>% mutate(g = detected_by_both)
  
non_overlap %>% 
  filter(grepl("ABR",tool_a) | grepl("Res",tool_a), grepl("ABR",tool_b) | grepl("Res",tool_b)) %>% 
  ungroup() %>% 
  group_by(g) %>% 
  summarise(mn = mean(p), sd = sd(p))

non_overlap %>% 
  filter((grepl("RGI",tool_a) & grepl("Dee",tool_b)) | (grepl("RGI",tool_b) & grepl("Dee",tool_a))) %>% 
  ungroup() %>% 
  group_by(g)

non_overlap %>% 
  filter((grepl("RGI",tool_a) & grepl("CAR",tool_b)) | (grepl("RGI",tool_b) & grepl("CAR",tool_a))) %>% 
  ungroup() %>% 
  group_by(g)

non_overlap %>%  
  filter(grepl("RGI",tool_a) | grepl("RGI",tool_b)) %>% 
  filter(!grepl("fARG",tool_a), !grepl("fARG",tool_b), !grepl("AMR",tool_a), !grepl("AMR",tool_b)) %>% 
  group_by(g) %>% summarise(mn = mean(p), sd = sd(p))

non_overlap %>%  
  filter(grepl("ResF",tool_a) & grepl("ResF",tool_b), grepl("ResF",tool_b) & grepl("ResF",tool_a)) %>% 
  group_by(g) %>% summarise(mn = mean(p), sd = sd(p))

non_overlap %>% 
  filter(grepl("fARG",tool_a) | (grepl("fARG",tool_b))) %>% 
  ungroup() %>% 
  group_by(g) %>% summarise(mn = mean(p), sd = sd(p))

non_overlap %>% 
  filter(grepl("AMR",tool_a) | (grepl("AMR",tool_b))) %>% 
  ungroup() %>% 
  group_by(g) %>% summarise(mn = mean(p), sd = sd(p))


decompose_difference_plot <- all_pairs_summary_2 %>%
  bind_rows(all_pairs_summary_2 %>% rename(tool_a = tool_b, tool_b = tool_a)) %>% 
  mutate(tool_a = factor(tools_labels[tool_a], levels = tools_labels_factor),
         tool_b = factor(tools_labels[tool_b], levels = tools_labels_factor)) %>%
  mutate(detected_by_both = fct_rev(detected_by_both)) %>% 
  filter(as.numeric(tool_a) < as.numeric(tool_b)) %>%
  ggplot(aes(y = tool_b, x = p, fill = detected_by_both)) +
  geom_col(width = 0.8, color = "black", linewidth = 0.2) +
  facet_grid(. ~ tool_a, switch = "both") +
  labs(fill = "") +
  xlab("") +
  ylab("") +
  scale_fill_manual(values = rev(c("#e6ab02","#666666","#7570b3")), 
                    guide = guide_legend(nrow = 1, reverse = TRUE)) +
  theme(
    axis.text.x  = element_blank(),
    axis.text.y  = element_text(size = general_size),
    axis.ticks = element_blank(),
    strip.text = element_text(size = general_size, margin = margin(0, 0, 0, 0)),
    strip.background.x = element_blank(),
    strip.text.y.left = element_text(angle = 0),
    strip.text.x.bottom = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "bottom",
    plot.margin = margin(0, 0, 0, 0, unit = "pt"),
    legend.text = element_text(size = general_size),
    title = element_text(size = general_size + 2, face = "bold"),
    panel.background = element_blank()
  )

diff_cols <- c("#bf812d",
               "#f5f5f5",
            "#35978f") 

#diff_cols <- c("#e6ab02","#666666","#7570b3")

lv <- c("Overlap\n(Jaccard Index)", "Pipeline-driven\n difference", "Database-driven\n difference")

decompose_difference_plot3 <- all_pairs_summary_2 %>% ungroup() %>%
  bind_rows(rename(ungroup(all_pairs_summary_2), tool_a = tool_b, tool_b = tool_a)) %>%
  mutate(across(c(tool_a, tool_b), ~ factor(tools_labels[.x], levels = tools_labels_factor)),
         detected_by_both = factor(detected_by_both, levels = lv[c(2, 1, 3)])) %>%
  filter(as.integer(tool_a) < as.integer(tool_b)) %>%
  group_by(tool_a, tool_b) %>% arrange(detected_by_both, .by_group = TRUE) %>%
  mutate(xmax = cumsum(p) - p[detected_by_both == lv[2]] - p[detected_by_both == lv[1]] / 2) %>%
  ungroup() %>%
  mutate(s = 0.45 / max(abs(c(xmax, xmax - p))),
         cx = as.numeric(tool_a) + (xmax - p / 2) * s,
         cy = as.numeric(tool_b),
         w = p * s) %>%
  ggplot(aes(x = cx, y = cy, width = w, fill = detected_by_both)) +
  geom_tile(height = 0.7, color = "black", linewidth = 0.2) +
  scale_x_continuous(breaks = 1:9, labels = tools_labels_factor[1:9], expand = expansion(add = 0.5)) +
  scale_y_continuous(breaks = 2:10, labels = tools_labels_factor[2:10], expand = expansion(add = 0.5)) +
  scale_fill_manual(values = diff_cols, breaks = lv[c(2, 1, 3)], guide = guide_legend(nrow = 1)) +
  labs(fill = "", x = NULL, y = NULL) +
  theme1 +
  theme(axis.text.x = element_text(size = general_size, angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = general_size),
        axis.ticks = element_blank(), legend.position = "bottom",
        legend.text = element_text(size = general_size),
        plot.margin = margin(0, 0, 0, 0, unit = "pt"), panel.background = element_blank()) +
  theme(panel.border = element_blank(), panel.grid.minor  = element_blank())

# decompose_difference_plot3 <- all_pairs_summary_2 %>% ungroup() %>%
#   bind_rows(rename(ungroup(all_pairs_summary_2), tool_a = tool_b, tool_b = tool_a)) %>%
#   mutate(across(c(tool_a, tool_b), ~ factor(tools_labels[.x], levels = tools_labels_factor)),
#          detected_by_both = factor(detected_by_both, levels = lv[c(2, 1, 3)])) %>%
#   filter(as.integer(tool_a) < as.integer(tool_b)) %>%
#   group_by(tool_a, tool_b) %>% arrange(detected_by_both, .by_group = TRUE) %>%
#   mutate(xmax = cumsum(p) - p[detected_by_both == lv[2]] - p[detected_by_both == lv[1]] / 2,
#          xmid = xmax - p / 2) %>%
#   ggplot(aes(y = tool_b, x = xmid, width = p, fill = detected_by_both)) +
#   geom_tile(height = 0.7, color = "black", linewidth = 0.2) +
#   geom_vline(xintercept = 0, linewidth = 0.2, color = "grey") +
#   facet_grid(. ~ tool_a, switch = "both") +
#   labs(fill = "") +
#   xlab("") +
#   ylab("") +
#   scale_fill_manual(values = diff_cols, breaks = lv,
#                     guide = guide_legend(nrow = 1)) +
#   theme(
#     axis.text.x  = element_blank(),
#     axis.text.y  = element_text(size = general_size),
#     axis.ticks = element_blank(),
#     strip.text = element_text(size = general_size, margin = margin(0, 0, 0, 0)),
#     strip.background.x = element_blank(),
#     strip.text.y.left = element_text(angle = 0),
#     strip.text.x.bottom = element_text(angle = 90, vjust = 0.5, hjust = 1),
#     legend.position = "bottom",
#     plot.margin = margin(0, 0, 0, 0, unit = "pt"),
#     legend.text = element_text(size = general_size),
#     title = element_text(size = general_size + 2, face = "bold"),
#     panel.background = element_blank(),
#     panel.spacing.x = unit(3, "pt"),
#     panel.spacing.y = unit(3, "pt")) +
#   #geom_vline(xintercept = 0, linewidth = 0.15, colour = "grey") +
#   geom_tile(height = 0.7, color = "black", linewidth = 0.2) +
#   geom_segment(data = function(d) d %>% distinct(tool_a, tool_b) %>% group_by(tool_b) %>%
#                  mutate(xend = ifelse(as.integer(tool_a) == max(as.integer(tool_a)), 0, Inf)),
#                aes(x = -Inf, xend = Inf, y = tool_b, yend = tool_b),
#                inherit.aes = FALSE, linewidth = 0.2, colour = "grey") +
#   #geom_segment(data = function(d) d %>% distinct(tool_a, tool_b) %>% group_by(tool_a) %>%
#   #               summarise(y = min(as.numeric(tool_b)), yend = max(as.numeric(tool_b)), .groups = "drop"),
#   #             aes(x = 0, xend = 0, y = y, yend = yend),
#   #             inherit.aes = FALSE, linewidth = 0.15, colour = "grey") #+
#   geom_tile(height = 0.7, color = "black", linewidth = 0.2) 
# 


decompose_difference_plot2 <- all_pairs_summary_2 %>%
  bind_rows(all_pairs_summary_2 %>% rename(tool_a = tool_b, tool_b = tool_a)) %>% 
  mutate(tool_a = factor(tools_labels[tool_a], levels = tools_labels_factor),
         tool_b = factor(tools_labels[tool_b], levels = tools_labels_factor)) %>%
  mutate(detected_by_both = fct_rev(detected_by_both)) %>% 
  filter(as.numeric(tool_a) < as.numeric(tool_b)) %>%
  ggplot(aes(y = tool_b, x = p, fill = detected_by_both)) +
  geom_col(width = 0.8, color = "black", linewidth = 0.2, position = position_dodge2()) +
  facet_grid(. ~ tool_a, switch = "both") +
  labs(fill = "") +
  xlab("") +
  ylab("") +
  scale_fill_manual(values = rev(c("#e6ab02","#666666","#7570b3")), 
                    guide = guide_legend(nrow = 1, reverse = TRUE)) +
  theme(
    axis.text.x  = element_blank(),
    axis.text.y  = element_text(size = general_size),
    axis.ticks = element_blank(),
    strip.text = element_text(size = general_size, margin = margin(0, 0, 0, 0)),
    strip.background.x = element_blank(),
    strip.text.y.left = element_text(angle = 0),
    strip.text.x.bottom = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "bottom",
    plot.margin = margin(0, 0, 0, 0, unit = "pt"),
    legend.text = element_text(size = general_size),
    title = element_text(size = general_size + 2, face = "bold"),
    panel.background = element_blank()
  )


d <- all_pairs_summary_2 %>%
  mutate(tool_a = tools_labels[as.character(tool_a)],
         tool_b = tools_labels[as.character(tool_b)]) %>% 
  mutate(tool_a = factor(tool_a, levels = tools_labels_factor[1:10]),
         tool_b = factor(tool_b, levels = tools_labels_factor[1:10])) %>% 
  ungroup() %>% 
  bind_rows(., rename(., tool_a = tool_b, tool_b = tool_a)) %>% 
  mutate(x = as.numeric(tool_a), y = as.numeric(tool_b),
         g = fct_recode(detected_by_both,
                        jac  = "Overlap\n(Jaccard Index)",
                        pipe = "Pipeline-driven\n difference",
                        db   = "Database-driven\n difference")) %>%
  filter(x < y)


w <- d %>% 
  select(tool_a, tool_b, g, p) %>% 
  pivot_wider(names_from = g, values_from = p) %>% 
  mutate(db_share = db / (db + pipe))


decompose_B <- ggplot() +
  geom_tile(data = w, aes(tool_a, tool_b, fill = jac), colour = "grey", linewidth = 0.8) +
  geom_text(data = filter(w, jac >= 0.05),
            aes(tool_a, tool_b, color = jac > 0.5, label = scales::percent(jac, accuracy = 1)),
            size = general_size / .pt) +
  scale_color_manual(values = c( "black", "white")) +
  scale_fill_gradientn( colors = brewer.pal(9, "YlOrBr"),
                       labels = scales::percent, name = "Overlap\n(Jaccard index)",
                       guide = guide_colourbar(order = 1)) +
  new_scale_fill() +
  geom_tile(data = w, aes(tool_b, tool_a, fill = db_share), colour = "grey", linewidth = 0.8) +
  #geom_text(data = w, aes(tool_b, tool_a, label = scales::percent(db_share, accuracy = 1)),
  #          size = general_size / .pt) +
  scale_fill_gradientn(colours = brewer.pal(9, "BrBG"), limits = c(0, 1),
                       breaks = c(0, 1),
                       labels = c("100%\nPipeline-\ndriven", "100%\nDatabase-\ndriven"),
                       name = "Difference",
                       guide = guide_colourbar(order = 2)) +
  scale_x_discrete(drop = FALSE) +
  scale_y_discrete(drop = FALSE) +
  coord_fixed() +
  labs(x = NULL, y = NULL) +
  guides(color = "none") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = general_size),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.title = element_text(size = general_size, vjust = 0.8, hjust = 0),
        legend.text = element_text(size = general_size),
        legend.position = "bottom",
        plot.margin = margin(0, 0, 0, 0))





###

levels_abundance_div <- abundance_class %>% 
  group_by(habitat, tool, gene) %>% 
  summarise(a = sum(abundance), d = sum(richness)) %>% 
  arrange(desc(a), desc(d)) %>% 
  ungroup()


# classes to show in fig2
top_abundance <- c("efflux pump", "van" , "class A beta-lactamase", 
                   "tet RPG",  "cell wall charge", 
                   "rpoB", "erm", "aph")

# classes to show in fig4, the other classes are shown in fig_csc_sup_class

top_cso <- c("van", "efflux pump",  "tet RPG", "class A beta-lactamase", 
             "class B beta-lactamase","class C beta-lactamase", "class D beta-lactamase",
             "aph", "erm", "aac")


# aggregate the abundance in human gut by class and tool, add the gene class "other"
# if the gene class is not present in top_abundance

df_abundance_class_human <- abundance_class %>% 
  filter(habitat %in% "human gut") %>% 
  mutate(total = abundance) %>%
  mutate(gene = ifelse(gene %in% top_abundance, gene, "other")) %>% 
  mutate(texture = ifelse(tool %in% tools_texture, "yes", "no")) %>%
  mutate(gene = factor(gene, levels = c(top_abundance, "other"))) %>%
  ungroup() %>% 
  group_by(sample, tool, gene) %>%
  summarise(total = sum(total), texture = texture[1]) %>%
  mutate(tools_labels = factor(tools_labels[tool], levels = tools_labels[!duplicated(tools_labels)]),
         texture = ifelse(tool %in% tools_texture, "yes", "no"),
         tools_db = factor(tools_db[tool], levels = tools_db[!duplicated(tools_db)]))




# plot the number of unigenes considered ARG by tool

p2a <- unigenes %>%
  filter(tool %in% basic_tools) %>% 
  ungroup() %>% 
  group_by(tools_labels, texture, tools_db) %>% 
  summarise(n = n_distinct(query)) %>%
  ggplot(aes (x = tools_labels, y = n, fill = tools_db, pattern = texture )) +
  geom_col_pattern(position = position_dodge2(preserve = "single", width = 0.8), 
                   width = 0.8, pattern_color = "black", pattern_fill = pattern_fill, 
                   pattern_size =  0.12, color = "black") +
  scale_pattern_manual(values = c('no' = 'none', 'yes' = 'stripe')) +
  facet_grid( . ~ tools_db, scales = "free_x", space = "free") +
  scale_fill_manual(values = pal_7) +
  xlab("") + ylab("Number of ARGs") + 
  scale_y_continuous(expand = c(0.01, 0.01), 
                     breaks = c(25000,50000,75000,100000,125000), 
                     labels = scales::comma) + 
  guides(fill = guide_legend(
    override.aes = list(
      pattern = rep("none", 7),
      fill  = pal_7)), pattern = "none") + 
  theme1 +
  theme(panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        plot.margin = margin(0, 0, 0, 0, unit = "pt"),
        strip.text.x = element_text(size = general_size, angle = 0, vjust = 0, hjust = 0.5))

# Create the ECDF data for the plot, merging ABRicate and ResFinder, 
# and removing AMRFinderPlus with HMM and fARGene

id_plot_data <- unigenes %>% 
  filter(tool %in% basic_tools) %>% 
  filter(!tool %in% "fARGene") %>% 
  ungroup() %>%
  mutate(tool_2 = ifelse(tool %in% c("ABRicate-ARGANNOT", "ABRicate-MEGARes", "ABRicate-CARD",
                                     "ABRicate-ResFinder", "ABRicate-NCBI", "ResFinder"), "ABRicate/\nResFinder", 
                         as.character(tools_labels))) %>% 
  filter(!((tool_2 %in% "AMRFinder-\nPlus") & (query %in% lst$amrfinder.norm.prot$query[lst$amrfinder.norm.prot$Method %in% "HMM"]))) %>% 
  mutate(tool_2 = ifelse(tool_2 %in% "AMRFinder-\nPlus", "AMRFinder-\nPlus (aa)",tool_2)) %>% 
  mutate(tool_2 = factor(tool_2, levels = c("DeepARG", "RGI",  "ABRicate/\nResFinder", "AMRFinder-\nPlus (aa)"))) %>% 
  mutate(data_type = ifelse(tool_2 %in% c("DeepARG", "RGI", "ABRicate/\nResFinder"), "nucleotide", "amino acid")) %>% 
  mutate(data_type = factor(data_type, levels = c( "nucleotide", "amino acid"))) %>% 
  ungroup() %>% 
  group_by(tool_2) %>% mutate(N_obs = n()) %>% 
  mutate(tool_3 = paste0(tool_2, "\n(n = ", scales::comma(N_obs), ")"))


id_plot_data %>% group_by(tool_2) %>% summarise(n = n())

# plot the ECDF

id_plot <- id_plot_data %>% 
  mutate(tool_3 = factor(tool_3, levels = c("DeepARG\n(n = 100,075)","RGI\n(n = 65,937)",
                                            "ABRicate/\nResFinder\n(n = 39,170)", "AMRFinder-\nPlus (aa)\n(n = 3,141)"))) %>% 
  mutate(data_number = factor(data_type, levels = c( "nucleotide", "amino acid"))) %>% 
  ggplot(aes(x = id , color = tool_3, fill = tool_3)) +
  stat_ecdf(geom = "step", linewidth = 1) + 
  scale_color_manual(values = pal_7[c(7,1,4,2)], 
                     guide = guide_legend(nrow = 1)) +
  ggtitle("") + 
  ylab("Empirical Cumulative \nDistribution Function") + 
  xlab("% Identity") +
  labs(color = "", fill = "") +
  scale_y_continuous(limits = c(0,1), expand = c(0.01, 0.01)) + 
  scale_x_continuous(limits = c(0,100), expand = c(0, 0.5)) + 
  theme_minimal() + 
  theme1 +
  theme(legend.position = "bottom", panel.border = element_blank(),
        plot.margin = margin(0, 0, 0, 0, unit = "pt"))


# prepare the data for the proportion of gene class in each tool

df_plot <- unigenes %>% 
  filter(tool %in% basic_tools) %>% 
  group_by(tool, tools_labels, tools_db, gene_class) %>% 
  summarise(n = n()) %>% 
  mutate(p = n / sum(n)) %>% 
  ungroup() %>% 
  group_by(gene_class) %>%
  filter(max(p, na.rm = TRUE) >= 0.03) %>%  
  ungroup() %>%
  # change the name of gene classes to fit in the plot
  mutate(gene_class = gsub(" beta-lactamase","", gene_class)) %>%
  mutate(gene_class = gsub("rifampin inactivation enzyme","RIF-inact.", gene_class)) %>%
  mutate(gene_class = gsub("MFS efflux pump","MFS efflux", gene_class)) %>%
  mutate(gene_class = gsub("efflux pump","efflux", gene_class)) %>%
  mutate(gene_class = gsub("beta-lactam modulation resistance","beta-lactam\nmod.", gene_class)) %>%
  mutate(gene_class = gsub("target-modifying enzyme","target-modif.\nenzyme", gene_class)) %>%
  #mutate(gene_class = gsub("RIF-inact.", "RIF-inact.\nenz.",gene_class)) %>%
  mutate(gene_class = gsub("self-resistance", "self-\nresistance",gene_class)) %>%
  mutate(gene_class = gsub("cell wall charge","cell wall\ncharge", gene_class)) %>%
  mutate(gene_class = gsub("self-resistance","self-resistance", gene_class)) 

# plot the heatmap of the proportion of each gene class by tool
df_plot1 <-  df_plot %>% 
  mutate(gene_class = fct_reorder(gene_class, p, .fun = sum, .desc = FALSE)) %>%
  ungroup() %>% 
  filter(tool %in% basic_tools) %>% 
  ggplot(aes(x = tools_labels, y = gene_class, fill = p)) + 
  geom_tile(color = "grey") + 
  scale_fill_gradientn( colors = brewer.pal(9, "YlOrBr"),
                        labels = percent_format(accuracy = 1),
                        breaks = c(0.0001, 0.2, 0.4, 0.6)) + 
  facet_grid(.~tools_db, scales = "free_x", space = "free") +
  scale_y_discrete(labels = function(x) {
    out <- ifelse(
      x %in% c("rpoB", "van", "fos", "erm", "cat", "aph", "ant", "aac", "lnu", "nim","vat","mph","qnr"),
      paste0("italic('", x, "')"),
      ifelse(x %in% "abcF", "ABC-F",
             paste0("'", x, "'"))
    )
    parse(text = out)
  }) + 
  geom_text(
    data = df_plot %>% filter(tool %in% basic_tools, p >= 0.03),  # only these get labels
    aes(label = scales::percent(p, accuracy = 1),
        color = p >= 0.40), 
    size = general_size / .pt, show.legend = F) +
  scale_color_manual(values = c( "black", "white")) +
  labs(fill = "") + 
  ylab("Proportion of ARG class") + 
  xlab("") + 
  theme1 +
  theme(legend.position = "bottom", panel.grid = element_blank(),
        panel.border =  element_blank(),
        axis.text.y = element_text(size = general_size, angle = 0, vjust = 0.5, hjust = 1),
        axis.text.x = element_text(size = general_size, angle = 90, vjust = 0.5, hjust = 1),
        plot.margin = margin(0, 0, 0, 0, unit = "pt"),
        strip.text.x = element_text(size = general_size, angle = 0, vjust = 0.5, hjust = 0.5))




df_plot1_alt <-  df_plot %>% 
  mutate(gene_class = fct_reorder(gene_class, p, .fun = sum, .desc = FALSE)) %>%
  mutate(gene_class = fct_rev(fct_reorder(gene_class, p, .fun = sum))) %>%
  mutate(tools_labels = fct_rev(tools_labels)) %>%
  mutate(tools_db = fct_rev(tools_db)) %>%
  ungroup() %>% 
  filter(tool %in% basic_tools) %>% 
  ggplot(aes(x = gene_class, y = tools_labels, fill = p)) + 
  geom_tile(color = "grey") + 
  scale_fill_gradientn( colors = brewer.pal(9, "YlOrBr"),
                        labels = percent_format(accuracy = 1),
                        breaks = c(0.0001, 0.2, 0.4, 0.6)) + 
  facet_grid(tools_db ~ ., scales = "free_y", space = "free") +
  scale_x_discrete(labels = function(x) {
    out <- ifelse(
      x %in% c("rpoB", "van", "fos", "erm", "cat", "aph", "ant", "aac", "lnu", "nim","vat","mph","qnr"),
      paste0("italic('", x, "')"),
      ifelse(x %in% "abcF", "ABC-F",
             paste0("'", x, "'"))
    )
    parse(text = out)
  }) + 
  geom_text(
    data = df_plot %>% filter(tool %in% basic_tools, p >= 0.03),  # only these get labels
    aes(label = scales::percent(p, accuracy = 1),
        color = p >= 0.40), 
    size = general_size / .pt, show.legend = F) +
  scale_color_manual(values = c( "black", "white")) +
  labs(fill = "") + 
  xlab("Proportion of ARG class") + 
  ylab("") + 
  theme1 +
  theme(legend.position = "bottom", panel.grid = element_blank(),
        panel.border =  element_blank(),
        axis.text.x = element_text(size = general_size, angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = general_size, angle = 0, vjust = 0.5, hjust = 0),
        plot.margin = margin(0, 0, 0, 0, unit = "pt"),
        strip.text.y = element_text(size = general_size, angle = 90, vjust = 0.5, hjust = 0.5))


#top_abundance_JI


#JI_all_plot <- JI_all %>% 
#  filter(tool_ref %in% basic_tools, tool_comp %in% basic_tools) %>% 
#  mutate(tool_lab_ref = factor(tools_labels[tool_ref], levels = tools_labels[ordered_tools]),
#         tool_lab_comp = factor(tools_labels[tool_comp], levels = tools_labels[ordered_tools]))


# plot the jaccard index 

tri0 <- delta_jaccard %>% 
  filter(rank %in% c("Reference database (ID99%, COV90%)", "Reported ARGs")) %>% 
  select(tool_ref, tool_comp, rank, jaccard)

tri <- bind_rows(tri0,
                 tri0 %>% rename(tool_ref = tool_comp,
                                tool_comp = tool_ref)) %>% 
  distinct(tool_ref, tool_comp, rank, .keep_all = TRUE) %>% 
  mutate(tool_comp = factor(tools_labels[tool_comp], levels = tools_labels_factor[1:10])) %>% 
  mutate(tool_ref = factor(tools_labels[tool_ref], levels = tools_labels_factor[1:10])) %>% 
  mutate(
         x = as.numeric(tool_ref),
         y = as.numeric(tool_comp)) %>% 
  filter((rank == "Reference database (ID99%, COV90%)" & y > x) |
           (rank == "Reported ARGs"     & y < x))

tri2 <- bind_rows(tri0,
                 tri0 %>% rename(tool_ref = tool_comp,
                                tool_comp = tool_ref)) %>% 
  distinct(tool_ref, tool_comp, rank, .keep_all = TRUE) %>% 
  mutate(tool_comp = factor(tools_labels[tool_comp], levels = tools_labels_factor[1:10])) %>% 
  mutate(tool_ref = factor(tools_labels[tool_ref], levels = tools_labels_factor[1:10])) %>% 
  mutate(
    x = as.numeric(tool_ref),
    y = as.numeric(tool_comp)) %>% 
  filter((rank == "Reported ARGs" & y > x) |
           (rank ==   "Reference database (ID99%, COV90%)"   & y < x))



jaccard_plot1 <-   tri2 %>% 
  filter(rank == "Reported ARGs") %>% 
  ggplot(aes(x = tool_ref, y = tool_comp, fill = jaccard)) + 
  geom_tile(color = "grey") +
  scale_fill_gradientn( colors = brewer.pal(9, "YlOrBr"),
                        labels = percent_format(accuracy = 1),
                        breaks = c(0.25, 0.5, 0.75)) + 
  geom_text(
    data = tri2 %>% filter(rank == "Reported ARGs") %>% 
      filter(jaccard >= 0.1),
    aes(label = scales::percent(jaccard, accuracy = 1),
        color = jaccard >= .5), 
    size = general_size / .pt, show.legend = F) +
  scale_color_manual(values = c( "black", "white")) +
  labs(fill = "") + 
  ylab("") + 
  xlab("") + 
  theme1 +
  theme(legend.position = "bottom", panel.grid = element_blank(),
        panel.border =  element_blank(),
        plot.margin = margin(0, 0, 0, 0, unit = "pt"),
        strip.text.x = element_text(size = general_size, angle = 0, vjust = 0, hjust = 0.5))

jaccard_plot1

delta_jaccard2 <- delta_jaccard %>% bind_rows(
  delta_jaccard %>% rename(tool_comp = tool_ref, tool_ref = tool_comp)
) %>% 
  ungroup() %>% 
  group_by(tool_ref, rank) %>% 
  summarise(md = median(jaccard), q1 = quantile(jaccard, 0.25), q3 = quantile(jaccard, 0.75))


 
jc2 <- ggplot(delta_jaccard2 %>% 
                ungroup() %>% 
                filter(rank %in% c("Reference database (ID99%, COV90%)", "Reported ARGs")), 
              aes(x = tool_ref, y = md, color = rank, fill = rank, group = rank)) +
  geom_ribbon(aes(ymin = q1, ymax = q3), colour = NA) +
  geom_line(linewidth = 1) +
  geom_point(size = 1, color = "black") +
  scale_color_manual(values = c("#bf812d", "#35978f")) + 
  scale_fill_manual(values = c("#f6e8c3", "#c7eae5")) +
  scale_x_discrete(expand = c(0, 0.5)) +
  labs(x = "", y = "Jaccard index", color = "", fill = "") + 
  theme1 +
  theme(panel.border = element_blank(),
        legend.position = "bottom")

#jc2 <- jc2 + coord_flip()


jaccard_plot <-   tri %>% 
  ggplot(aes(x = tool_ref, y = tool_comp, fill = jaccard)) + 
  geom_tile(color = "grey") +
  scale_fill_gradientn( colors = brewer.pal(9, "YlOrBr"),
                        labels = percent_format(accuracy = 1),
                        breaks = c(0.25, 0.5, 0.75)) + 
  geom_text(
    data = tri %>% 
      filter((rank == "Reported ARGs" & jaccard >= 0.05) | 
               (rank != "Reported ARGs" & jaccard >= 0.3)),
    aes(label = scales::percent(jaccard, accuracy = 1),
        color = jaccard > 0.65), 
    size = general_size / .pt, show.legend = F) +
  scale_color_manual(values = c( "black", "white")) +
  labs(fill = "") + 
  ylab("") + 
  xlab("") + 
  theme1 +
  theme(legend.position = "bottom", panel.grid = element_blank(),
        panel.border =  element_blank(),
        panel.background = element_rect(fill = "white", colour = NA),
        plot.margin = margin(0, 0, 0, 0, unit = "pt"),
        strip.text.x = element_text(size = general_size, angle = 0, vjust = 0, hjust = 0.5)) +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) 

jaccard_plot



p2_1 <- ((p2a + ggtitle("a")) /
           (jaccard_plot1 + ggtitle("c")) /
           (id_plot + ggtitle("d")) +
           patchwork::plot_layout(heights = c(1,2.2,1)) |
           (df_plot1 + theme(plot.margin = margin(0, 0, 0, 20, unit = "pt")) + ggtitle("b"))) + 
  patchwork::plot_layout(widths = c(1,1))


top_p2 <- (p2a + ggtitle("a") | id_plot + ggtitle("b")+ theme(plot.margin = margin(0, 0, 0, 10, unit = "pt")))

left_bottom_p2 <- (decompose_difference_plot3 + ggtitle("c")) 

bottom_p2 <- (decompose_difference_plot3 + ggtitle("c")) | 
  (df_plot1 + 
     theme(plot.margin = margin(0, 0, 0, 10, unit = "pt"),
           strip.text = element_text(size = general_size, margin = margin(0, 0, 0, 0))) + 
  ggtitle("d"))


bottom_p2_2 <- (jaccard_plot + ggtitle("c")) | 
  (df_plot1 + 
   theme(plot.margin = margin(0, 0, 0, 10, unit = "pt"),
         strip.text = element_text(size = general_size, margin = margin(0, 0, 0, 0))) + 
   ggtitle("d"))

left2 <- ((p2a + ggtitle("a")) / 
  (id_plot + ggtitle("c")) /
  (decompose_difference_plot3 + ggtitle("d"))) + 
  patchwork::plot_layout(heights = c(1,1,3))
  
  
p2_1 <- (left2 | 
           (df_plot1 +ylab("")+ ggtitle("b") + theme(plot.margin = margin(0, 0, 0, 20, unit = "pt")))) +
  patchwork::plot_layout(widths = c(1,1)) 

p2_1_1 <- (top_p2 / bottom_p2_2) + 
  patchwork::plot_layout(heights = c(1,6))


plot_tools_per_unigene <- unigenes %>%
  filter(tool %in% basic_tools) %>%
  group_by(query) %>%
  summarise(n = n_distinct(tool)) %>%
  count(n, name = "unigenes") %>%
  mutate(pct = unigenes / sum(unigenes),
         n = factor(n, levels = 1:10)) %>%
  ggplot(aes(x = n, y = unigenes)) +
  geom_col(color = "#666666", width = 0.7) +
  geom_text(aes(label = paste0(scales::percent(pct, accuracy = 0.1))),
            vjust = -0.3, size = general_size / .pt, lineheight = 0.9) +
  scale_x_discrete(drop = FALSE) +
  labs(x = "Number of pipelines\nreporting an ARG", y = "ARGs") +
  scale_y_continuous(labels = scales::comma,
                     expand = expansion(mult = c(0, 0.15))) +   # more room above the bars
  coord_cartesian(clip = "off") + 
  theme1 +
  theme(axis.text.x = element_text(size = general_size, angle = 0),
        panel.border = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank())

p2_alt <- (p2a + ggtitle("a")| df_plot1_alt + ggtitle("b") + 
             theme(legend.position = "bottom",
                   legend.text  = element_text(size = general_size),
                   legend.title = element_text(size = general_size)) + 
             labs(x = "", fill = "")) + 
  patchwork::plot_layout(widths = c(1,3))

p2_alt2 <- ((p2a + ggtitle("a")) / 
              (id_plot + ggtitle("e") + scale_color_manual(values = pal_7[c(7,1,4,2)], 
                                                        guide = guide_legend(nrow = 2))))
              
p2_alt3 <- (df_plot1_alt + ggtitle("b") + 
  theme(legend.position = "bottom",
        legend.text  = element_text(size = general_size),
        legend.title = element_text(size = general_size)) + 
  labs(x = "", fill = "")) / 
  (decompose_difference_plot3 + ggtitle("d") + theme(plot.margin = margin(0, 0, 0, 5, unit = "pt")))

p3_alt <- ((id_plot + ggtitle("c") + 
              scale_color_manual(values = pal_7[c(7,1,4,2)], 
                                 guide = guide_legend(nrow = 2))) |
             (decompose_difference_plot3 + ggtitle("d") + theme(plot.margin = margin(0, 0, 0, 5, unit = "pt")))) +
  patchwork::plot_layout(widths = c(1,3))


p2_flip <- (p2a + ggtitle("a") + theme(strip.clip = "off",
            strip.text = element_text(hjust = 0))) +
  (df_plot1_alt + ggtitle("b") +
     theme(legend.position = "bottom",
           legend.text  = element_text(size = general_size),
           legend.title = element_text(size = general_size)) +
     labs(x = "", fill = "")) +
  (id_plot + ggtitle("c") +
     scale_color_manual(values = pal_7[c(7,1,4,2)],
                        guide = guide_legend(nrow = 2))) +
  (decompose_difference_plot3 + ggtitle("d") +
     theme(plot.margin = margin(0, 0, 0, 5, unit = "pt"))) +
  patchwork::plot_layout(ncol = 2, widths = c(1.4, 3), heights = c(1.3, 1))

ggsave("code_R_analysis/output_plots/fig2_flip.svg", p2_flip, 
       width = 180, height = 180, unit = "mm")

ggsave("code_R_analysis/output_plots/fig2.svg", p2_1, width = 180, height = 180, unit = "mm")
ggsave("code_R_analysis/output_plots/fig2_1.svg", p2_1_1, width = 180, height = 180, unit = "mm")
ggsave("code_R_analysis/output_plots/sup_jaccards.svg", all_jaccard_plot, width = 180, height = 180, unit = "mm")
ggsave("code_R_analysis/output_plots/sup_inclusion_ratio.svg", all_csc_plot, width = 180, height = 180, unit = "mm")


# tools 

# pre-calculate the information in the box plots abundance

lims_abundance <- abundance_tool_sample %>%
  group_by(tool, habitat,tools_db, tools_labels, texture) %>% 
  summarise(
    median = ifelse(quantile(abundance, 0.5) < 0, 0, quantile(abundance, 0.5)),
    q25 = ifelse(quantile(abundance, 0.25) < 0, 0, quantile(abundance, 0.25)),
    q75 = ifelse(quantile(abundance, 0.75) < 0, 0, quantile(abundance, 0.75)),
    w1 = ifelse(quantile(abundance, 0.25) - 1.5*IQR(abundance)<0,0,
                quantile(abundance, 0.25) - 1.5*IQR(abundance)),
    w2 = ifelse(quantile(abundance, 0.75) + 1.5*IQR(abundance)<0,0,
                quantile(abundance, 0.75) + 1.5*IQR(abundance))) 

# pre-calculate the information in the box plots richness

lims_richness <- abundance_tool_sample %>%
  group_by(tool, habitat,tools_db, tools_labels, texture) %>% 
  summarise(
    median = ifelse(quantile(richness, 0.5) < 0, 0, quantile(richness, 0.5)),
    q25 = ifelse(quantile(richness, 0.25) < 0, 0, quantile(richness, 0.25)),
    q75 = ifelse(quantile(richness, 0.75) < 0, 0, quantile(richness, 0.75)),
    w1 = ifelse(quantile(richness, 0.25) - 1.5*IQR(richness)<0,0,
                quantile(richness, 0.25) - 1.5*IQR(richness)),
    w2 = ifelse(quantile(richness, 0.75) + 1.5*IQR(richness)<0,0,
                quantile(richness, 0.75) + 1.5*IQR(richness))) 


# Create a list to where to keep the total abundance plot of each habitat
abu_plots <- list()

# set seed for selecting random observations to plot as dots within each tool and habitat 
set.seed(2026)

for(j in 1:length(EN)){
  
  # calculate total number of observations per habitat, and per tools and habitat
  df_plot <- abundance_tool_sample %>%
    filter(tool %in% basic_tools) %>% 
    filter(habitat %in% EN[j]) %>% 
    group_by(habitat) %>% 
    mutate(N = n_distinct(sample),
           abundance = abundance)
  
  # select the observations to plot as dots, keeping them below 1.5IQR for visualization
  df_jitter <- df_plot %>% 
    ungroup() %>%
    group_by(habitat, tool) %>% 
    filter(abundance < quantile(abundance, 0.75) + 1.5*IQR(abundance))
  
  # filter the limits for the tools and habitat 
  
  lims_abundance_max = lims_abundance %>%
    filter(habitat %in% EN[j]) %>%
    filter(tool %in% basic_tools) %>% 
    ungroup() %>% summarise(max(w2)) %>% pull()
  
  # plot 
  abu_plots[[j]] <- lims_abundance %>%
    filter(habitat %in% EN[j]) %>%
    filter(tool %in% basic_tools) %>% 
    mutate(habitat_label = paste0(habitat, "\n(n = ", scales::comma(df_plot$N[1]), ")")) %>%
    ggplot(aes(x = tools_labels, y=median, fill = tools_db))  + 
    geom_boxplot_pattern(aes(ymin = w1, lower = q25, middle = median, upper = q75, ymax = w2, pattern = texture), stat = "identity",
                         position = position_dodge2(preserve = "single", width = 0.3, padding = 0), 
                         width = .8, pattern_color = "black", pattern_fill = "black", pattern_density = 0.00001,
                         pattern_spacing = 0.15,
                         pattern_size =  0.4, color = "black", outliers = FALSE, outlier.shape = NA, linewidth = 0.15) +
    scale_x_discrete(expand = expansion(add = 1)) +
    geom_jitter(data = df_jitter %>%
                  slice_sample(n = min(100, df_jitter$N[1])), 
                aes(y = abundance),
                width = 0.35, size = 0.4, alpha = 0.1) + 
    facet_grid(habitat_label ~ tools_db  , scales = "free_x", space = "free_x") +
    scale_y_continuous(labels = scales::comma, limits = c(0, lims_abundance_max + 2)) + 
    scale_fill_manual(values = pal_7) +
    scale_pattern_manual(values = c('no' = 'none', 'yes' = 'stripe')) +
    theme1 + 
    theme(panel.border = element_blank(),
          strip.text.x = element_text(size = general_size, angle = 90, vjust = 0.5, hjust = 0),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          plot.margin = margin(1, 1, 1, 1, unit = "pt")) 
}

# plot figure 2 panel a in the paper

main_abundance_left <- 
  ((abu_plots[[1]] + theme(axis.text.x = element_blank(), strip.text.x = element_blank()) + xlab("") + ylab("") + ggtitle("a")) / 
     (abu_plots[[10]] + theme(axis.text.x = element_blank(), strip.text.x = element_blank()) + xlab("") +  ylab("Relative abundance\n(aligned reads per million)")) /
     (abu_plots[[11]] + theme(strip.text.x = element_blank()) + xlab("") + ylab("")))  + patchwork::plot_layout(heights = c(1,1,1))


main_abundance_right <- 
  ((abu_plots[[9]] + theme(axis.text.x = element_blank(), strip.text.x = element_blank()) + xlab("") + ylab("")) / 
     (abu_plots[[12]] + theme(axis.text.x = element_blank(), strip.text.x = element_blank()) + xlab("") +  ylab("")) /
     (abu_plots[[13]] + theme(strip.text.x = element_blank()) + xlab("") + ylab(""))) + patchwork::plot_layout(heights = c(1,1,1)) 

# this is actually panel a

left1 <- main_abundance_left | main_abundance_right

# get the gene levels for human gut abundance to use as factors 
# change names for plotting

gene_levels <- levels(df_abundance_class_human$gene)
gene_levels <- gsub(" beta-lactamase","", gene_levels)
gene_levels <- gsub("cell wall ","cell\nwall\n", gene_levels)
gene_levels <- gsub("MFS efflux pump","MFS\nefflux", gene_levels)
gene_levels <- gsub("efflux pump","efflux", gene_levels)


# calculate the class "other" comprising all other classes not in top_abundance

df_abundance_class_all <- abundance_class %>% 
  mutate(gene = ifelse(gene %in% top_abundance, gene, "other")) %>% 
  mutate(texture = ifelse(tool %in% tools_texture, "yes", "no")) %>%
  mutate(gene = factor(gene, levels = c(top_abundance, "other"))) %>%
  ungroup() %>% 
  group_by(habitat, sample, tool, gene) %>%
  summarise(abundance = sum(abundance), richness = sum(richness), texture = texture[1]) %>%
  mutate(tools_labels = factor(tools_labels[tool], levels = tools_labels[!duplicated(tools_labels)]),
         texture = ifelse(tool %in% tools_texture, "yes", "no"),
         tools_db = factor(tools_db[tool], levels = tools_db[!duplicated(tools_db)]))

# modify the name of the gene classes for plotting

df_a0 <- df_abundance_class_all  %>% 
  mutate(gene = gsub(" beta-lactamase","", gene)) %>%
  mutate(gene = gsub("cell wall ","cell\nwall\n", gene)) %>%
  mutate(gene = gsub("MFS efflux pump","MFS\nefflux", gene)) %>%
  mutate(gene = gsub("efflux pump","efflux", gene)) %>%
  mutate(gene = factor(gene, levels = gene_levels)) 

# create a new palette using as guide tool_levels 

pal_10_q_2 <- pal_10_q
names(pal_10_q_2) <- tools_levels


# plot panel b in figure 2

a00 <- df_a0 %>% 
  filter(tool %in% basic_tools) %>% 
  filter(habitat %in% c("human gut"), gene %in% "van") %>% 
  ungroup() %>% group_by(habitat) %>%  
  mutate(N_samples = length(unique(sample))) %>% 
  ungroup() %>%
  mutate(habitat_label = "") %>%
  ggplot(aes(x = abundance, y = fct_rev(tool), fill = tool, pattern = texture)) +
  geom_boxplot_pattern(coef = 0, position = position_dodge2(preserve = "single", width = 0.5, padding = 0), 
                       width = 1, pattern_color = "black", pattern_fill = pattern_fill, pattern_spacing = 0.03,
                       pattern_density = 0.12,
                       pattern_size =  0.05, color = "black", outliers = FALSE, outlier.shape = NA,
                       linewidth = 0.05) +
  scale_fill_manual(values = pal_10_q[1:10], labels = 
                      gsub("ABRicate-", "ABRicate-\n", 
                           gsub("AMRFinderPlus", "AMRFinder-\nPlus", 
                                gsub("RGI-DIAMOND", "RGI", basic_tools)))) +
  scale_pattern_manual(values = c('no' = 'none', 'yes' = 'stripe')) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
  facet_grid(gene ~ habitat_label, scales = "free", space = "free",
             labeller = labeller(
               gene = as_labeller(function(x) {
                 out <- ifelse(
                   x %in% c("rpoB", "van", "fos", "erm", "cat", 
                            "aph", "ant", "aac", "lnu", "nim", 
                            "vat", "mph", "qnr"),
                   paste0("italic('", x, "')"),
                   ifelse(x == "abcF", "ABC-F",
                          paste0("'", x, "'")))
                 return(out)}, default = label_parsed))) + 
  xlab("Relative abundance\n(aligned reads per million)") +
  scale_x_continuous(labels = scales::comma, expand = c(0,0)) +
  ylab("") + 
  labs(fill = "") + 
  ggtitle("") +
  theme1 + 
  theme(
    strip.text.x = element_text(size = general_size, vjust = 0, hjust = 0.5),    
    strip.text.y = element_text(size = general_size, vjust = 0, hjust = 0.5),
    legend.position = "bottom",
    panel.grid.major.y = element_blank(),
    axis.text.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.border =  element_blank(),
    plot.margin = margin(0, 0, 0, 0, unit = "pt")) +
  guides(fill = guide_legend(
    override.aes = list(
      pattern = c("none", "stripe","none","stripe","none","stripe",
                  rep("none", 4)),
      fill  = pal_10_q[1:10])), pattern = "none") +
  labs(title = NULL) +
  theme(strip.text.x = element_blank(),
        strip.background.x = element_blank(),
        plot.margin = margin(0, 0, 0, 0, unit = "pt"))

a01 <- df_a0 %>% 
  filter(tool %in% basic_tools) %>% 
  filter(habitat %in% c("human gut"), !gene %in% "van") %>% 
  ungroup() %>% group_by(habitat) %>%  
  mutate(N_samples = length(unique(sample))) %>% 
  ungroup() %>%
  mutate(habitat_label = paste0(habitat, " (n = ", scales::comma(N_samples[1]), ")")) %>%
  ggplot(aes(x = abundance, y = fct_rev(tool), fill = tool, pattern = texture)) +
  geom_boxplot_pattern(coef = 0, position = position_dodge2(preserve = "single", width = 0.5, padding = 0), 
                       width = 1, pattern_color = "black", pattern_fill = pattern_fill, pattern_spacing = 0.03,
                       pattern_density = 0.12,
                       pattern_size =  0.05, color = "black", outliers = FALSE, outlier.shape = NA,
                       linewidth = 0.05) +
  scale_fill_manual(values = pal_10_q[1:10], labels = 
                      gsub("ABRicate-", "ABRicate-\n", 
                           gsub("AMRFinderPlus", "AMRFinder-\nPlus", 
                                gsub("RGI-DIAMOND", "RGI", basic_tools)))) +
  scale_pattern_manual(values = c('no' = 'none', 'yes' = 'stripe')) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
  facet_grid(gene ~ habitat_label, scales = "free", space = "free",
             labeller = labeller(
               gene = as_labeller(function(x) {
                 out <- ifelse(
                   x %in% c("rpoB", "van", "fos", "erm", "cat", 
                            "aph", "ant", "aac", "lnu", "nim", 
                            "vat", "mph", "qnr"),
                   paste0("italic('", x, "')"),
                   ifelse(x == "abcF", "ABC-F",
                          paste0("'", x, "'")))
                 return(out)}, default = label_parsed))) + 
  xlab("Relative abundance\n(aligned reads per million)") +
  scale_x_continuous(labels = scales::comma, expand = c(0,0)) +
  ylab("") + 
  labs(fill = "") + 
  ggtitle("b") +
  theme1 + 
  theme(
    strip.text.x = element_text(size = general_size, vjust = 0, hjust = 0.5),    
    strip.text.y = element_text(size = general_size, vjust = 0, hjust = 0.5),
    legend.position = "bottom",
    panel.grid.major.y = element_blank(),
    axis.text.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.border =  element_blank(),
    plot.margin = margin(0, 0, 0, 0, unit = "pt")) +
  guides(fill = guide_legend(
    override.aes = list(
      pattern = c("none", "stripe","none","stripe","none","stripe",
                  rep("none", 4)),
      fill  = pal_10_q[1:10])), pattern = "none") +
  xlab(NULL) +
  theme(axis.text.x = element_text(angle = 0, vjust = 1, size = general_size),
        axis.ticks.length.x = unit(1, "pt"),
        plot.margin = margin(0, 0, 0, 0, unit = "pt"))

a0 <- ((a01 + theme(legend.position = "none") + xlab(""))/ 
         (a00 + theme(legend.position = "none"))) + 
  patchwork::plot_layout(heights = c(9,1))

#a0

# merge the figure 2 panel a with the appropriate legend


library(uwot)

hg <- abundance_class %>% filter(habitat == "human gut", !is.na(new_level), tool %in% basic_tools)

joint_umap <- function(d, value, seed = 1) {
  w <- d %>% pivot_wider(id_cols = c(sample, tools_labels, tool, tools_db, texture), names_from = gene,
                         values_from = {{ value }}, values_fn = sum, values_fill = 0)
  m <- log10(as.matrix(select(w, -sample, -tools_labels, -tool, -tools_db, -texture)) + 1)
  keep <- rowSums(m) > 0
  set.seed(seed)
  u <- umap(m[keep, apply(m[keep, ], 2, var) > 0], n_neighbors = 30, min_dist = 0.3)
  tibble(w[keep, c("sample", "tools_labels", "tool", "tools_db", "texture")], u1 = u[, 1], u2 = u[, 2])
}

umap1 <- joint_umap(hg, abundance)
umap2 <- joint_umap(hg, richness)

umap3 <- joint_umap(abundance_class, abundance)
umap3$habitat <- abundance_class$habitat[match(umap3$sample, abundance_class$sample)]
umap3 <- umap3 %>% 
  mutate(place = ifelse(habitat %in% c("soil", "marine", "freshwater", "soil","wastewater"), "External", 
  ifelse(grepl("gut", habitat), "Gut-associated", "Human (non-gut)"))) 

set.seed(2)
umap1.1 <- umap1 %>% filter(sample %in% sample(unique(sample), 500))
set.seed(2)
umap2.1 <- umap2 %>% filter(sample %in% sample(unique(sample), 500))

umap3 <- umap3 %>% 
  mutate(place = ifelse(habitat %in% c("soil", "marine", "freshwater", "soil","wastewater"), "External", 
                        ifelse(grepl("gut", habitat), "Gut-associated", "Human (non-gut)"))) 

set.seed(2030)
ids <- umap3 %>%
  distinct(sample, place, habitat) %>%
  group_by(habitat) %>%
  slice_sample(n = 200) %>% 
  ungroup() %>% 
  distinct(sample, place) %>%
  group_by(place) %>%
  slice_sample(n = 300) %>% 
  ungroup()

umap3.1 <- umap3 %>% semi_join(ids, by = c("sample", "place"))


draw_key_point_line <- function(data, params, size) {
  pt <- ggplot2::draw_key_point(data, params, size)
  if (isTRUE(as.numeric(data$shape) == 24)) {
    ln <- grid::segmentsGrob(0.15, 0.15, 0.85, 0.85,
                             gp = grid::gpar(col = "black", lwd = 1))
    grid::grobTree(pt, ln)
  } else pt
}


u1_plot <- umap1.1 %>%
  ggplot(aes(u1, u2, colour = tool, fill = tool, shape = texture)) +
  geom_point(size = 1.5, color = "black", stroke = .2, key_glyph = draw_key_point_line) +
  scale_colour_manual(values = pal_10_q, name = "") +
  scale_fill_manual(values = pal_10_q, name = "") +
  guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1), nrow = 2)) +
  theme1 + 
  xlab("UMAP1") +
  ylab("UMAP2") +
  scale_shape_manual(values = c("no" = 21, "yes" = 24)) + 
  guides(fill = 
         guide_legend(nrow = 2, 
         override.aes = list(shape = shape_tools,color = "black", stroke = 0.4)), 
         shape = "none")  


u2_plot <- umap2.1 %>%
  ggplot(aes(u1, u2, colour = tool, fill = tool, shape = texture)) +
  geom_point(size = 1.5, color = "black", stroke = .2, key_glyph = draw_key_point_line) +
  scale_colour_manual(values = pal_10_q, name = "") +
  scale_fill_manual(values = pal_10_q, name = "") +
  guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1), nrow = 2)) +
  theme1 + 
  xlab("UMAP1") +
  ylab("UMAP2") +
  scale_shape_manual(values = c("no" = 21, "yes" = 24)) + 
  guides(fill = 
           guide_legend(nrow = 2, override.aes = list(shape = shape_tools,
           color = "black", stroke = 0.4)), shape = "none")  


u3_plot1 <- umap3.1 %>%
  mutate(al = place == "External") %>% 
  ggplot(aes(u1, u2, colour = tools_labels, fill = tools_labels, shape = texture, size = al, alpha = al)) +
  geom_point(color = "black", stroke = .2, key_glyph = draw_key_point_line) +
  scale_size_manual(values = c("TRUE" = 0.7, "FALSE" = 0.05), guide = "none") +
  scale_colour_manual(values = pal_tool, name = "") +
  scale_fill_manual(values = pal_tool, name = "") +
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.1), guide = "none") +
  guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1), nrow = 2)) +
  theme1 + 
  theme(plot.subtitle = element_text(size = general_size, color = "black")) +
  xlab("UMAP1") +
  ylab("UMAP2") +
  labs(subtitle = "External habitats") + 
  scale_shape_manual(values = c("no" = 21, "yes" = 24)) + 
  guides(fill = 
           guide_legend(nrow = 2, override.aes = list(shape = shape_tools,
                                                      color = "black", stroke = 0.4)), shape = "none")  

u3_plot2 <- umap3.1 %>%
  mutate(al = place == "Gut-associated") %>% 
  ggplot(aes(u1, u2, colour = tool, fill = tool, shape = texture, size = al, alpha = al)) +
  geom_point(color = "black", stroke = .2, key_glyph = draw_key_point_line) +
  scale_size_manual(values = c("TRUE" = 0.7, "FALSE" = 0.05), guide = "none") +
  scale_colour_manual(values = pal_10_q, name = "") +
  scale_fill_manual(values = pal_10_q, name = "") +
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.1), guide = "none") +
  guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1), nrow = 2)) +
  theme1 + 
  theme(plot.subtitle = element_text(size = general_size, color = "black")) +
  xlab("UMAP1") +
  ylab("UMAP2") +
  labs(title = "c", subtitle = "Gut-associated") + 
  scale_shape_manual(values = c("no" = 21, "yes" = 24)) + 
  guides(fill = 
           guide_legend(nrow = 2, override.aes = list(shape = shape_tools,
                                                      color = "black", stroke = 0.4)), shape = "none")  

u3_plot3 <- umap3.1 %>%
  mutate(al = place == "Human (non-gut)") %>% 
  ggplot(aes(u1, u2, colour = tool, fill = tool, shape = texture, size = al, alpha = al)) +
  geom_point(color = "black", stroke = .2, key_glyph = draw_key_point_line) +
  scale_size_manual(values = c("TRUE" = 0.7, "FALSE" = 0.05), guide = "none") +
  scale_colour_manual(values = pal_10_q, name = "") +
  scale_fill_manual(values = pal_10_q, name = "") +
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.1), guide = "none") +
  guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1), nrow = 2)) +
  theme1 + 
  theme(plot.subtitle = element_text(size = general_size, color = "black")) +
  xlab("UMAP1") +
  ylab("UMAP2") +
  labs(subtitle = "Human (non-gut)") + 
  scale_shape_manual(values = c("no" = 21, "yes" = 24)) + 
  guides(fill = 
           guide_legend(nrow = 2, override.aes = list(shape = shape_tools,
                                                      color = "black", stroke = 0.4)), shape = "none")  

#u3_plot2 + xlab("") | u3_plot3 + ylab("") | u3_plot1 + ylab("") + xlab("")

main_abundance_left <- 
  ((abu_plots[[1]] + theme(axis.text.x = element_blank(), strip.text.x = element_blank()) + xlab("") + ylab("") + ggtitle("a")) / 
     (abu_plots[[10]] + theme(axis.text.x = element_blank(), strip.text.x = element_blank()) + xlab("") +  ylab("Relative abundance\n(aligned reads per million)")) /
     (abu_plots[[11]] + theme(axis.text.x = element_blank(),strip.text.x = element_blank()) + xlab("") + ylab("")))  + patchwork::plot_layout(heights = c(1,1,1))


main_abundance_right <- 
  ((abu_plots[[9]] + theme(axis.text.x = element_blank(), strip.text.x = element_blank()) + xlab("") + ylab("")) / 
     (abu_plots[[12]] + theme(axis.text.x = element_blank(), strip.text.x = element_blank()) + xlab("") +  ylab("")) /
     (abu_plots[[13]] + theme(axis.text.x = element_blank(), strip.text.x = element_blank()) + xlab("") + ylab(""))) + patchwork::plot_layout(heights = c(1,1,1)) 

# this is actually panel a

left1 <- main_abundance_left | main_abundance_right

leg <- g_legend(u3_plot1 + ylab("") + xlab("") + theme(legend.position = "bottom"))



#umap_together <- wrap_elements(full = leg) / 
#  (((u2_plot + ggtitle("c") + theme(plot.margin = margin(10, 0, 0, 0, unit = "pt")))) | 
#     (u1_plot + xlab("") + ylab("") + ggtitle("d") + theme(plot.margin = margin(10, 0, 0, 0, unit = "pt")))) +
#  patchwork::plot_layout(heights = c(1,15))

umap_together <- (wrap_elements(full = leg) + theme(plot.margin = margin(5, 0, 0, 0, unit = "pt"))) / 
  (u3_plot2 + xlab("") | u3_plot3 + ylab("") + theme(plot.margin = margin(5, 0, 0, 0, unit = "pt")) | 
  u3_plot1 + ylab("") + xlab("") + theme(plot.margin = margin(5, 0, 0, 0, unit = "pt"))) +
  patchwork::plot_layout(heights = c(1,6))

left_block <- (left1 /
                 umap_together) + 
  patchwork::plot_layout(heights = c(1,1))

# merge all panels in figure 2 

p30 <- (left_block | (a0 + theme(legend.position = "none"))) + patchwork::plot_layout(widths = c(4,1))
ggsave("code_R_analysis/output_plots/fig3.svg", p30, width = 180, height = 180, unit = "mm")


set.seed(2)
ids <- sample(unique(umap1$sample), 500)

w <- hg %>% filter(!is.na(new_level), sample %in% ids) %>%
  pivot_wider(id_cols = c(sample, tools_labels), names_from = new_level,
              values_from = abundance, values_fn = sum, values_fill = 0)

m <- log10(as.matrix(select(w, -sample, -tools_labels)) + 1)
d <- vegan::vegdist(m, method = "euclidean")     # about 5,000 rows, about 12 million distances

# pipeline effect, permutations restricted within sample
vegan::adonis2(d ~ tools_labels, data = w,
               permutations = permute::how(blocks = w$sample, nperm = 199))



sup_abundance_left <- 
  (abu_plots[[2]] + theme(axis.text.x = element_blank(), strip.text.x = element_blank()) + xlab("") + ylab("")) / 
  (abu_plots[[3]] + theme(axis.text.x = element_blank(), strip.text.x = element_blank()) + xlab("") + ylab("Relative abundance\n(aligned reads per million)")) /
  (abu_plots[[4]] + theme(axis.text.x = element_blank(), strip.text.x = element_blank()) + xlab("") + ylab("")) /
  (abu_plots[[5]] + theme(strip.text.x = element_blank()) + xlab("") + ylab(""))

sup_abundance_right <- 
  (abu_plots[[6]] + theme(axis.text.x = element_blank(), strip.text.x = element_blank()) + xlab("") + ylab("")) / 
  (abu_plots[[7]] + theme(axis.text.x = element_blank(), strip.text.x = element_blank()) + xlab("") + ylab("")) /
  (abu_plots[[8]] + theme(strip.text.x = element_blank()) + xlab("") + ylab("")) 

sup_abundance <- sup_abundance_left | sup_abundance_right

ggsave("code_R_analysis/output_plots/sup_abundance.svg", sup_abundance, width = 180, height = 180, unit = "mm")



###




################################################################################################

# create vector for shaping the geom_points in the core-resistome and pan-resistome plots


shape_tools <- rep(21, length(tools_labels[tools_levels %in% basic_tools]))
shape_tools[tools_levels[tools_levels %in% basic_tools] %in% tools_texture] <- 24

# for the pan and core, define only one metric name for both and add the type of data (core or pan)

pan_core <- bind_rows(pan %>% mutate(metric = "Pan-resistome") %>% 
                        group_by(metric, tool, habitat, texture) %>% 
                        summarise(size = sum(p_pan), .groups = "drop"),
                      
                      core %>% mutate(metric = "Core-resistome") %>% 
                        group_by(metric, tool, habitat, texture) %>% 
                        summarise(size = n_distinct(centroid), .groups = "drop"))

# plot the pan resistome

p4a1 <- pan_core %>% filter(metric %in% "Pan-resistome") %>% 
  filter(habitat %in% c("human gut", "pig gut", "wastewater", "marine", "freshwater", "soil")) %>% 
  ggplot(aes(y = fct_rev(tool), x =  size)) +
  geom_point(aes(fill = tool, shape = texture),  color = "black", stroke = 0.3, size = 1.7) + 
  facet_grid(habitat ~ metric, scales = "free") +
  scale_fill_manual(values = pal_10_q[1:10], labels = 
                      gsub("ABRicate-", "ABRicate-\n", 
                           gsub("AMRFinderPlus", "AMRFinder-\nPlus", 
                                gsub("RGI-DIAMOND", "RGI", basic_tools)))) +
  scale_shape_manual(values = c("no" = 21, "yes" = 24)) +
  theme_minimal() +
  ylab("") +
  xlab("Number of ARGs") +
  labs(fill = "") +
  ggtitle("a") + 
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
  theme1 +
  geom_segment(aes(x = 0, xend = size, y = tool, yend = tool, color = tool), linewidth = 0.2, show.legend = F) +
  scale_color_manual(values = pal_10_q) +
  guides(fill = guide_legend(
    override.aes = list(
      shape = shape_tools,
      color  = pal_10_q[1:10])), shape = "none") + 
  theme(axis.text.y = element_blank(), 
        strip.text.x = element_text(size = general_size, vjust = 0, hjust = 0.5),
        strip.text.y = element_text(size = general_size, vjust = 0.5, hjust = 0.5, angle = 0),
        legend.position = "bottom",
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.border = element_blank(),
        plot.margin = margin(5, 0, 0, 0, unit = "pt")) +
  scale_x_reverse(labels = label_comma()) 

# plot the core resistome

p4a2 <- pan_core %>% filter(!metric %in% "Pan-resistome") %>% 
  filter(habitat %in% c("human gut", "pig gut", "wastewater", "marine", "freshwater", "soil")) %>% 
  filter(tool %in% basic_tools) %>%
  ggplot(aes(y = fct_rev(tool), x =  size)) +
  geom_point(aes(fill = tool, shape = texture),  color = "black", stroke = 0.3, size = 1.7) + 
  facet_grid(habitat ~ metric, scales = "free") +
  scale_fill_manual(values = pal_10_q) +
  scale_shape_manual(values = c("no" = 21, "yes" = 24)) +
  geom_segment(aes(x = 0, xend = size, y = tool, yend = tool, color = tool), linewidth = 0.2, show.legend = F) +
  scale_color_manual(values = pal_10_q) +
  theme_minimal() +
  ylab("") +
  xlab("Number of ARGs") + 
  ggtitle("b") + 
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
  theme1 +
  theme(axis.text.y = element_blank(), 
        strip.text.x = element_text(size = general_size, vjust = 0.5, hjust = 0.5),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.text.y = element_blank(),
        panel.border = element_blank(),
        plot.margin = margin(5, 0, 0, 0, unit = "pt")) +
  scale_x_continuous(labels = label_comma()) 

# merge pan and core resistome and add the legend of the plot 

panel1 <- (p4a1 + theme(legend.position = "none") | p4a2) + patchwork::plot_layout(widths = c(1,1))
p4 <- (panel1 / (g_legend(p4a1) )) + patchwork::plot_layout(heights = c(12,1))

#ggsave("code_R_analysis/output_plots/fig4.svg", p4, width = 120, height = 140, unit = "mm")

# 

pan_core_wide <- bind_rows(pan %>% mutate(metric = "Pan") %>% 
                        group_by(metric, tool, habitat, texture, new_level) %>% 
                        summarise(size = sum(p_pan), .groups = "drop"),
                      
                      core %>% mutate(metric = "Core") %>% 
                        group_by(metric, tool, habitat, texture, new_level) %>% 
                        summarise(size = n_distinct(centroid), .groups = "drop")) %>% 
  pivot_wider( names_from  = metric, values_from = size) %>% 
  mutate(Pan = ifelse(is.na(Pan), 0, Pan),
         Core = ifelse(is.na(Core), 0, Core))
           
pan_core_wide %>% 
  filter(habitat %in% c("human gut", "pig gut", "wastewater", "marine", "freshwater", "soil")) %>% 
  ggplot(aes(x = Pan, y = Core,  shape = texture, fill = tool)) +
  geom_point() + 
  facet_grid(new_level ~ habitat, scales = "free") +
  scale_fill_manual(values = pal_10_q) +
  scale_shape_manual(values = c("no" = 21, "yes" = 24)) +
  theme_minimal() 
  

pc <- pan_core_wide %>% 
  filter(Pan >= 1) %>%                          # drop classes with < 1 expected centroid
  mutate(r = pmin(Core / Pan, 1))




pan_core_sum <- pan_core_wide %>% 
  group_by(tool, habitat, texture) %>% 
  summarise(Pan = sum(Pan), Core = sum(Core)) %>% 
  mutate(r = Core / Pan) %>% 
  group_by(place = ifelse(habitat %in% c("soil", "marine", "freshwater", "soil","wastewater"), 
                          "External", 
                          ifelse(grepl("gut", habitat), "Gut-associated", "Human (non-gut)"))) %>% 
  mutate(tools_db = factor(tools_db[tool], levels = tools_db_factor),
         tools_labels = factor(tools_labels[tool], levels = tools_labels_factor)) %>% 
  ungroup() %>% 
  group_by(place) %>% 
  mutate(hab_id = factor(dense_rank(as.character(habitat)))) %>%
  ungroup() 


lims_pan_core <- pan_core_sum %>% 
  ungroup() %>% 
  group_by(place, tools_labels, tool) %>% 
  summarise(minr = min(r), maxr = max(r),
            minpan = min(Pan), maxpan = max(Pan),
            mincore = min(Core), maxcore = max(Core))

p_ratio <- pan_core_wide %>% 
  group_by(tool, habitat, texture) %>% 
  summarise(Pan = sum(Pan), Core = sum(Core)) %>% 
  filter(Pan >= 1) %>%
  mutate(r = Core / Pan) %>% 
  group_by(place = ifelse(habitat %in% c("soil", "marine", "freshwater", "soil","wastewater"), 
                          "external", 
                          ifelse(grepl("gut", habitat), "gut-associated", "other human"))) %>%
  mutate(tools_db = factor(tools_db[tool], levels = tools_db_factor),
         tools_labels = factor(tools_labels[tool], levels = tools_labels_factor)) %>% 
  ggplot(aes(y = tools_labels, x = r, fill = tool, pattern = texture)) +
  geom_boxplot_pattern(position = position_dodge2(preserve = "single", width = 0.3, padding = 0), 
                       width = 0.7, pattern_color = "black", 
                       pattern_fill = "black", 
                       pattern_density = 0.000000001,
                       pattern_spacing = 0.1,
                       pattern_size =  0.3, 
                       color = "black", outliers = F) +
  facet_grid(place ~ ., scales = "free_y", space = "free_y") +
  scale_fill_manual(values = pal_10_q[1:10]) +
  scale_pattern_manual(values = c('no' = 'none', 'yes' = 'stripe')) +
  labs(fill = "") +
  theme1 +
  theme(legend.position = "bottom") + 
  guides(fill = guide_legend(
    override.aes = list(
      pattern = c("none", "stripe","none","stripe","none","stripe",
                  rep("none", 4)),
      fill  = pal_10_q[1:10])), pattern = "none")


hab_lab <- pan_core_sum %>% 
  distinct(place, hab_id, habitat) %>% 
  arrange(place) %>% 
  group_by(hab_id) %>% 
  summarise(lab = paste(habitat, collapse = "\n")) %>%   # or collapse = " / "
  { setNames(.$lab, .$hab_id) }

tool_lv   <- levels(droplevels(factor(pan_core_sum$tools_labels)))
hatch_key <- tool_lv %in% pan_core_sum$tools_labels[pan_core_sum$texture == "yes"]

draw_key_hatch <- function(data, params, size) {
  base <- draw_key_point(data, params, size)
  if (isTRUE(data$hatch)) {
    x <- data
    x$shape <- 4; x$colour <- "black"; x$fill <- NA; x$size <- 1.8; x$stroke <- 0.35
    grid::grobTree(base, draw_key_point(x, params, size))
  } else base
}


p_core_ratio <- list()

places <- unique(pan_core_sum$place)

for(j in 1:length(places)){
  p_core_ratio[[j]] <- 
    pan_core_sum %>% 
    filter(r > 0,
           place %in% places[j]) %>% 
    ggplot(aes(y = tools_labels, x = r, colour = tool, shape = habitat)) +
    geom_segment(data = pan_core_sum %>% 
                   filter(r > 0,
                          place %in% places[j]) %>% 
                   ungroup() %>% 
                   group_by(place, tools_labels, tool) %>% 
                   summarise(minr = min(r), maxr = max(r),
                             minpan = min(Pan), maxpan = max(Pan),
                             mincore = min(Core), maxcore = max(Core)) %>%
                   ungroup() , 
                 aes(x = minr, xend = maxr, y = tools_labels, yend = tools_labels, colour = tool), 
                 linewidth = 0.7, show.legend = FALSE, inherit.aes = FALSE) +
    geom_point(aes(fill = tool), stroke = 0.5, size = 2.5, key_glyph = draw_key_hatch) +
    geom_point(data = ~ filter(.x, texture == "yes"),
               shape = 4, colour = "black", size = 1.8, stroke = 0.35,
               show.legend = FALSE) +
    facet_grid(place ~ ., scales = "free_x", space = "free_x") +
    scale_colour_manual(values = pal_10_q[seq_along(tool_lv)]) +
    scale_fill_manual(values = pal_10_q[seq_along(tool_lv)], guide = "none") +
    guides(colour = guide_legend(override.aes = list(shape = 21,
                                                     fill  = pal_10_q[seq_along(tool_lv)],
                                                     hatch = hatch_key))) +
    labs(colour = "") +
    theme1 +
    theme(legend.position = "right") + 
    scale_x_continuous(labels = scales::comma, expand = c(0, 0.005)) +
    scale_shape_manual(values = c(21, 22, 23, 24, 25), labels = hab_lab,
                       guide = guide_legend(override.aes = list(fill = "grey60", colour = "black"),
                                            title = ""))   
}


(p_core_ratio[[3]]  + guides(colour = "none") + xlab("") + ylab(""))  /
  (p_core_ratio[[1]]  + guides(colour = "none")+ xlab("") + ylab(""))  /
  (p_core_ratio[[2]]  + guides(colour = "none") + xlab("Core/Pan")+ ylab(""))



  
plot_core_pan_list <- list()

tool_lv <- c("RGI", "ABRicate-CARD", "AMRFinderPlus", "ABRicate-NCBI", "ResFinder",
             "ABRicate-ResFinder", "ABRicate-ARGANNOT", "ABRicate-MEGARes", "fARGene", "DeepARG")  # your order
pal_tool <- setNames(pal_10_q[seq_along(tool_lv)], tool_lv)
names(pal_tool) <- gsub("^(ABRicate|AMRFinder)-", "\\1-\n", names(pal_tool))
names(pal_tool)[names(pal_tool) == "AMRFinderPlus"] <- "AMRFinder-\nPlus"

for(j in 1:length(places)){
  plot_core_pan_list[[j]] <- 
    pan_core_sum %>% 
    rename(tool_name = tools_labels) %>%
    filter(habitat %in% c("human gut", "pig gut", "dog gut", "human nose", "wastewater",
                          "marine", "soil", "human oral", "human skin")) %>% 
    filter(
           place %in% places[j]) %>% 
    ggplot(aes(y = log10(Core+1), x = log10(Pan+1), shape = habitat)) +
    geom_line(aes(group = habitat), colour = "#666666", linewidth = 0.3, show.legend = FALSE) +
    geom_point(aes(fill = tool_name), color = "black", stroke = 0.15, size = 2.5, key_glyph = draw_key_hatch) +
    geom_point(data = ~ filter(.x, texture == "yes"), shape = 4, colour = "black",
               size = 1.8, stroke = 0.35, show.legend = FALSE) +
    facet_grid(place ~ ., scales = "free_x", space = "free_x") +
    scale_fill_manual(values = pal_tool, guide = "none") +
    scale_shape_manual(values = c(21, 22, 23, 24, 25), labels = hab_lab,
                       guide = guide_legend(override.aes = list(fill = "#666666", colour = "black"),
                                            title = "")) +
    theme1 + 
    theme(legend.position = "inside",
          legend.position.inside = c(0.001, 1.1),
          legend.justification = c(0, 1),
          legend.background = element_blank(),
          legend.key = element_blank()) + 
    xlim(c(min(log10(pan_core_sum$Pan+1)), max(log10(pan_core_sum$Pan+1)))) +
    ylim(c(0, max(log10(pan_core_sum$Core+1)))) +
    xlab(expression(log[10](Pan-resistome + 1))) + 
    ylab(expression(log[10](Core-resistome + 1)))
}
 
plot_core_pan_list_sup <- list()

for(j in 1:length(places)){
  plot_core_pan_list_sup[[j]] <- 
    pan_core_sum %>% 
    rename(tool_name = tools_labels) %>%
    filter(!habitat %in% c("human gut", "pig gut", "dog gut", "human nose", "wastewater",
                          "marine", "soil", "human oral", "human skin")) %>% 
    filter(
      place %in% places[j]) %>% 
    ggplot(aes(y = log10(Core+1), x = log10(Pan+1), shape = habitat)) +
    geom_line(aes(group = habitat), colour = "#666666", linewidth = 0.3, show.legend = FALSE) +
    geom_point(aes(fill = tool_name), color = "black", stroke = 0.15, size = 2.5, key_glyph = draw_key_hatch) +
    geom_point(data = ~ filter(.x, texture == "yes"), shape = 4, colour = "black",
               size = 1.8, stroke = 0.35, show.legend = FALSE) +
    facet_grid(place ~ ., scales = "free_x", space = "free_x") +
    scale_fill_manual(values = pal_tool, guide = "none") +
    scale_shape_manual(values = c(21, 22, 23, 24, 25), labels = hab_lab,
                       guide = guide_legend(override.aes = list(fill = "#666666", colour = "black"),
                                            title = "")) +
    theme1 + 
    theme(legend.position = "inside",
          legend.position.inside = c(0.001, 1.1),
          legend.justification = c(0, 1),
          legend.background = element_blank(),
          legend.key = element_blank()) + 
    xlim(c(min(log10(pan_core_sum$Pan+1)), max(log10(pan_core_sum$Pan+1)))) +
    ylim(c(0, max(log10(pan_core_sum$Core+1)))) +
    xlab(expression(log[10](Pan-resistome + 1))) + 
    ylab(expression(log[10](Core-resistome + 1)))
}

pan_core_new <- bind_rows(pan %>% 
                            mutate(metric = "Pan-resistome") %>% 
                            group_by(metric, habitat, tool, tools_labels, texture, tools_db, new_level) %>% 
                            summarise(size = sum(p_pan)) %>% ungroup(), 
          core %>% mutate(metric = "Core-resistome") %>% 
            group_by(metric, habitat, tool, tools_labels, texture, tools_db, new_level) %>% 
            summarise(size = n_distinct(centroid), .groups = "drop"))


pan_core_new <- pan_core_new %>% 
  group_by(metric,  tool, tools_labels, tools_db, texture, habitat) %>% 
  arrange(desc(size)) %>% 
  mutate(p = size / sum(size)) %>% 
  ungroup() %>% 
  mutate(metric = factor(metric, levels = c("Pan-resistome", "Core-resistome")))

pan_core_new <- pan_core_new %>% mutate(new_level = gsub(" beta-lactamase","", new_level))
#pan_core_new <- pan_core_new %>% mutate(new_level = gsub("cell wall ","cell wall\n", new_level))
pan_core_new <- pan_core_new %>% mutate(new_level = gsub("efflux pump", "efflux", new_level))
pan_core_new <- pan_core_new %>% mutate(new_level = gsub("rifampin inactivation enzyme","RIF-inact.", new_level))
pan_core_new <- pan_core_new %>% mutate(new_level = gsub("beta-lactam modulation resistance","beta-lactam mod.", new_level)) 
pan_core_new <- pan_core_new %>% mutate(new_level = gsub("MFS efflux pump","MFS\nefflux", new_level))


composition_pan_core_list <-  list()
thres_plot_pan <- 0.15

for(h in 1:length(EN)){
  cl_pan <- pan_core_new %>% 
    filter(habitat %in% EN[h]) %>% 
    group_by(habitat, metric, tool) %>% 
    filter(p > thres_plot_pan) %>%
    ungroup() %>% 
    distinct(new_level) %>% 
    pull(new_level)   
  
  composition_pan_core_list[[h]] <- pan_core_new %>% 
    group_by(habitat, metric) %>% 
    filter(new_level %in% cl_pan, habitat %in% EN[h]) %>% 
    ggplot(aes(x = tools_labels, y = new_level, fill = p)) +
    geom_tile(color = "white") + 
    geom_text(
      data = pan_core_new %>% 
        group_by(habitat, metric) %>% 
        filter(new_level %in% cl_pan, habitat %in% EN[h], p > 0.15),  # only these get labels
      aes(label = scales::percent(p, accuracy = 1),
          color = p >= 0.6),
      size = (general_size -1) / .pt, show.legend = F) +
    scale_color_manual(values = c( "black", "white")) +
    scale_y_discrete(labels = function(x) {
      out <- ifelse(
        x %in% c("rpoB", "sul", "van", "fos", "erm", "cat", 
                 "aph", "ant", "aac", "lnu", "nim","vat","mph","qnr"),
        paste0("italic('", x, "')"),
        ifelse(x %in% "abcF", "ABC-F",
               paste0("'", x, "'"))
      )
      parse(text = out)
    }) + 
    facet_grid(habitat ~ metric) +
    scale_fill_gradientn(
      colors = brewer.pal(9, "YlOrBr"),
      labels = percent_format(accuracy = 1),
      breaks = c(0, 0.33, 0.66, 1),
      limits = c(0, 1)) +
    theme1 +
    theme(legend.position = "bottom", panel.grid = element_blank())
} 



  composition_pan_core_list_sup <- pan_core_new %>% 
    filter(habitat %in% c("human oral", "human nose", "human vagina", "human skin",
                          "dog gut", "cat gut", "mouse gut", "wastewater", "freshwater")) %>% 
    group_by(habitat, metric) %>% 
    filter(p > thres_plot_pan) %>%
    ggplot(aes(x = tools_labels, y = new_level, fill = p)) +
    geom_tile(color = "white") + 
    geom_text(
      data = pan_core_new %>% 
        filter(habitat %in% c("human oral", "human nose", "human vagina", "human skin",
                              "dog gut", "cat gut", "mouse gut", "wastewater", "freshwater")) %>% 
        group_by(habitat, metric) %>% 
        filter(p > thres_plot_pan),  # only these get labels
      aes(label = scales::percent(p, accuracy = 1),
          color = p >= 0.6),
      size = (general_size -1) / .pt, show.legend = F) +
    scale_color_manual(values = c( "black", "white")) +
    scale_y_discrete(labels = function(x) {
      out <- ifelse(
        x %in% c("rpoB", "sul", "van", "fos", "erm", "cat", 
                 "aph", "ant", "aac", "lnu", "nim","vat","mph","qnr"),
        paste0("italic('", x, "')"),
        ifelse(x %in% "abcF", "ABC-F",
               paste0("'", x, "'"))
      )
      parse(text = out)
    }) + 
    facet_grid(habitat ~ metric, scales = "free_y", space = "free_y") +
    scale_fill_gradientn(
      colors = brewer.pal(9, "YlOrBr"),
      labels = percent_format(accuracy = 1),
      breaks = c(0, 0.33, 0.66, 1),
      limits = c(0, 1)) +
    theme1 +
    theme(legend.position = "bottom", panel.grid = element_blank())



composition_pan_core <- (composition_pan_core_list[[1]] + 
                           ggtitle("b") + ylab ("") + xlab("") + 
                           theme(axis.text.x = element_blank(), legend.position = "none")) / 
  (composition_pan_core_list[[9]] + xlab("") + ylab ("")  + 
     theme(axis.text.x = element_blank() ,  
                 legend.position = "none", 
   strip.text.x = element_blank())) / 
  (composition_pan_core_list[[13]] + ylab ("") + theme(axis.text.x = element_blank(), 
                   legend.position = "none", 
                   strip.text.x = element_blank())) / 
  (composition_pan_core_list[[11]] + xlab("") + ylab ("") + theme( strip.text.x = element_blank(), legend.position = "none")) +
  patchwork::plot_layout(heights = c(4,4,4,4))

composition_pan_core_sup1 <- (composition_pan_core_list[[2]] + 
                           ggtitle("") + ylab ("") + xlab("") + 
                           theme(axis.text.x = element_blank(), legend.position = "none")) / 
  (composition_pan_core_list[[3]] + xlab("") + ylab ("")  + 
     theme(axis.text.x = element_blank() ,  
           legend.position = "none", 
           strip.text.x = element_blank())) / 
  (composition_pan_core_list[[4]] + xlab("")  + ylab ("") + theme(axis.text.x = element_blank(), 
                                                       legend.position = "none", 
                                                       strip.text.x = element_blank())) / 
  (composition_pan_core_list[[5]] + xlab("") +ylab ("") + theme( strip.text.x = element_blank(), legend.position = "none")) +
  patchwork::plot_layout(heights = c(4,4,4,4))

  
composition_pan_core_sup2 <- (composition_pan_core_list[[6]] + 
                                ggtitle("") + ylab ("") + xlab("") + 
                                theme(axis.text.x = element_blank(), legend.position = "none")) / 
  (composition_pan_core_list[[7]] + xlab("") + ylab ("")  + 
     theme(axis.text.x = element_blank() ,  
           legend.position = "none", 
           strip.text.x = element_blank())) / 
  (composition_pan_core_list[[8]] + xlab("") + ylab ("") + theme(axis.text.x = element_blank(), 
                                                      legend.position = "none", 
                                                      strip.text.x = element_blank())) / 
  (composition_pan_core_list[[10]] + xlab("") + ylab ("") + theme(axis.text.x = element_blank(), 
                                                       legend.position = "none", 
                                                       strip.text.x = element_blank())) /
  (composition_pan_core_list[[12]] + xlab("") + ylab ("") + theme( strip.text.x = element_blank(), legend.position = "none")) +
  patchwork::plot_layout(heights = c(4,4,4,4, 4))

p_core_pan_1 <- (plot_core_pan_list[[1]] + xlab("") + ylab("") + ggtitle("a")) /
  (plot_core_pan_list[[2]] + xlab("") + ylab(expression(log[10]("Core-resistome" + 1)))) /
  (plot_core_pan_list[[3]] + xlab(expression(log[10]("Pan-resistome" + 1))) + ylab("")) 


plot_core_pan_list_sup_1 <- (plot_core_pan_list_sup[[1]] + xlab("") + ylab("") + ggtitle("")) /
  (plot_core_pan_list_sup[[2]] + xlab("") + ylab(expression(log[10]("Core-resistome" + 1)))) /
  (plot_core_pan_list_sup[[3]] + xlab(expression(log[10]("Pan-resistome" + 1))) + ylab("")) 


#  / 
#(wrap_elements(full = leg)) patchwork::plot_layout(heights = c(5,5,5,1))

p4_new <- ((p_core_pan_1 | composition_pan_core) +
    patchwork::plot_layout(widths = c(3,5)))/
(((wrap_elements(full = leg)) |  wrap_elements(g_legend(co_marine))) +
   patchwork::plot_layout(widths = c(5,3))) +
   patchwork::plot_layout(heights = c(9,1))

ggsave("code_R_analysis/output_plots/fig4.svg", p4_new, width = 180, height = 180, unit = "mm")

#ggsave("code_R_analysis/output_plots/sup_all_pan_core_other.svg", plot_core_pan_list_sup_1, width = 90, height = 90, unit = "mm")

composition_pan_core_sup1 |composition_pan_core_sup2 | plot_core_pan_list_sup_1
#ggsave("code_R_analysis/output_plots/sup_comp_pan_core_human.svg", composition_pan_core_sup1, width = 130, height = 180, unit = "mm")
#ggsave("code_R_analysis/output_plots/sup_comp_pan_core_other.svg", composition_pan_core_sup2, width = 130, height = 180, unit = "mm")


u3_legend <- umap3.1 %>%
  mutate(al = place == "External") %>% 
  ggplot(aes(u1, u2, colour = tools_labels, fill = tools_labels, shape = texture, size = al, alpha = al)) +
  geom_point(color = "black", stroke = .2, key_glyph = draw_key_point_line) +
  scale_size_manual(values = c("TRUE" = 0.7, "FALSE" = 0.05), guide = "none") +
  scale_colour_manual(values = pal_tool, name = "") +
  scale_fill_manual(values = pal_tool, name = "") +
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.1), guide = "none") +
  guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1), nrow = 2)) +
  theme1 + 
  theme(plot.subtitle = element_text(size = general_size, color = "black")) +
  xlab("UMAP1") +
  ylab("UMAP2") +
  labs(subtitle = "External habitats") + 
  scale_shape_manual(values = c("no" = 21, "yes" = 24)) + 
  guides(fill = 
           guide_legend(nrow = 4, override.aes = list(shape = shape_tools,
                                                      color = "black", stroke = 0.4)), shape = "none")  


leg2 <- g_legend(u3_legend + ylab("") + xlab("") + theme(legend.position = "bottom"))

ggsave("code_R_analysis/output_plots/sup_pan_core.svg",
       ((((plot_core_pan_list_sup_1) / wrap_elements(full = leg2)) +
           patchwork::plot_layout(heights = c(5/3, 5/3, 5/3, 1))) |
          (composition_pan_core_list_sup + labs(fill = "") + xlab("") + ylab(""))) +
         patchwork::plot_layout(widths = c(3, 5)), width = 180, height = 180, unit = "mm")

#########

# basic theme for the overlap plot 

theme5 <- theme(
  legend.position = "none",
  legend.text = element_text(size = general_size ),
  panel.grid.major.y = element_blank(),
  panel.grid.minor.y = element_blank(),
  panel.grid.minor.x = element_blank(),
  plot.margin = margin(0, 0, 0, 0, unit = "pt"),
  legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
  legend.margin = margin(0, 0, 0, 0, unit = "pt"),
  panel.spacing = unit(0, "pt"),
  title = element_text(size = general_size + 2, face = "bold"),
  axis.title = element_text(size = general_size + 1, face = "bold"),
  axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = general_size),
  axis.text.y = element_text(size = general_size),
  panel.background = element_rect(colour = "black", fill = NA),
  strip.text = element_text(size = general_size ),
  strip.background = element_blank(),
  panel.grid.major.x = element_line(colour = "black", linewidth = 0.1))


# modify the recall_fnr object to add texture, change the tools' name 

d1 <- recall_fnr %>% 
  mutate(tool_ref = factor(as.character(tool_ref), 
                           levels = tools_levels)) %>% 
  mutate(texture = ifelse(tool_comp %in% tools_texture, "yes", "no")) %>%
  mutate(texture2 = ifelse(tool_ref %in% tools_texture, "yes", "no")) %>%
  mutate(facet_var = gsub("-", "-\n", tool_ref)) %>%
  mutate(facet_var = gsub(" ", "\n", facet_var)) %>%
  mutate(facet_var = gsub("Plus", "\n Plus", facet_var)) %>%
  mutate(facet_var = fct_reorder(facet_var, as.numeric(tool_ref)))  %>% 
  mutate(val = csc, d = "csc") %>% 
  mutate(tools_labels_comp = factor(tools_labels[tool_comp], levels = tools_labels_factor),
         texture_comp = ifelse(tool_comp %in% tools_texture, "yes", "no"),
         tools_db_comp = factor(tools_db[tool_comp], levels = tools_db_factor)) %>% 
  mutate(tools_labels_comp = factor(gsub("-\n","",tools_labels_comp), levels = gsub("-\n","", levels(tools_labels_comp)))) %>%
  mutate(tools_labels_ref = factor(tools_labels[tool_ref], levels = tools_labels_factor))


# change the name of some classes for plotting

d1 <- d1 %>% 
  filter(tool_ref %in% basic_tools, tool_comp %in% basic_tools) %>%
  mutate(gene_class = gsub(" beta-lactamase","", gene_class)) %>%
  mutate(gene_class = gsub("MFS efflux pump","efflux", gene_class)) %>%
  mutate(gene_class = gsub("efflux pump","efflux", gene_class)) %>% 
  group_by(tool_ref,tools_labels_ref, gene_class) %>% mutate(n_obs = n()) %>% 
  mutate(n_obs = paste0('n = ',n_obs))


# plot the csc of top_cso classes

cs11 <- d1 %>% 
  filter(gene_class %in% c(top_cso,"class A", "class B", "class C", "class D", "efflux")) %>% 
  ggplot(aes(x = csc*100, y = 0)) + 
  geom_boxplot_pattern(aes(fill = tool_ref, pattern = texture2),
                       position = position_dodge2(preserve = "single", width = 0.3, padding = 0), 
                       width = 1.3, pattern_color = "black", pattern_fill = "black", pattern_density = 0.000000001,
                       pattern_spacing = 0.2,
                       pattern_size =  0.3, color = "black", outliers = FALSE, outlier.shape = NA, linewidth = 0.15) +
  scale_pattern_manual(values = c('no' = 'none', 'yes' = 'stripe')) +
  geom_text(aes(x = 50, y = 0.9, label = n_obs), size = general_size / .pt) + 
  facet_grid(gene_class ~ tools_labels_ref, scales = "free_y", space = "free",
             labeller = labeller(
               gene_class = as_labeller(function(x) {
                 out <- ifelse(
                   x %in% c("rpoB", "van", "fos", "erm", "cat", "aph", "ant", 
                            "aac", "lnu", "nim", "vat", "mph", "qnr"),
                   paste0("italic('", x, "')"),
                   ifelse(x == "abcF", "ABC-F",
                          paste0("'", x, "'"))
                 )
                 return(out)
               }, default = label_parsed)
             )) +
  scale_fill_manual(values = pal_10_q)  + 
  xlab("Class-specific coverage (%)") +
  ylab("ARG class") + 
  theme_minimal() +
  theme5 +
  theme( panel.grid = element_blank(),
         strip.text.x = element_text(size = general_size, vjust = 0, hjust = 0.5 , angle = 0),
         strip.text.y = element_text(size = general_size, vjust = 0.5, hjust = 0, angle = 0),
         panel.spacing = unit(5, "pt")) +
  scale_y_discrete(labels = function(x) {
    out <- ifelse(
      x %in% c("rpoB", "van", "fos", "erm", "cat", 
               "aph", "ant", "aac", "lnu", "nim","vat","mph","qnr"),
      paste0("italic('", x, "')"),
      ifelse(x %in% "abcF", "ABC-F",
             paste0("'", x, "'"))
    )
    parse(text = out)
  })


ggsave("code_R_analysis/output_plots/fig5.svg", cs11, width = 180, height = 100, unit = "mm")


#######################################################################
#######################################################################

# calculate the intra-tool overalap and the overlap between tools with the same database

# basic theme for the plots 

theme_overlap <- theme_minimal() + 
  theme(legend.position = "bottom",
        panel.border = element_blank(),
        axis.text.x = element_text(size = general_size),
        axis.text.y = element_text( size = general_size),
        strip.text = element_text(size = general_size, face = "bold"),
        plot.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.margin = margin(0, 0, 0, 0, unit = "pt"),
        panel.spacing = unit(0, "pt"),
        title = element_text(size = general_size + 2, face = "bold"),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title = element_text(size = general_size , face = "bold"),
        panel.background = element_blank(),
        legend.text = element_text(size = general_size))


# for RGI, the overlap is between amino acid and nucleotide and DIAMOND and BLAST

# create the overlap for RGI

lst2 <- readRDS("code_R_analysis/output_abundance_diversity_resistome/results_tools_all_GMGC.rds")

rgi_query <- data.frame(gene = unique(c(lst2$rgi.blast$query, lst2$rgi.diamond$query, lst2$rgi.diamond.prot$query)))
rgi_query <- rgi_query %>% 
  mutate(dataset = 
           ifelse(gene %in% intersect(intersect(lst2$rgi.blast$query, lst2$rgi.diamond$query), lst2$rgi.diamond.prot$query), "all setups",
                  ifelse(gene %in% intersect(lst2$rgi.blast$query, lst2$rgi.diamond$query), "BLAST nt and DIAMOND (nt or aa)",
                         ifelse(gene %in% intersect(lst2$rgi.diamond.prot$query, lst2$rgi.diamond$query), "only\nDIAMOND\nnt and/or aa",
                                ifelse(gene %in% intersect(lst2$rgi.diamond.prot$query, lst2$rgi.blast$query), "BLAST nt and DIAMOND (nt or aa)",
                                       ifelse(gene %in% lst2$rgi.diamond.prot$query, "only\nDIAMOND\nnt and/or aa",
                                              ifelse(gene %in% lst2$rgi.diamond$query, "only\nDIAMOND\nnt and/or aa",
                                                     ifelse(gene %in% lst2$rgi.blast$query, "only BLAST\nnt", NA)))))))) %>%
  mutate(dataset = factor(dataset, levels = c("all setups", "only\nDIAMOND\nnt and/or aa", 
                                              "BLAST nt and DIAMOND (nt or aa)",  "only BLAST\nnt")))

# modify to plot in the right order and with the right label 

rgi_query2 <- rgi_query %>% 
  group_by(dataset) %>% summarise(n = n()) %>% mutate(p = n/sum(n)) %>% 
  ungroup() %>%
  arrange(max(as.numeric(dataset)) - as.numeric(dataset)) %>%
  mutate(N = cumsum(n)) %>%
  mutate(diff_by1 = (N - lag(N, n = 1, default = 0))/2)

# plot the overlap for RGI 

rgi_over <- rgi_query2%>%
  ggplot(aes(y = "RGI", x = n, fill = dataset)) + 
  geom_col() + 
  scale_fill_manual(values = brewer.pal(8, "Dark2")[2:8], 
                    labels = function(x) str_replace(x, " and ", " and\n"))  +
  geom_text(data = rgi_query2 %>% filter(p>0.02),
            aes(x = N - diff_by1, label = scales::percent(p, accuracy = 1)), 
            size = general_size / .pt) +
  xlab("ARGs") + 
  guides(fill = guide_legend(nrow = 2)) +
  scale_x_continuous(labels = scales::comma) +
  ylab("") +
  ggtitle("RGI") + 
  labs(fill = "" ) +
  scale_y_discrete(labels = function(x) {
    x <- gsub(" ", "\n", x)
    x}, expand = 0) +
  theme_overlap + 
  theme(plot.margin = margin(5, 5, 5, 5, unit = "pt"))

# create a function for overlap between only two pipelines

plot_overlaps <- function(x, y, nxy, nx, ny, nam=""){
  db_query <- data.frame(gene = unique(c(x, y)))
  db_query <- db_query %>% mutate(dataset = ifelse(gene %in% intersect(x, y), nxy,
                                                   ifelse(gene %in% x, nx, ny))) %>% 
    mutate(dataset = factor(dataset, levels = c(nxy, nx, ny)))
  
  db_query2 <- db_query %>% 
    group_by(dataset) %>% summarise(n = n()) %>% mutate(p = n/sum(n)) %>% 
    ungroup() %>%
    arrange(max(as.numeric(dataset)) - as.numeric(dataset)) %>%
    mutate(N = cumsum(n)) %>%
    mutate(diff_by1 = (N - lag(N, n = 1, default = 0))/2)
  
  db_query_plot <- db_query2 %>% 
    ggplot(aes(y = nam, x = n, fill = dataset)) + 
    geom_col() + 
    scale_fill_manual(values = brewer.pal(8, "Dark2")[c(2,5,3,6)]) +#, 
    geom_text(data = db_query2 %>% filter(p>0.02),
              aes(x = N - diff_by1, label = scales::percent(p, accuracy = 1)), 
              size = general_size / .pt) +
    xlab("") + 
    guides(fill = guide_legend(nrow = 2)) +
    scale_x_continuous(labels = scales::comma) +
    ylab("") +
    ggtitle("") + 
    labs(fill = "" ) +
    scale_y_discrete(labels = function(x) {
      x <- gsub(" ", "\n", x)
      x}, expand = 0) +
    theme_overlap
  
  return(db_query_plot)
}

# overlap for deeparg

deeparg_over <- plot_overlaps(lst2$deeparg.norm$query, lst2$deeparg.norm.prot$query, "nt and aa", "only nt", "only aa", nam = "DeepARG")

# overlap for fargene

fargene_over <- plot_overlaps(lst2$fargene$query, lst2$fargene.prot$query, "nt and aa", "only nt", "only aa", nam = "fARGene")

# overlap for amrfinder

amrfinder_over <- plot_overlaps(lst2$amrfinder.norm$query, lst2$amrfinder.norm.prot$query, "nt and aa", "only nt", "only aa", nam = "AMRFinder-\nPlus")


#  merge rgi, fargene, deeparg, amrfinder

plot_overlaps_same_tool <- 
  ((rgi_over + xlab("") + ggtitle("a")) / 
     (deeparg_over + theme(legend.position = "none") + ggtitle("c")) /
     (fargene_over + theme(legend.position = "none")) /
     (amrfinder_over + xlab("ARGs"))) + 
  patchwork::plot_layout(heights = c(1, 1, 1, 1)) & 
  theme(plot.margin = margin(5, 5, 5, 5, unit = "pt"))


#  plot card, ncbi, resfinder

card_over <- plot_overlaps(lst2$rgi.diamond.prot$query, lst2$abricate.card.norm$query, "RGI and\nABRicate", "only RGI", "only ABRicate", nam = "CARD")
ncbi_over <- plot_overlaps(lst2$amrfinder.norm.prot$query, lst2$abricate.ncbi.norm$query, "AMRFinderPlus\nand ABRicate", "only AMRFinderPlus", "only ABRicate", nam = "NCBI")
resfinder_over <- plot_overlaps(lst2$resfinder.norm$query, lst2$abricate.resfinder.norm$query, "ResFinder\nand ABRicate", "only ResFinder", "only ABRicate", nam = "ResFinder")

#  merge

plot_db <- 
  (card_over + ggtitle("b")) / 
  ncbi_over / 
  resfinder_over + xlab("ARGs") & 
  theme(plot.margin = margin(5, 5, 5, 5, unit = "pt"))


# plot for megares vs other abricate 

other_abricate <- unique(c(lst2$abricate.argannot.norm$query, lst2$abricate.card.norm$query, 
                           lst2$abricate.ncbi.norm$query, lst2$abricate.resfinder.norm$query))

abricate_query <- data.frame(rbind(lst2$abricate.megares.norm[,c("tool", "query")], 
                                   lst2$abricate.argannot.norm[,c("tool", "query")], 
                                   lst2$abricate.card.norm[,c("tool", "query")], 
                                   lst2$abricate.ncbi.norm[,c("tool", "query")], 
                                   lst2$abricate.resfinder.norm[,c("tool", "query")]))

ntools_abricate <- abricate_query %>% 
  group_by(query) %>% 
  summarise(n_tools = n_distinct(tool))

# modify the labels for intercepts with megares 

ntools_abricate <- ntools_abricate %>% mutate(s = ifelse(query %in% lst2$abricate.megares.norm$query, "MEGARes","other"))
ntools_abricate <- ntools_abricate %>% mutate(s = ifelse(n_tools == 2 & s == "MEGARes" & query %in% lst2$abricate.argannot.norm$query, "MEGARes and\nARGANNOT", s))
ntools_abricate <- ntools_abricate %>% mutate(s = ifelse(n_tools == 2 & s == "MEGARes" & query %in% lst2$abricate.ncbi.norm$query, "MEGARes and\nNCBI", s))
ntools_abricate <- ntools_abricate %>% mutate(s = ifelse(n_tools == 2 & s == "MEGARes" & query %in% lst2$abricate.resfinder.norm$query, "MEGARes and\nResFinder", s))
ntools_abricate <- ntools_abricate %>% mutate(s = ifelse(n_tools == 2 & s == "MEGARes" & query %in% lst2$abricate.card.norm$query, "MEGARes and\nCARD", s))
ntools_abricate <- ntools_abricate %>% mutate(s = ifelse(n_tools == 1 & s != "MEGARes", "other \n single\n tool", s))
ntools_abricate <- ntools_abricate %>% mutate(s = ifelse(n_tools > 1 & !grepl("MEGARes", s), "other \n>1 tools", s))
ntools_abricate <- ntools_abricate %>% mutate(s = ifelse(n_tools > 2 & grepl("MEGARes", s), "MEGARes and\n >1 tools", s))
ntools_abricate <- ntools_abricate %>% mutate(s = ifelse(s == "MEGARes", "only MEGARes", s))

ntools_abricate %>% ungroup %>% group_by(n_tools, s) %>% summarise(n = n())

# factorize the labels 

ntools_abricate1 <- ntools_abricate %>% mutate(s = factor(s, levels = 
                                                            c("only MEGARes","MEGARes and\nCARD","MEGARes and\nResFinder",
                                                              "MEGARes and\nNCBI","MEGARes and\nARGANNOT",
                                                              "MEGARes and\n >1 tools", "other \n>1 tools","other \n single\n tool")))

ntools_abricate2 <- ntools_abricate %>% 
  mutate(s = ifelse(!grepl("MEGARes",s),"not in MEGARes",s)) %>%
  mutate(s = ifelse(grepl("MEGARes",s) & n_tools > 1,"MEGARes and \n other tools",s)) %>%
  mutate(s = factor(s, levels = 
                      c("only MEGARes","MEGARes and \n other tools","not in MEGARes")))

abricate_p1 <- ntools_abricate2 %>% ungroup %>% group_by(s) %>% 
  summarise(n = n()) %>% mutate(p = n/sum(n)) %>%
  ungroup() %>%
  arrange(max(as.integer(s)) - as.integer(s)) %>%
  mutate(N = cumsum(n)) %>%
  mutate(diff_by1 = (N - lag(N, n = 1, default = 0))/2) 

# abricate plot 

abricate_plot <- abricate_p1 %>%
  ggplot(aes(y = "ABRicate", x = n, fill = s)) + 
  geom_col() + 
  scale_fill_manual(values = as.vector(pal_10_complete[-1]), 
                    labels = function(x) str_replace(x, " and ", " and\n"))  +
  xlab("ARGs") + 
  geom_text( data = abricate_p1,
             aes(x = N - diff_by1, label = scales::percent(p, accuracy = 1)),  
             size = general_size / .pt, show.legend = F) +
  scale_color_manual(values = c( "black", "white")) +
  guides(fill = guide_legend(nrow = 3)) +
  ylab("") +
  ggtitle("") + 
  labs(fill = "" ) +
  scale_x_continuous(labels = scales::comma) + 
  scale_y_discrete(labels = function(x) {
    x <- gsub(" ", "\n", x)
    x}, expand = 0) +
  theme_overlap + 
  theme(plot.margin = margin(5, 5, 5, 5, unit = "pt"))


# megares plot 

meg_plot_1 <- ntools_abricate1 %>% filter(grepl("MEGARes and", s)) %>% 
  ungroup %>% group_by(s) %>% 
  summarise(n = n()) %>% mutate(p = n/sum(n)) %>%
  ungroup() %>%
  arrange(max(as.integer(s)) - as.integer(s)) %>%
  mutate(N = cumsum(n)) %>%
  mutate(diff_by1 = (N - lag(N, n = 1, default = 0))/2) 

megares_plot <- meg_plot_1 %>%
  ggplot(aes(y = "MEGARes", x = n, fill = s)) + 
  geom_col() + 
  scale_fill_manual(values = as.vector(pal_10_complete[-1]), 
                    labels = function(x) str_replace(x, " and ", " and\n"))  +
  xlab("ARGs") + 
  guides(fill = guide_legend(nrow = 3)) +
  geom_text( data = meg_plot_1 %>% filter(p > 0.025),
             aes(x = N - diff_by1, label = scales::percent(p, accuracy = 1)),  
             size = general_size / .pt, show.legend = F) +
  scale_color_manual(values = c( "black", "white")) +
  ylab("") +
  ggtitle("") + 
  labs(fill = "" ) +
  scale_x_continuous(labels = scales::comma) + 
  scale_y_discrete(labels = function(x) {
    x <- gsub(" ", "\n", x)
    x}, expand = 0) +
  theme_overlap + 
  theme(plot.margin = margin(5, 5, 5, 5, unit = "pt"))

# merge megares and abricate plots

plot_megares <- (abricate_plot + ggtitle("c")) / megares_plot


overlap_all_plots <- (plot_overlaps_same_tool / abricate_plot + ggtitle("d")) + patchwork::plot_layout(heights = c(1,1,1,1,1)) |
  (plot_db / megares_plot + ggtitle("e")) + patchwork::plot_layout(heights = c(1,1,1,1))

ggsave("code_R_analysis/output_plots/overlaps_all_plots.svg", overlap_all_plots, width = 180, height = 180, unit = "mm")



####


# create the plots of richness

rich_plots <- list()
set.seed(2026)

for(j in 1:length(EN)){
  
  # get the number of samples per habitat
  df_plot <- abundance_tool_sample %>%
    filter(tool %in% basic_tools) %>% 
    filter(habitat %in% EN[j]) %>% 
    group_by(habitat) %>% 
    mutate(N = n_distinct(sample),
           richness = richness)
  
  # get the limits for the jitter 
  df_jitter <- df_plot %>% 
    ungroup() %>%
    group_by(habitat, tool) %>% 
    filter(richness < quantile(richness, 0.75) + 1.5*IQR(richness))
  
  # get the maximum for the limits of the plot
  
  lims_richness_max = lims_richness %>%
    filter(habitat %in% EN[j]) %>%
    filter(tool %in% basic_tools) %>% 
    ungroup() %>% summarise(max(w2)) %>% pull()
  
  # plot
  
  rich_plots[[j]] <- lims_richness %>%
    mutate(habitat_label = paste0(habitat, "\n(n = ", scales::comma(df_plot$N[1]), ")")) %>%
    filter(habitat %in% EN[j]) %>%
    filter(tool %in% basic_tools) %>% 
    ggplot(aes(x = tools_labels, y=median, fill = tools_db))  + 
    geom_boxplot_pattern(aes(ymin = w1, lower = q25, middle = median, upper = q75, ymax = w2, pattern = texture), stat = "identity",
                         position = position_dodge2(preserve = "single", width = 0.3, padding = 0), 
                         width = .8, pattern_color = "black", pattern_fill = "black", pattern_density = 0.00001,
                         pattern_spacing = 0.15,
                         pattern_size =  0.4, color = "black", outliers = FALSE, outlier.shape = NA, linewidth = 0.15) +
    scale_x_discrete(expand = expansion(add = 1)) +
    geom_jitter(data = df_jitter %>%
                  slice_sample(n = min(100, df_jitter$N[1])), 
                aes(y = abundance),
                width = 0.35, size = 0.4, alpha = 0.1) + 
    facet_grid(habitat_label ~ tools_db , scales = "free_x", space = "free_x") +
    scale_y_continuous(labels = scales::comma, limits = c(0, lims_richness_max + 2)) + 
    scale_fill_manual(values = pal_7) +
    scale_pattern_manual(values = c('no' = 'none', 'yes' = 'stripe')) +
    theme1 + 
    theme(panel.border = element_blank(),
          strip.text.x = element_text(size = general_size, angle = 90, vjust = 0.5, hjust = 0),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          plot.margin = margin(5.5, 5.5, 5.5, 5.5, unit = "pt")) 
}


# create the three panels of the figures 

sup_rich_left <- 
  ((rich_plots[[1]] + theme(axis.text.x = element_blank(), strip.text.x = element_blank()) + xlab("") + ylab("")) / 
     (rich_plots[[2]] + theme(axis.text.x = element_blank(), strip.text.x = element_blank()) + xlab("") + ylab("")) / 
     (rich_plots[[3]] + theme(axis.text.x = element_blank(), strip.text.x = element_blank()) + xlab("") + ylab("Richness")) /
     (rich_plots[[4]] + theme(axis.text.x = element_blank(), strip.text.x = element_blank()) + xlab("") + ylab("")) /
     (rich_plots[[5]] + theme(strip.text.x = element_blank()) + xlab("") + ylab(""))) +
  patchwork::plot_layout(heights = c(1, 1, 1, 1, 1))

sup_rich_mid <- 
  (patchwork::plot_spacer() /
     (rich_plots[[6]] + theme(axis.text.x = element_blank(), strip.text.x = element_blank()) + xlab("") + ylab("")) / 
     (rich_plots[[7]] + theme(axis.text.x = element_blank(), strip.text.x = element_blank()) + xlab("") + ylab("")) /
     (rich_plots[[8]] + theme(axis.text.x = element_blank(), strip.text.x = element_blank()) + xlab("") + ylab("")) / 
     (rich_plots[[9]] + theme(strip.text.x = element_blank()) + xlab("") + ylab(""))) +
  patchwork::plot_layout(heights = c(1, 1, 1, 1, 1))


sup_rich_right <- 
  (patchwork::plot_spacer() /
     (rich_plots[[10]] + theme(axis.text.x = element_blank(), strip.text.x = element_blank()) + xlab("") + ylab("")) / 
     (rich_plots[[11]] + theme(axis.text.x = element_blank(), strip.text.x = element_blank()) + xlab("") + ylab("")) /
     (rich_plots[[12]] + theme(axis.text.x = element_blank(), strip.text.x = element_blank()) + xlab("") + ylab("")) / 
     (rich_plots[[13]] + theme(strip.text.x = element_blank()) + xlab("") + ylab(""))) +
  patchwork::plot_layout(heights = c(1, 1, 1, 1, 1))


# merge the three panels

sup_rich <- sup_rich_left | sup_rich_mid | sup_rich_right 

ggsave("code_R_analysis/output_plots/sup_richness.svg", sup_rich, width = 180, height = 180, unit = "mm")



###

# create the plots of the abundance per class and pipeline for supplementary

abu_class_plots <- list()
df_class_summary_plots <- list()

for( j in 1:length(EN)){
  
  # get the most relevant gene classes
  all_top_abundance <- unique(
    levels_abundance_div %>% 
      filter(habitat %in% EN[j]) %>% 
      select(gene) %>% distinct() %>% slice_head(n=10) %>%
      ungroup() %>% pull())
  
  # add the "other" class
  df_abundance_class_all_sup <- abundance_class %>% 
    filter(habitat %in% EN[j]) %>%
    mutate(gene = ifelse(gene %in% all_top_abundance, gene, "other")) %>% 
    mutate(texture = ifelse(tool %in% tools_texture, "yes", "no")) %>%
    mutate(gene = factor(gene, levels = c(all_top_abundance, "other"))) %>%
    ungroup() %>% 
    group_by(habitat, sample, tool, gene) %>%
    summarise(abundance = sum(abundance), richness = sum(richness), texture = texture[1]) %>%
    mutate(tools_labels = factor(tools_labels[tool], levels = tools_labels[!duplicated(tools_labels)]),
           texture = ifelse(tool %in% tools_texture, "yes", "no"),
           tools_db = factor(tools_db[tool], levels = tools_db[!duplicated(tools_db)]))
  
  # change the name of the labels
  gene_levels_all <- all_top_abundance
  gene_levels_all <-  gsub(" beta-lactamase","", gene_levels_all)
  gene_levels_all <- gsub("cell wall ","cell\nwall\n", gene_levels_all)
  gene_levels_all <- gsub("MFS efflux pump","MFS\nefflux", gene_levels_all)
  gene_levels_all <- gsub("efflux pump","efflux", gene_levels_all)
  gene_levels_all <- gsub("target-modifying enzyme","target-\n modifying\n enzyme", gene_levels_all)
  gene_levels_all <- gsub("molecular bypass","molecular\n bypass", gene_levels_all)
  gene_levels_all <- gsub("self-resistance","self-\n resistance", gene_levels_all)
  gene_levels_all <- gsub("variant or mutant","variant or\n mutant", gene_levels_all)
  gene_levels_all <- gsub("rifampin inactivation enzyme","rifampin\n inactivation\n enzyme", gene_levels_all)
  gene_levels_all <- gsub("beta-lactam modulation resistance","beta-lactam\n modulation\n resistance", gene_levels_all)
  
  # change the name of the labels
  df_abundance_class_all_sup <- df_abundance_class_all_sup  %>% 
    mutate(gene = gsub(" beta-lactamase","", gene)) %>%
    mutate(gene = gsub("cell wall ","cell\nwall\n", gene)) %>%
    mutate(gene = gsub("MFS efflux pump","MFS\nefflux", gene)) %>%
    mutate(gene = gsub("efflux pump","efflux", gene)) %>%
    mutate(gene = gsub("target-modifying enzyme","target-\n modifying\n enzyme", gene)) %>%
    mutate(gene = gsub("molecular bypass","molecular\n bypass", gene)) %>%
    mutate(gene = gsub("self-resistance","self-\n resistance", gene)) %>%
    mutate(gene = gsub("variant or mutant","variant or\n mutant", gene)) %>%
    mutate(gene = gsub("rifampin inactivation enzyme","rifampin\n inactivation\n enzyme", gene)) %>%
    mutate(gene = gsub("beta-lactam modulation resistance","beta-lactam\n modulation\n resistance", gene)) %>%
    mutate(gene = factor(gene, levels = c(gene_levels_all,"other"))) 
  
  # get the limits of the plots
  df_class_summary_plots[[j]] <- df_abundance_class_all_sup %>%
    ungroup() %>%
    group_by(tool, habitat, gene,texture, tools_labels, tools_db) %>% 
    summarise(
      q50 = ifelse(quantile(abundance, 0.5) < 0, 0, quantile(abundance, 0.5)),
      q25 = ifelse(quantile(abundance, 0.25) < 0, 0, quantile(abundance, 0.25)),
      q75 = ifelse(quantile(abundance, 0.75) < 0, 0, quantile(abundance, 0.75)),
      w1 = ifelse(quantile(abundance, 0.25) - 1.5*IQR(abundance)<0,0,
                  quantile(abundance, 0.25) - 1.5*IQR(abundance)),
      w2 = ifelse(quantile(abundance, 0.75) + 1.5*IQR(abundance)<0,0,
                  quantile(abundance, 0.75) + 1.5*IQR(abundance)),
      n = n_distinct(sample)) 
  
  # plot
  abu_class_plots[[j]] <- df_class_summary_plots[[j]] %>% 
    filter(tool %in% basic_tools) %>% 
    mutate(habitat_label = paste0(habitat, "\n(n = ", scales::comma(n[1]), ")")) %>%
    ggplot(aes(y = fct_rev(tool), fill = tool, pattern = texture)) +
    geom_boxplot_pattern( 
      aes(xmin = q25, xlower = q25, xmiddle = q50, xupper = q75, xmax = q75, pattern = texture),
      stat = "identity", 
      position = position_dodge2(preserve = "single", width = 0.5, padding = 0), 
      width = 1, pattern_color = "black", pattern_fill = pattern_fill, pattern_spacing = 0.07,
      pattern_density = 0.15,
      pattern_size =  0.07, color = "black", outliers = FALSE, outlier.shape = NA,
      linewidth = 0.1) +
    scale_fill_manual(values = pal_10_q[1:10], labels = 
                        gsub("ABRicate-", "ABRicate-\n", 
                             gsub("AMRFinderPlus", "AMRFinder-\nPlus", 
                                  gsub("RGI-DIAMOND", "RGI", basic_tools)))) +
    scale_pattern_manual(values = c('no' = 'none', 'yes' = 'stripe')) +
    facet_grid(gene ~ habitat_label, scales = "free", space = "free",
               labeller = labeller(
                 gene = as_labeller(function(x) {
                   out <- ifelse(
                     x %in% c("rpoB", "van", "fos", "erm", "cat", "aph", "ant", "aac", "lnu", "nim", "vat", "mph", "qnr"),
                     paste0("italic('", x, "')"),
                     ifelse(x == "abcF", "ABC-F",
                            paste0("'", x, "'"))
                   )
                   return(out)
                 }, default = label_parsed)
               )) + 
    xlab("Relative abundance\n(aligned reads per million)") +
    ylab("") + 
    labs(fill = "") + 
    scale_x_continuous(labels = scales::comma, expand = c(0,0)) +
    ggtitle("") +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
    theme1 + 
    theme(strip.text.y = element_text(size = general_size, vjust = 0.5, hjust = 0, angle = -90),
          legend.position = "none",
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.border =  element_blank(),
          axis.text.y = element_blank(),
          plot.margin = margin(5, 0, 0, 0, unit = "pt")) +
    guides(fill = guide_legend(
      override.aes = list(
        pattern = c("none", "stripe","none","stripe","none","stripe",
                    rep("none", 4)),
        fill  = pal_10_q[1:10])), pattern = "none")
  
}

# create the three figures

# human
sup_class_1 <- 
  (((abu_class_plots[[2]] +  xlab("") + ylab("")) |
      (abu_class_plots[[3]] +   ylab("")) |  
      (abu_class_plots[[4]] +  xlab("") + ylab("")) |  
      (abu_class_plots[[5]] +  xlab("") + ylab(""))) / 
     (patchwork::wrap_elements(
       full = g_legend(a00 + 
                         guides(fill = guide_legend(
                           override.aes = list(
                             pattern = c("none", "stripe","none","stripe","none","stripe",
                                         rep("none", 4)),
                             fill  = pal_10_q[1:10]),
                           nrow = 2), pattern = "none"))))) + 
  patchwork::plot_layout(heights = c(10,1))

# other hosts 

sup_class_2 <- (((abu_class_plots[[6]] +  xlab("")) |
                   (abu_class_plots[[7]] +  xlab("") + ylab("")) |
                   (abu_class_plots[[8]] +   ylab("")) |  
                   (abu_class_plots[[9]] +  xlab("") + ylab(""))) / 
                  (patchwork::wrap_elements(
                    full = g_legend(a00 + 
                                      guides(fill = guide_legend(
                                        override.aes = list(
                                          pattern = c("none", "stripe","none","stripe","none","stripe",
                                                      rep("none", 4)),
                                          fill  = pal_10_q[1:10]),
                                        nrow = 2), pattern = "none"))))) + 
  patchwork::plot_layout(heights = c(10,1))


# external habitats 

sup_class_3 <- (((abu_class_plots[[10]] +  xlab("")) |
                   (abu_class_plots[[11]] +  xlab("") + ylab("")) |
                   (abu_class_plots[[12]] +   ylab("")) |  
                   (abu_class_plots[[13]] +  xlab("") + ylab(""))) / 
                  (patchwork::wrap_elements(
                    full = g_legend(a00 + 
                                      guides(fill = guide_legend(
                                        override.aes = list(
                                          pattern = c("none", "stripe","none","stripe","none","stripe",
                                                      rep("none", 4)),
                                          fill  = pal_10_q[1:10]),
                                        nrow = 2), pattern = "none"))))) + 
  patchwork::plot_layout(heights = c(10,1))


ggsave("code_R_analysis/output_plots/sup_abundance_class_human.svg", sup_class_1, width = 180, height = 200, unit = "mm")
ggsave("code_R_analysis/output_plots/sup_abundance_class_other_animal.svg", sup_class_2, width = 180, height = 200, unit = "mm")
ggsave("code_R_analysis/output_plots/sup_abundance_class_external_env.svg", sup_class_3, width = 180, height = 200, unit = "mm")





# create the CSC for the supplement, gene classes not in the main results 

d2 <- d1 %>% mutate(texture = unigenes$texture[match(tool_ref, unigenes$tool)])

cs11_sup <- d1 %>% 
  filter(!gene_class %in% top_cso) %>%
  filter(!gene_class %in% c("efflux", "class A", "class B", "class C", "class D")) %>%
  ggplot(aes(x = csc*100, y = 0)) + 
  geom_boxplot_pattern(aes(fill = tool_ref, pattern = texture2),
                       position = position_dodge2(preserve = "single", width = 0.3, padding = 0), 
                       width = 1.3, pattern_color = "black", pattern_fill = "black", pattern_density = 0.000000001,
                       pattern_spacing = 0.2,
                       pattern_size =  0.3, color = "black", outliers = FALSE, outlier.shape = NA, linewidth = 0.15) +
  scale_pattern_manual(values = c('no' = 'none', 'yes' = 'stripe')) +
  geom_text(aes(x = 50, y = 0.8, label = n_obs), size = general_size / .pt) + 
  facet_grid(gene_class ~ tools_labels_ref, scales = "free_y", space = "free",
             labeller = labeller(
               gene_class = as_labeller(function(x) {
                 out <- ifelse(
                   x %in% c("rpoB", "van", "fos", "erm", "cat", "aph", "ant", "bah", "cpa", 
                            "cpt", "dfr", "fai", "mel", "mgt","sat","sul","vph",
                            "aac", "lnu", "nim", "vat", "mph", "qnr"),
                   paste0("italic('", x, "')"),
                   ifelse(x == "abcF", "ABC-F",
                          paste0("'", x, "'"))
                 )
                 return(out)
               }, default = label_parsed)
             )) +
  scale_fill_manual(values = pal_10_q)  + 
  #scale_y_discrete(drop = FALSE) +
  xlab("Class-specific coverage (%)") +
  ylab("ARG class") + 
  theme_minimal() +
  theme5 +
  theme( panel.grid = element_blank(),
         strip.text.x = element_text(size = general_size, vjust = 0, hjust = 0.5 , angle = 0),
         strip.text.y = element_text(size = general_size, vjust = 0, hjust = 0, angle = 0),
         panel.spacing = unit(5, "pt")) +
  scale_y_discrete(labels = function(x) {
    out <- ifelse(
      x %in% c("rpoB", "van", "fos", "erm", "cat", "aph", "ant", "aac", "lnu", "nim","vat","mph","qnr"),
      paste0("italic('", x, "')"),
      ifelse(x %in% "abcF", "ABC-F",
             paste0("'", x, "'"))
    )
    parse(text = out)
  })

ggsave("code_R_analysis/output_plots/fig_csc_sup_class.svg", cs11_sup, width = 180, height = 200, unit = "mm")



### 

# create the three panels of the composition of the human gut core+resistome

core_class_human <- core %>%
  rename(gene_class = new_level) %>% 
  filter(tool %in% basic_tools, habitat %in% "human gut") %>%
  ungroup() %>%
  group_by(centroid) %>% mutate(n_tool = ifelse(n_distinct(tool)==1,"Exclusive to the pipeline", "Found also in other pipelines")) %>% 
  ungroup() %>% 
  group_by(tool, n_tool, gene_class, tools_labels) %>% 
  summarise(n = n()) %>% 
  mutate(r = ifelse(tool %in% basic_tools[c(2,3,4,6)], "1",
                    ifelse(tool %in% basic_tools[c(7,8,9,5)], "2", "3"))) %>% 
  mutate(gene_class = gsub(" beta-lactamase","", gene_class)) %>%
  mutate(gene_class = gsub("rifampin inactivation enzyme","RIF-inact.", gene_class)) %>%
  #mutate(gene_class = gsub("cell wall ","cell wall\n", gene_class)) %>%
  mutate(gene_class = gsub("MFS efflux pump","MFS efflux", gene_class)) %>%
  mutate(gene_class = gsub("efflux pump","efflux", gene_class)) %>%
  mutate(gene_class = gsub("beta-lactam modulation resistance","beta-lactam\nmod.", gene_class)) %>%
  mutate(gene_class = gsub("target-modifying enzyme","target-modif.\nenzyme", gene_class)) %>%
  mutate(gene_class = gsub("self-resistance","self-resistance", gene_class)) 

pcore1 <- ggplot(core_class_human %>% filter(r %in% "1"), aes(x = gene_class, y = n, fill=n_tool)) +
  geom_bar(stat = "identity", width = 0.8) +
  facet_grid(r~tools_labels, scales = "free_x") +
  scale_y_continuous(breaks = scales::pretty_breaks()) +
  scale_x_discrete(labels = function(x) {
    out <- ifelse(
      x %in% c(c("rpoB", "van", "fos", "erm", "cat", "aph", "ant", "aac", "lnu", "nim","vat","mph","qnr"),
               c("rpoB", "van", "fos", "erm", "cat", "aph", "ant", "bah", "cpa", 
                 "cpt", "dfr", "fai", "mel", "mgt","sat","sul","vph",
                 "aac", "lnu", "nim", "vat", "mph", "qnr")),
      paste0("italic('", x, "')"),
      ifelse(x %in% "abcF", "ABC-F",
             paste0("'", x, "'"))
    )
    parse(text = out)
  }) +
  ylab("Number of ARGs") +
  xlab("") +
  scale_fill_manual(values = brewer.pal(8, "Dark2")[c(1,2)]) +#, 
  theme1 +
  theme(
    strip.text.x = element_text(size = general_size, angle = 0, vjust = 0, hjust = 0.5),
    axis.text.x = element_text(size = general_size, angle = 90, vjust = 0.5, hjust = 1),
    strip.text.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    plot.margin = margin(5.5, 5.5, 5.5, 5.5, unit = "pt")) 


pcore2 <- ggplot(core_class_human %>% filter(r %in% "2"), aes(x = gene_class, y = n, fill=n_tool)) +
  geom_bar(stat = "identity", width = 0.8) +
  facet_grid(r~tools_labels, scales = "free_x") +
  scale_y_continuous(breaks = scales::pretty_breaks()) +
  scale_x_discrete(labels = function(x) {
    out <- ifelse(
      x %in% c(c("rpoB", "van", "fos", "erm", "cat", "aph", "ant", "aac", "lnu", "nim","vat","mph","qnr"),
               c("rpoB", "van", "fos", "erm", "cat", "aph", "ant", "bah", "cpa", 
                 "cpt", "dfr", "fai", "mel", "mgt","sat","sul","vph",
                 "aac", "lnu", "nim", "vat", "mph", "qnr")),
      paste0("italic('", x, "')"),
      ifelse(x %in% "abcF", "ABC-F",
             paste0("'", x, "'"))
    )
    parse(text = out)
  }) +
  ylab("Number of ARGs") +
  xlab("") +
  scale_fill_manual(values = brewer.pal(8, "Dark2")[c(1,2)]) +#, 
  theme1 +
  theme(
    axis.text.x = element_text(size = general_size, angle = 90, vjust = 0.5, hjust = 1),
    strip.text.x = element_text(size = general_size, angle = 0, vjust = 0, hjust = 0.5),
    strip.text.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    plot.margin = margin(5.5, 5.5, 5.5, 5.5, unit = "pt")) 


pcore3 <- ggplot(core_class_human %>% filter(r %in% "3"), aes(x = gene_class, y = n, fill=n_tool)) +
  geom_bar(stat = "identity", width = 0.8) +
  facet_grid(r~tools_labels, scales = "free_x") +
  scale_y_continuous(breaks = scales::pretty_breaks()) +
  scale_x_discrete(labels = function(x) {
    out <- ifelse(
      x %in% c(c("rpoB", "van", "fos", "erm", "cat", "aph", "ant", "aac", "lnu", "nim","vat","mph","qnr"),
               c("rpoB", "van", "fos", "erm", "cat", "aph", "ant", "bah", "cpa", 
                 "cpt", "dfr", "fai", "mel", "mgt","sat","sul","vph",
                 "aac", "lnu", "nim", "vat", "mph", "qnr")),
      paste0("italic('", x, "')"),
      ifelse(x %in% "abcF", "ABC-F",
             paste0("'", x, "'"))
    )
    parse(text = out)
  }) +
  ylab("Number of ARGs") +
  xlab("") +
  scale_fill_manual(values = brewer.pal(8, "Dark2")[c(1,2)]) +#, 
  theme1 +
  theme(
    axis.text.x = element_text(size = general_size, angle = 90, vjust = 0.5, hjust = 1),
    strip.text.x = element_text(size = general_size, angle = 0, vjust = 0, hjust = 0.5),
    strip.text.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    plot.margin = margin(5.5, 5.5, 5.5, 5.5, unit = "pt")) 


# merge the ponels

core_human_all <- ((pcore1)  / 
                     (pcore2) /
                     (pcore3 + labs(fill = "") + theme(legend.position = "bottom"))) + patchwork::plot_layout(widths = c(4,4,2))


ggsave("code_R_analysis/output_plots/fig_human_core.svg", core_human_all, width = 180, height = 180, unit = "mm")


###  Supplementary tables 

# core-resistance genes human gut

tools_core <- core %>%
  filter(tool %in% basic_tools, habitat %in% "human gut") %>% 
  group_by(centroid, new_level) %>% 
  summarise(
    tools = paste(sort(unique(tool)), collapse = ", "),
    number_tools = n_distinct(tool),
    .groups = "drop"
  ) %>% 
  arrange(desc(number_tools), tools) %>% 
  rename(unigene = centroid)

gut <- core %>%
  filter(tool %in% basic_tools, habitat == "human gut") %>%
  mutate(tool = as.character(tool)) %>%
  distinct(centroid, tool, new_level)

u <- gut %>%
  group_by(centroid) %>%
  summarise(n_tools  = n_distinct(tool),
            tools    = paste(sort(unique(tool)), collapse = ", "),
            only_RD  = all(tool %in% c("RGI-DIAMOND", "DeepARG")),
            n_levels = n_distinct(new_level),
            level    = paste(sort(unique(new_level)), collapse = " | "),
            .groups = "drop")

nrow(u)
count(u, n_levels)

u %>% filter(n_tools == 1) %>% count(tools)       # only one tool, by tool
sum(u$only_RD)                                   # only RGI and/or DeepARG (any combination)
sum(u$only_RD) / nrow(u)

u %>% filter(only_RD) %>% count(level, sort = TRUE)

u %>% filter(n_tools >= 5) %>% summarise(n = n())
u %>% filter(n_tools >= 5) %>% count(level, sort = TRUE)   # tet RPG 10, erm 5, class A 4?

gut %>% filter(new_level == "tet RPG") %>% count(tool, sort = TRUE)


core %>%
  filter(tool == "RGI-DIAMOND", habitat == "human gut") %>%
  distinct(centroid, new_level) %>%
  count(new_level, sort = TRUE) %>%
  mutate(pct = round(100 * n / sum(n), 1))

pan %>%
  filter(tool == "RGI-DIAMOND", habitat == "human gut") %>%
  group_by(new_level) %>%
  summarise(pan = sum(p_pan), .groups = "drop") %>%
  mutate(pct = round(100 * pan / sum(pan), 1)) %>%
  arrange(desc(pct))

# ARO and gene classes
ARO <- ARO %>% rename(gene_class = new_level, ARO_Term_ID = Term_ID) %>% 
  mutate(abbreviation = gene_class) %>% 
  mutate(abbreviation = gsub(" beta-lactamase","", abbreviation)) %>%
  mutate(abbreviation = gsub("rifampin inactivation enzyme","RIF-inact.", abbreviation)) %>%
  mutate(abbreviation = gsub("MFS efflux pump","efflux", abbreviation)) %>%
  mutate(abbreviation = gsub("efflux pump","efflux", abbreviation)) %>%
  mutate(abbreviation = gsub("beta-lactam modulation resistance","beta-lactam mod.", abbreviation)) %>%
  mutate(abbreviation = gsub("target-modifying enzyme","target-modif. enzyme", abbreviation)) %>%
  mutate(abbreviation = gsub("abcF","ABC-F", abbreviation)) 

write.csv(ARO, file = "code_R_analysis/output_plots/TableS2.csv") 
write.csv(tools_core, file = "code_R_analysis/output_plots/TableS3.csv") 
write.csv(metadata %>% select(sample_id, habitat) %>% rename(sample = sample_id), "code_R_analysis/output_plots/TableS1.csv")



###

disc <- unigenes %>%
  mutate(new_level = ifelse(new_level == "MFS efflux pump", "efflux pump", new_level)) %>% 
  filter(tool %in% basic_tools) %>%
  mutate(tool = as.character(tool)) %>%
  distinct(query, tool, new_level) %>%
  group_by(query) %>%
  filter(n_distinct(new_level) > 1) %>%
  ungroup()

# All pairs of different levels within the same unigene (each pair once)
pairs <- disc %>%
  inner_join(disc, by = "query", suffix = c("_1", "_2"),
             relationship = "many-to-many") %>%
  filter(new_level_1 < new_level_2)

# (a) level vs level: number of unigenes, plus the tools behind each level
ambiguity_aro <- pairs %>%
  group_by(new_level_1, new_level_2) %>%
  summarise(n_unigenes = n_distinct(query),
            tools_1 = paste(sort(unique(tool_1)), collapse = ", "),
            tools_2 = paste(sort(unique(tool_2)), collapse = ", "),
            .groups = "drop") %>%
  arrange(desc(n_unigenes))

ambiguity_aro <- ambiguity_aro %>% rename(class_A=new_level_1,
                         class_B=new_level_2,
                         pipeline_A = tools_1,
                         pipeline_B = tools_2)


n_by_pair <- pairs %>%
  group_by(new_level_1, new_level_2) %>%
  summarise(n_unigenes = n_distinct(query), .groups = "drop")

side_1 <- pairs %>%
  distinct(query, new_level_1, new_level_2, tool_1) %>%
  count(new_level_1, new_level_2, tool_1) %>%
  arrange(new_level_1, new_level_2, tool_1) %>%
  group_by(new_level_1, new_level_2) %>%
  summarise(pipeline_A = paste0(tool_1, " (", n, ")", collapse = ", "),
            .groups = "drop")

side_2 <- pairs %>%
  distinct(query, new_level_1, new_level_2, tool_2) %>%
  count(new_level_1, new_level_2, tool_2) %>%
  arrange(new_level_1, new_level_2, tool_2) %>%
  group_by(new_level_1, new_level_2) %>%
  summarise(pipeline_B = paste0(tool_2, " (", n, ")", collapse = ", "),
            .groups = "drop")

ambiguity_aro <- n_by_pair %>%
  left_join(side_1, by = c("new_level_1", "new_level_2")) %>%
  left_join(side_2, by = c("new_level_1", "new_level_2")) %>%
  rename(class_A = new_level_1, class_B = new_level_2) %>%
  arrange(desc(n_unigenes))

#write.csv(ambiguity_aro, "TableS5.csv", row.names = FALSE)

write.csv(ambiguity_aro, file = "code_R_analysis/output_plots/TableS5.csv", row.names = FALSE) 

# (b) tool vs tool: how often two tools gave different levels for the same unigene
disc %>%
  inner_join(disc, by = "query", suffix = c("_1", "_2"),
             relationship = "many-to-many") %>%
  filter(tool_1 < tool_2, new_level_1 != new_level_2) %>%
  count(tool_1, tool_2, sort = TRUE)

# (c) which tools give each level in the discordant unigenes
disc %>% count(new_level, tool, sort = TRUE)


#

funnel <- unigenes %>%
  mutate(new_level = ifelse(new_level == "MFS efflux pump", "efflux pump", new_level)) %>% 
  filter(tool %in% basic_tools, !is.na(ARO)) %>%
  mutate(tool = as.character(tool)) %>%
  distinct(query, tool, ARO, new_level) %>%
  group_by(query) %>%
  summarise(n_tools = n_distinct(tool),
            n_aro   = n_distinct(ARO),
            n_level = n_distinct(new_level), .groups = "drop") %>%
  filter(n_tools > 1)                       # unigenes detected by >= 2 pipelines

funnel %>% summarise(
  unigenes_in_2plus_pipelines = n(),
  different_aro               = sum(n_aro > 1),
  different_aro_same_level    = sum(n_aro > 1 & n_level == 1),  # removed by hierarchy
  different_new_level         = sum(n_level > 1)                # still discordant
)

##


k_df <- read.csv("code_R_analysis/output_abundance_diversity_resistome/number_samples_per_centroid_per_habitat.csv.gz")


set.seed(202610)
cuts <- 0.5

core_accumulation <- list()
pan_accumulation <- list()

for(j in (1:20)*5){
  df <- k_df %>% 
    mutate(m = j) %>% 
    tidyr::crossing(cut = cuts) %>% 
    mutate(t      = ceiling(cut * m - 1e-9),
           p_core = phyper(t - 1, k, N - k, m, lower.tail = FALSE),
           cnt    = p_core * 500) %>% 
    filter(p_core > 0) %>% 
    select(centroid, new_level_centroid, tool, habitat, cut, p_core, cnt)
  core_accumulation[[j/5]] <- df %>% filter(cnt > 450)   %>% 
    mutate(m = j) %>% 
    filter(tool %in% basic_tools) %>% 
    filter(habitat %in% EN) %>% 
    mutate(habitat = factor(habitat, levels = EN), 
           tool = factor(tool, levels =  tools_levels)) %>% 
    mutate(tool = factor(tool, levels = tools_levels)) %>%
    mutate(tools_labels = factor(tools_labels[tool], levels = tools_labels_factor),
           texture = ifelse(tool %in% tools_texture, "yes", "no"),
           tools_db = factor(tools_db[tool], levels = tools_db_factor))
  
  pan_accumulation[[j/5]] <- k_df %>% 
    mutate(m = j) %>% 
    mutate(p_pan = phyper(0, k, N - k, m, lower.tail = FALSE),   # P(at least 1 of m samples)
           cnt   = p_pan * 500) %>% 
    filter(tool %in% basic_tools) %>% 
    filter(habitat %in% EN) %>% 
    mutate(habitat = factor(habitat, levels = EN), 
           tool = factor(tool, levels =  tools_levels)) %>% 
    mutate(tool = factor(tool, levels = tools_levels)) %>%
    mutate(tools_labels = factor(tools_labels[tool], levels = tools_labels_factor),
           texture = ifelse(tool %in% tools_texture, "yes", "no"),
           tools_db = factor(tools_db[tool], levels = tools_db_factor))
  
}


pan_accumulation_summary <- map_dfr(pan_accumulation, function(df) {
  df %>%
    mutate(metric = "Pan-resistome") %>%
    group_by(metric, tool, tools_labels, texture, tools_db, habitat, new_level_centroid, m) %>%
    summarise(size = sum(p_pan), .groups = "drop")
})


cor_accumulation_summary <- map_dfr(core_accumulation, function(df) {
  df %>% mutate(metric = "Core-resistome") %>% 
    group_by(metric, tool, tools_labels, texture, tools_db, habitat, new_level_centroid, m) %>% 
    summarise(size = n_distinct(centroid), .groups = "drop")
})

pan_acc_plot <- pan_accumulation_summary %>% 
  group_by(tool, tools_labels, texture, tools_db, metric, habitat, m) %>% 
  filter(habitat %in% c("human gut", "pig gut","soil","marine")) %>% 
  summarise(size = sum(size)) %>%
  ggplot(aes(x = m, y = size, fill = tools_labels)) +
  geom_line(aes(linetype = texture), color = "#666666", linewidth = 0.3) +
  facet_wrap(habitat~ .) +
  geom_point(aes(shape = texture), color = "black", stroke = 0.15, size = 2.5, key_glyph = draw_key_hatch) +
  geom_point(data = ~ filter(.x, texture == "yes"), shape = 4, colour = "black",
           size = 1.8, stroke = 0.35, show.legend = FALSE) +
  scale_shape_manual(values = c(21, 24, 22, 23, 25), labels = hab_lab,
                     guide = guide_legend(override.aes = list(fill = "#666666", colour = "black"),
                                          title = ""))  +
  scale_fill_manual(values = pal_tool, guide = "none") +
  xlab("Number of samples") +
  ylab("Pan-resistome") + 
  theme1 

core_acc_plot <- cor_accumulation_summary %>% 
  group_by(tool, tools_labels, texture, tools_db, metric, habitat, m) %>% 
  filter(habitat %in% c("human gut", "pig gut","soil","marine")) %>% 
  summarise(size = sum(size)) %>%
  ggplot(aes(x = m, y = size, fill = tools_labels)) +
  geom_line(aes(linetype = texture), color = "#666666", linewidth = 0.3) +
  facet_wrap(habitat~ .) +
  geom_point(aes(shape = texture), color = "black", stroke = 0.15, size = 2.5, key_glyph = draw_key_hatch) +
  geom_point(data = ~ filter(.x, texture == "yes"), shape = 4, colour = "black",
             size = 1.8, stroke = 0.35, show.legend = FALSE) +
  scale_shape_manual(values = c(21, 24, 22, 23, 25), labels = hab_lab,
                     guide = guide_legend(override.aes = list(fill = "#666666", colour = "black"),
                                          title = ""))  +
  scale_fill_manual(values = pal_tool, guide = "none") +
  xlab("Number of samples") +
  ylab("Core-resistome") + 
  theme1 

pan_acc_plot + core_acc_plot

############################################################
############################################################
############################################################
############################################################
# response to point 3 reviewer 3

# part 1

d_rpob_deeparg <- abundance_class %>% filter(gene == "rpoB", tool == "DeepARG")


d_rpob_deeparg %>%
  group_by(habitat) %>%
  summarise(n = n(), pct_zero = mean(richness == 0) * 100,
            mean_richness = mean(richness), median_richness = median(richness),
            .groups = "drop")

m_rpob_habitat <- glm.nb(richness ~ habitat, data = d_rpob_deeparg)
m_rpob_null    <- glm.nb(richness ~ 1, data = d_rpob_deeparg)

anova(m_rpob_null, m_rpob_habitat, test = "Chisq")   # does habitat explain anything at all?
summary(m_rpob_habitat)

emmeans::emmeans(m_rpob_habitat, ~ habitat, type = "response")

###

reply1 <- abundance_class %>% filter(gene == "rpoB", tool %in% c("DeepARG")) %>% ggplot(aes(x = habitat, y = richness)) +
  geom_boxplot(fill = "#1b9e77", outliers = F) +
  theme1 + theme(legend.position = "none") +
  labs(
    title = "DeepARG - rpoB",
    subtitle = "Distribution of richness per sample across habitats")

reply2 <- abundance_class %>% filter(gene == "van", tool %in% c("RGI-DIAMOND")) %>% ggplot(aes(x = habitat, y = richness)) +
  geom_boxplot(fill = "#e6ab02", outliers = F) +
  theme1 + theme(legend.position = "none") +
  labs(
    title = "RGI - van",
    subtitle = "Distribution of richness per sample across habitats")

####
d_van_deeparg <- abundance_class %>% filter(gene == "van", tool == "RGI-DIAMOND")


d_van_deeparg %>%
  group_by(habitat) %>%
  summarise(n = n(), pct_zero = mean(richness == 0) * 100,
            mean_richness = mean(richness), median_richness = median(richness),
            .groups = "drop")

m_van_habitat <- glm.nb(richness ~ habitat, data = d_van_deeparg)
m_van_null    <- glm.nb(richness ~ 1, data = d_van_deeparg)

anova(m_van_null, m_van_habitat, test = "Chisq")   # does habitat explain anything at all?
summary(m_van_habitat)

emmeans::emmeans(m_van_habitat, ~ habitat, type = "response")


####

reply3 <- abundance_class %>% filter(gene == "class A beta-lactamase", tool %in% basic_tools) %>% ggplot(aes(x = habitat, y = richness)) +
  geom_boxplot(aes(fill = tool), outliers = F) +
  theme1 +
  labs(
    title = "Distribution of class A beta-lactamase richness per sample across habitats") +
  scale_fill_manual(values = pal_10_complete[c(1,1,2,2,3,3,4:7)])

reply4 <- unigenes %>% filter(new_level == "class A beta-lactamase", tool %in% basic_tools) %>% 
  group_by(tool, tools_labels) %>% summarise(n = n_distinct(query)) %>% 
  ggplot(aes(x = tools_labels, y = n)) +
  geom_col(aes(fill = tool), color = "black") +
  theme1 +
  labs(
    title = "Number of unigenes classified as class A beta-lactamase") +
  scale_fill_manual(values = pal_10_complete[c(1,1,2,2,3,3,4:7)]) +
  ylab("Number of unigenes")


d_classA <- abundance_class %>% filter(gene == "class A beta-lactamase", tool %in% basic_tools)
d_classA <- d_classA %>% 
  left_join(unigenes %>% 
              filter(new_level == "class A beta-lactamase") %>% 
              group_by(tool) %>% 
              summarise(n = n_distinct(query)), by = "tool")

m_null    <- glm.nb(richness ~ n, data = d_classA)
m_habitat <- glm.nb(richness ~ habitat + n, data = d_classA)
m_full    <- glm.nb(richness ~ habitat * n, data = d_classA)

anova(m_null, m_habitat, test = "Chisq")    
anova(m_habitat, m_full, test = "Chisq")    

emmeans::emmeans(m_full, ~ n, type = "response")
emmeans::emmeans(m_full, ~ habitat, type = "response")


(reply1 + xlab("")| reply2 + xlab("")| reply4+ xlab("")) / reply3 + xlab("")

#######

# part 2


d_rpob_abund <- abundance_class %>%
  filter(gene == "rpoB", tool == "DeepARG") %>%
  mutate(lrich = log10(richness + 1), labund = log10(abundance + 1))

mean(d_rpob_abund$richness == 0) * 100   # zero check before trusting lm

m_rpob_null    <- lm(labund ~ lrich, data = d_rpob_abund)
m_rpob_habitat <- lm(labund ~ lrich * habitat, data = d_rpob_abund)
anova(m_rpob_null, m_rpob_habitat)       # does the richness-abundance slope vary by habitat?

# --- van, AMRFinderPlus / RGI-DIAMOND / ABRicate-CARD only ---
d_van_abund <- abundance_class %>%
  filter(gene == "van", tool %in% c("AMRFinderPlus", "RGI-DIAMOND", "ABRicate-CARD")) %>%
  mutate(lrich = log10(richness + 1), labund = log10(abundance + 1))

d_van_abund %>% count(tool, habitat) %>% filter(n == 0)   # coverage check first

m_van_null    <- lm(labund ~ lrich, data = d_van_abund)
m_van_habitat <- lm(labund ~ lrich * habitat, data = d_van_abund)
m_van_full    <- lm(labund ~ lrich * habitat + lrich * tool, data = d_van_abund)
anova(m_van_null, m_van_habitat)    # slope varies by habitat?
anova(m_van_habitat, m_van_full)    # slope varies by tool, beyond habitat?

# --- tet RPG, all tools ---
d_tetRPG_abund <- abundance_class %>%
  filter(gene == "tet RPG", tool %in% basic_tools) %>%
  mutate(lrich = log10(richness + 1), labund = log10(abundance + 1))

d_tetRPG_abund %>% count(tool, habitat) %>% filter(n == 0)

m_tetRPG_null    <- lm(labund ~ lrich, data = d_tetRPG_abund)
m_tetRPG_habitat <- lm(labund ~ lrich * habitat, data = d_tetRPG_abund)
m_tetRPG_full    <- lm(labund ~ lrich * habitat + lrich * tool, data = d_tetRPG_abund)
anova(m_tetRPG_null, m_tetRPG_habitat)
anova(m_tetRPG_habitat, m_tetRPG_full)

# --- class A beta-lactamase, all tools ---
d_classA_abund <- abundance_class %>%
  filter(gene == "class A beta-lactamase", tool %in% basic_tools) %>%
  mutate(lrich = log10(richness + 1), labund = log10(abundance + 1))

d_classA_abund %>% count(tool, habitat) %>% filter(n == 0)

m_classA_null    <- lm(labund ~ lrich, data = d_classA_abund)
m_classA_habitat <- lm(labund ~ lrich * habitat, data = d_classA_abund)
m_classA_full    <- lm(labund ~ lrich * habitat + lrich * tool, data = d_classA_abund)
anova(m_classA_null, m_classA_habitat)
anova(m_classA_habitat, m_classA_full)


# --- efflux, all tools ---
d_efflux_abund <- abundance_class %>%
  filter(gene == "efflux pump", tool %in% basic_tools) %>%
  mutate(lrich = log10(richness + 1), labund = log10(abundance + 1))

d_efflux_abund %>% count(tool, habitat) %>% filter(n == 0)

m_efflux_null    <- lm(labund ~ lrich, data = d_efflux_abund)
m_efflux_habitat <- lm(labund ~ lrich * habitat, data = d_efflux_abund)
m_efflux_full    <- lm(labund ~ lrich * habitat + lrich * tool, data = d_efflux_abund)
anova(m_efflux_null, m_efflux_habitat)
anova(m_efflux_habitat, m_efflux_full)

##
emmeans::emtrends(m_tetRPG_full, ~ tool, var = "lrich", at = list(habitat = "human gut"))   # slopes by tool, one habitat
emmeans::emtrends(m_classA_full, ~ tool, var = "lrich", at = list(habitat = "human gut"))   # slopes by tool, one habitat
emmeans::emtrends(m_tetRPG_full, ~ habitat, var = "lrich", at = list(tool = "RGI-DIAMOND"))  # slopes by habitat, one tool
emmeans::emtrends(m_tetRPG_full, ~ habitat, var = "lrich", at = list(tool = "DeepARG"))  # slopes by habitat, one tool
emmeans::emtrends(m_tetRPG_full, ~ habitat, var = "lrich", at = list(tool = "fARGene"))  # slopes by habitat, one tool
emmeans::emtrends(m_tetRPG_full, ~ habitat, var = "lrich", at = list(tool = "AMRFinderPlus"))  # slopes by habitat, one tool
emmeans::emtrends(m_tetRPG_full, ~ habitat, var = "lrich", at = list(tool = "ResFinder"))  # slopes by habitat, one tool

emmeans::emtrends(m_van_full, ~ tool, var = "lrich")
emmeans::emtrends(m_van_full, pairwise ~ tool, var = "lrich")$contrasts

emmeans::emtrends(m_tetRPG_full, ~ tool, var = "lrich")
emmeans::emtrends(m_tetRPG_full, pairwise ~ tool, var = "lrich")$contrasts
emmeans::emtrends(m_tetRPG_full, ~ habitat, var = "lrich")

emmeans::emtrends(m_classA_full, ~ tool, var = "lrich")
emmeans::emtrends(m_classA_full, pairwise ~ tool, var = "lrich")$contrasts
emmeans::emtrends(m_classA_full, ~ habitat, var = "lrich")

emmeans::emtrends(m_rpob_habitat, ~ habitat, var = "lrich")

emmeans::emtrends(m_efflux_full, ~ tool, var = "lrich")
emmeans::emtrends(m_efflux_full, pairwise ~ tool, var = "lrich")$contrasts
emmeans::emtrends(m_efflux_full, ~ habitat, var = "lrich")




# 2. Linear Model for Log-Abundance
d_abundance <- d_sub %>% 
  filter(richness > 0, abundance > 0) %>% 
  mutate(log_abundance = log10(abundance))

m_abundance <- lm(
  log_abundance ~ richness * tool + richness * habitat + richness * gene,
  data = d_abundance
)

# Reduced Model (Drops tool completely)
m_abundance_notool <- lm(
  log_abundance ~ richness + richness * habitat + richness * gene,
  data = d_abundance
)

# Likelihood Ratio Test / ANOVA comparison
anova(m_abundance_notool, m_abundance)

s_rich <- summary(m_richness)$coefficients
res_richness <- data.frame(
  term           = rownames(s_rich),
  estimate       = s_rich[, 1],
  std.error      = s_rich[, 2],
  IRR_FoldChange = exp(s_rich[, 1]),
  IRR_lower      = exp(s_rich[, 1] - 1.96 * s_rich[, 2]),
  IRR_upper      = exp(s_rich[, 1] + 1.96 * s_rich[, 2]),
  p_value        = s_rich[, 4]
)

# Extract Step 2: Linear Abundance Coefficients
s_abun <- summary(m_abundance)$coefficients
res_abundance <- data.frame(
  term      = rownames(s_abun),
  estimate  = s_abun[, 1],
  std.error = s_abun[, 2],
  t_value   = s_abun[, 3],
  p_value   = s_abun[, 4]
)

anova_abun <- Anova(m_abundance, type = "II")
anova_abun$Eta_Sq <- anova_abun$`Sum Sq` / sum(anova_abun$`Sum Sq`)
print(anova_abun)

r2_full   <- summary(m_abundance)$r.squared
r2_notool <- summary(m_abundance_notool)$r.squared

# Variance explained specifically by tool choice
var_explained_by_tool <- (r2_full - r2_notool) * 100
cat("Tool choice explains", round(var_explained_by_tool, 2), "% of abundance variance.\n")


# ==============================================================================
# pan and core
# ==============================================================================

library(car); library(effectsize); library(emmeans)

d <- pan_core_sum %>% 
  ungroup() %>% 
  filter(!is.na(habitat)) %>% 
  mutate(lpan  = log10(Pan + 1),
         lcore = log10(Core + 1))

#  Is there one general rule linking pan and core across all
#     environments?

m_common  <- lm(lcore ~ lpan,           data = d)
m_habitat <- lm(lcore ~ lpan * habitat, data = d)

anova(m_common, m_habitat)

m <- lm(lcore ~ lpan * habitat + lpan * tool, data = d)

# 
#  What matters most for the size of the core-resistome?
# 
q2 <- eta_squared(Anova(m, type = 2), partial = FALSE) %>% 
  as_tibble() %>% 
  transmute(factor = Parameter,
            percent_of_variation = round(100 * Eta2, 1)) %>% 
  arrange(desc(percent_of_variation))

q2
Anova(m, type = 2)






hab_med <- abundance_tool_sample %>%
  filter(tool %in% basic_tools) %>%
  group_by(tool, habitat) %>%
  summarise(med = median(abundance), .groups = "drop") %>%   # use the same median as in your fold table (med10)
  group_by(tool) %>%
  mutate(rank = rank(-med, ties.method = "average")) %>%
  ungroup()

# a) bump chart
ggplot(hab_med, aes(tool, rank, group = habitat, colour = habitat)) +
  geom_line() + geom_point() +
  scale_y_reverse(breaks = 1:13) +
  theme_classic()

rk <- hab_med %>% select(tool, habitat, med) %>% pivot_wider(names_from = tool, values_from = med)
cor(rk[,-1], method = "spearman")



rng <- hab_med %>% group_by(habitat) %>% mutate(lo = min(rank), hi = max(rank), span = hi - lo) %>% ungroup()
ggplot(rng, aes(rank, reorder(habitat, span))) +
  geom_segment(aes(x = lo, xend = hi, yend = reorder(habitat, span)), colour = "grey70") +
  geom_point(aes(colour = tool)) + scale_x_reverse(breaks = 1:13) + theme_classic()


hl <- c("DeepARG", "RGI-DIAMOND")
rng <- hab_med %>% group_by(habitat) %>% mutate(lo = min(rank), hi = max(rank), span = hi - lo) %>% ungroup()
ggplot(rng, aes(y = reorder(habitat, span))) +
  geom_segment(aes(x = lo, xend = hi, yend = reorder(habitat, span)), colour = "grey80") +
  geom_point(data = filter(rng, !tool %in% hl), aes(x = rank), colour = "grey60", size = 2) +
  geom_point(data = filter(rng, tool %in% hl), aes(x = rank, colour = tool), size = 3) +
  scale_x_reverse(breaks = 1:13) +
  labs(x = "Habitat rank (1 = highest median relative abundance)", y = NULL, colour = NULL) +
  theme_classic()



cm <- cor(rk[,-1], method = "spearman")
ord <- hclust(as.dist(1 - cm))$order
lv <- colnames(cm)[ord]

as.data.frame(cm) %>%
  rownames_to_column("a") %>%
  pivot_longer(-a, names_to = "b", values_to = "rho") %>%
  mutate(a = factor(a, levels = lv), b = factor(b, levels = rev(lv))) %>%
  ggplot(aes(a, b, fill = rho)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.2f", rho)), size = 3) +
  scale_fill_gradient(low = "#f2f0f7", high = "#54278f", limits = c(0.6, 1), name = "Spearman ρ") +
  labs(x = NULL, y = NULL) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))