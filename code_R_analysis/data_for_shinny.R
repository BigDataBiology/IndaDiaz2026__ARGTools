library(dplyr)
library(ggplot2)
library(tidyverse)
library(RColorBrewer)


#setwd("~/Documents/GitHub/arg_compare/code_R_analysis") 
options(dplyr.summarise.inform = FALSE)
# source("helper.R")
source("code_R_analysis/helper.R")

# Sourced gene classes 

general_size <- 6
lab_fn <- function(x) {
  x <- gsub("-", "-\n", x)
  x <- gsub(" ", "\n", x)
  x <- gsub("/", "/\n", x)
  x
}

# FORMAT PLOTS OLD COLORS
pal_7 <- brewer.pal(8, "Dark2")
pal_7 <- pal_7[-7]
pal_7 <- pal_7[c(1,2,3,4,6,5,7)]
#pal_7 <- pal_7[c(1,7,2,6,3,5,4)]
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

basic_tools <- c(
  "DeepARG", "fARGene","ABRicate-ARGANNOT", "ABRicate-MEGARes",
  "RGI-DIAMOND", "ABRicate-CARD","AMRFinderPlus", "ABRicate-NCBI",
  "ResFinder", "ABRicate-ResFinder")

tools_levels <- c(
  "DeepARG", "DeepARG70","DeepARG80","DeepARG90",
  "fARGene","ABRicate-ARGANNOT", "ABRicate-MEGARes",
  "RGI-DIAMOND", "RGI-DIAMOND70","RGI-DIAMOND80","RGI-DIAMOND90",
  "ABRicate-CARD","AMRFinderPlus", "ABRicate-NCBI",
  "ResFinder", "ABRicate-ResFinder",
  "DeepARG-aa", "RGI-BLAST", "RGI-DIAMOND-aa", "fARGene-aa", "AMRFinderPlus-nt")

names(pal_10_q) <- basic_tools

da <- rep(pal_10_q["DeepARG"],3)
names(da) <- c("DeepARG70","DeepARG80","DeepARG90")
rgi <- rep(pal_10_q["RGI-DIAMOND"],3)
names(rgi) <- c("RGI-DIAMOND70","RGI-DIAMOND80","RGI-DIAMOND90")
other <- c(pal_10_q["DeepARG"], pal_10_q["RGI-DIAMOND"],
           pal_10_q["RGI-DIAMOND"],pal_10_q["fARGene"],pal_10_q["AMRFinderPlus"])
names(other) <- c("DeepARG-aa", "RGI-BLAST", 
                  "RGI-DIAMOND-aa", "fARGene-aa","AMRFinderPlus-nt")
pal_10_q <- c(pal_10_q[1],da,pal_10_q[2:5],rgi,pal_10_q[6:10],other)
#pal_10_q <- c(pal_10_q, da, rgi, other)
rm(da, rgi, other)

tools_labels <- c(
  "DeepARG", "DeepARG-70%","DeepARG-80%","DeepARG-90%", "fARGene", "ABRicate", "ABRicate",
  "RGI", "RGI-70%","RGI-80%","RGI-90%",
  "ABRicate", "AMRFinder-\nPlus", "ABRicate",
  "ResFinder", "ABRicate",
  "DeepARG-aa", "RGI/nBLAST", "RGI-aa", "fARGene-aa", "AMRFinder-\nPlus-nt")

names(tools_labels) <- tools_levels

tools_labels_factor <- c(
  "DeepARG","DeepARG-70%","DeepARG-80%","DeepARG-90%", "fARGene", "RGI",
  "RGI-70%","RGI-80%","RGI-90%", "AMRFinder-\nPlus", 
  "ResFinder", "ABRicate",
  "DeepARG-aa", "RGI/nBLAST", "RGI-aa", "fARGene-aa", "AMRFinder-\nPlus-nt")

# one space for DeepARG
# two spaces for fARGene

tools_db <- c(rep("",4), " ", "ARG-\nANNOT", "MEGA-\nRes", rep("CARD", 4),"CARD",
              "NCBI","NCBI", 
              "Res-\nFinder","Res-\nFinder","","CARD","CARD",
              " ","NCBI")

# one space for DeepARG
# two spaces for fARGene
tools_db_factor <- c("", " ", "ARG-\nANNOT", "MEGA-\nRes",
                     "CARD", "NCBI", "Res-\nFinder")

tools_texture <- c("ABRicate-CARD", "ABRicate-NCBI", "ABRicate-ResFinder")

## DATASETS 

# results per individual tool 
lst <- readRDS("code_R_analysis/output_abundance_diversity_resistome/results_tools.rds")

# metagenomes' metadata
metadata <- read.delim("data/metadata_GMGC10.sample.meta.tsv")

# abundance <- readRDS("output_abundance_diversity_resistome/abundance_diversity.rds")
abundance <- readRDS("code_R_analysis/output_abundance_diversity_resistome/abundance_diversity.rds")



abundance <- abundance %>%
  mutate(habitat = metadata$habitat[match(sample, metadata$sample)]) %>% 
  mutate(habitat = factor(habitat, levels = EN),
         tool = factor(tool, levels = tools_levels)) %>%
  mutate(tools_labels = factor(tools_labels[tool], levels = tools_labels_factor),
         texture = ifelse(tool %in% tools_texture, "yes", "no"),
         tools_db = factor(tools_db[tool], levels = tools_db_factor))


abundance_tool_sample <- abundance %>% 
  ungroup() %>%
  group_by(tool, sample, habitat) %>%  
  summarise(abundance = sum(abundance), richness = sum(richness)) %>% 
  ungroup() %>% 
  complete(sample, tool) %>% # complete with NAs
  left_join(abundance %>% select(sample, habitat) %>% 
              distinct(), by = "sample") %>% # get habitat and habitat2 
  mutate(habitat  = coalesce(habitat.x, habitat.y)) %>%
  select(-habitat.x, -habitat.y) %>% 
  mutate(abundance = replace_na(abundance, 0)) %>%  # change NAs to 0
  mutate(richness = replace_na(richness, 0)) %>% # change NAs to 0
  arrange(tool, sample) %>%
  mutate(tools_labels = factor(tools_labels[tool], levels = tools_labels_factor),
         texture = ifelse(tool %in% tools_texture, "yes", "no"),
         tools_db = factor(tools_db[tool], levels = tools_db_factor)) %>%
  group_by(tool, habitat) %>%
  mutate(N_samples = n_distinct(sample)) %>%
  ungroup()


abundance_class <- abundance %>% 
  ungroup() %>% 
  select(-c(richness_no_rarified)) %>%
  mutate(gene = factor(gene), sample = factor(sample), 
         tool = factor(tool, levels = tools_levels)) %>%
  complete( sample, gene, tool, 
            fill = list(abundance = 0,  richness = 0)) %>% ## complete for all samples, tools, and gene classes
  mutate(gene = as.character(gene),
         habitat = metadata$habitat[match(sample, metadata$sample_id)])  %>%
  mutate(tools_labels = factor(tools_labels[tool], levels = tools_labels_factor),
         texture = ifelse(tool %in% tools_texture, "yes", "no"),
         tools_db = factor(tools_db[tool], levels = tools_db_factor))

## CORE RESISTOME 

core <- readRDS(file = "code_R_analysis/output_abundance_diversity_resistome/core_resistome.rds")

core <- core %>% 
  rename(new_level = new_level_centroid, 
         X = centroid) %>% 
  filter(tool %in% tools_levels) %>% 
  mutate(habitat = factor(habitat, levels = EN), 
         tool = factor(tool, levels =  tools_levels)) %>% 
  mutate(tool = factor(tool, levels = tools_levels)) %>%
  mutate(tools_labels = factor(tools_labels[tool], levels = tools_labels_factor),
         texture = ifelse(tool %in% tools_texture, "yes", "no"),
         tools_db = factor(tools_db[tool], levels = tools_db_factor))

sumcore <- sum_core_adjust(core, 450, 0.5) %>%
  mutate(tools_labels = factor(tools_labels[tool], levels = tools_labels_factor),
         texture = ifelse(tool %in% tools_texture, "yes", "no"),
         tools_db = factor(tools_db[tool], levels = tools_db_factor),
         cut = 0.5, cnt = 450)

for(cnts in c(250,300,350,400,450)){
  for(cuts in c(0.3,0.4,0.5,0.6,0.7,0.8,0.9)){
    if(cnts == 450 & cuts == 0.5) {
      next} else{
        sumcore <- sumcore %>% bind_rows(
          sum_core_adjust(core, 450, 0.5) %>%
            mutate(tools_labels = factor(tools_labels[tool], levels = tools_labels_factor),
                   texture = ifelse(tool %in% tools_texture, "yes", "no"),
                   tools_db = factor(tools_db[tool], levels = tools_db_factor),
                   cut = 0.5, cnt = 450)) 
      }
  }
}



## PAN RESISTOME

pan <- readRDS(file = "code_R_analysis/output_abundance_diversity_resistome/pan_resistome.rds")
pan <- pan %>% 
  mutate(habitat = factor(habitat, levels = EN), 
         tool = factor(tool, levels = tools_levels)) %>% 
  mutate(tools_labels = factor(tools_labels[tool], levels = tools_labels_factor),
         texture = ifelse(tool %in% tools_texture, "yes", "no"),
         tools_db = factor(tools_db[tool], levels = tools_db_factor))

sumpan2 <- pan %>% ungroup() %>% 
  group_by(tool, habitat, epoch) %>% 
  summarise(s = sum(unigenes)) %>%
  ungroup() %>% 
  group_by(tool, habitat) %>% 
  summarise(md = median(s), mn = mean(s), sd = sd(s)) %>%
  mutate(tools_labels = factor(tools_labels[tool], levels = tools_labels_factor),
         texture = ifelse(tool %in% tools_texture, "yes", "no"),
         tools_db = factor(tools_db[tool], levels = tools_db_factor))




## unigenes identified by tool

unigenes <- readRDS(file = "code_R_analysis/output_abundance_diversity_resistome/unigenes_per_tool.rds") %>% 
  mutate(tool = factor(tool, levels = tools_levels)) %>%
  mutate(tools_labels = factor(tools_labels[tool], levels = tools_labels_factor),
         texture = ifelse(tool %in% tools_texture, "yes", "no"),
         tools_db = factor(tools_db[tool], levels = tools_db_factor))


# Jaccard index, recall/class specific concordance, fnr/class specific non-overlap
# per class and tool
recall_fnr <- create_class_overlaps(unigenes) 


lims_abundance <- abundance_tool_sample %>%
  group_by(tool, habitat,tools_db, tools_labels) %>% 
  summarise(
    median = ifelse(quantile(abundance, 0.5) < 0, 0, quantile(abundance, 0.5)),
    q25 = ifelse(quantile(abundance, 0.25) < 0, 0, quantile(abundance, 0.25)),
    q75 = ifelse(quantile(abundance, 0.75) < 0, 0, quantile(abundance, 0.75)),
    w1 = ifelse(quantile(abundance, 0.25) - 1.5*IQR(abundance)<0,0,
                quantile(abundance, 0.25) - 1.5*IQR(abundance)),
    w2 = ifelse(quantile(abundance, 0.75) + 1.5*IQR(abundance)<0,0,
                quantile(abundance, 0.75) + 1.5*IQR(abundance))) 

lims_richness <- abundance_tool_sample %>%
  group_by(tool, habitat,tools_db, tools_labels) %>% 
  summarise(
    median = ifelse(quantile(richness, 0.5) < 0, 0, quantile(richness, 0.5)),
    q25 = ifelse(quantile(richness, 0.25) < 0, 0, quantile(richness, 0.25)),
    q75 = ifelse(quantile(richness, 0.75) < 0, 0, quantile(richness, 0.75)),
    w1 = ifelse(quantile(richness, 0.25) - 1.5*IQR(richness)<0,0,
                quantile(richness, 0.25) - 1.5*IQR(richness)),
    w2 = ifelse(quantile(richness, 0.75) + 1.5*IQR(richness)<0,0,
                quantile(richness, 0.75) + 1.5*IQR(richness))) 


abundance_class_summary <- abundance_class %>%
  ungroup() %>%
  group_by(tool, habitat, gene,texture, tools_labels, tools_db) %>% 
  summarise(
    q50 = ifelse(quantile(abundance, 0.5) < 0, 0, quantile(abundance, 0.5)),
    q25 = ifelse(quantile(abundance, 0.25) < 0, 0, quantile(abundance, 0.25)),
    q75 = ifelse(quantile(abundance, 0.75) < 0, 0, quantile(abundance, 0.75)),
    w1 = ifelse(quantile(abundance, 0.25) - 1.5*IQR(abundance)<0,0,
                quantile(abundance, 0.25) - 1.5*IQR(abundance)),
    w2 = ifelse(quantile(abundance, 0.75) + 1.5*IQR(abundance)<0,0,
                quantile(abundance, 0.75) + 1.5*IQR(abundance))) 



lst_results <- list(abundance_tool_sample=abundance_tool_sample, 
     core=core, sumpan2=sumpan2, unigenes=unigenes, 
     recall_fnr=recall_fnr, abundance_class_summary=abundance_class_summary,
     sumcore = sumcore)

saveRDS(lst_results, file = "shiny/app/data.rds", compress = T)






