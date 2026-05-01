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

################################################################################
################################################################################
## DATASETS 

# gene classes 

ARO <- read.csv("code_R_analysis/output_abundance_diversity_resistome/conversion_ARO_parent_new_level.csv")

# results per individual tool 

lst <- readRDS("code_R_analysis/output_abundance_diversity_resistome/results_tools.rds")

# metagenomes' metadata

metadata <- read.delim("data/metadata_GMGC10.sample.meta.tsv")
metadata <- metadata %>% filter(!habitat %in% c("amplicon", "isolate", "built-environment"))

# metadata0 will be merged with abundance 
metadata0 <- metadata %>% select(sample_id, habitat) %>% rename(sample = sample_id)

abundance <- readRDS("code_R_analysis/output_abundance_diversity_resistome/abundance_diversity.rds")

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

core <- readRDS(file = "code_R_analysis/output_abundance_diversity_resistome/core_resistome.rds")

# convert MFS to efflux pump
core <- core %>% mutate(new_level_centroid = ifelse(new_level_centroid == "MFS efflux pump", "efflux pump", new_level_centroid))

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

# Extract the core-resistome with condition that a gene most be in
# 450 sub-core resistomes, and those are composed of genes present in at least 
# 50% of the metagenomic samples in a subsample

sumcore <- sum_core_adjust(core, 450, 0.5) %>%
  mutate(tools_labels = factor(tools_labels[tool], levels = tools_labels_factor),
         texture = ifelse(tool %in% tools_texture, "yes", "no"),
         tools_db = factor(tools_db[tool], levels = tools_db_factor))

## PAN RESISTOME

pan <- readRDS(file = "code_R_analysis/output_abundance_diversity_resistome/pan_resistome.rds")

# convert MFS to efflux pump
pan <- pan %>% mutate(gene_class = ifelse(gene_class == "MFS efflux pump", "efflux pump", gene_class))

pan <- pan %>% 
  mutate(habitat = factor(habitat, levels = EN), 
         tool = factor(tool, levels = tools_levels)) %>% 
  mutate(tools_labels = factor(tools_labels[tool], levels = tools_labels_factor),
         texture = ifelse(tool %in% tools_texture, "yes", "no"),
         tools_db = factor(tools_db[tool], levels = tools_db_factor))

# calculate the total number of unigenes in the pan-resistome of each sample
sumpan2 <- pan %>% ungroup() %>% 
  group_by(tool, habitat, epoch) %>% 
  summarise(s = sum(unigenes)) %>%
  ungroup() %>% 
  group_by(tool, habitat) %>% 
  summarise(md = median(s), mn = mean(s), sd = sd(s)) %>%
  mutate(tools_labels = factor(tools_labels[tool], levels = tools_labels_factor),
         texture = ifelse(tool %in% tools_texture, "yes", "no"),
         tools_db = factor(tools_db[tool], levels = tools_db_factor))

# merge the core and pan resistomes at the sample level
pan_core <- sumpan2 %>% 
  left_join((sumcore %>% 
               ungroup() %>% 
               group_by(tool, habitat) %>% 
               summarise(core = sum(unigenes))), by = c("tool", "habitat")) %>%
  mutate(core = ifelse(is.na(core), 0, core)) %>% 
  mutate(prop = core / md) %>%
  ungroup() %>% group_by(tool, habitat) %>% 
  mutate(texture = ifelse(tool %in% tools_texture, "yes", "no")) %>%
  mutate(tools_labels = factor(tools_labels[tool], levels = tools_labels_factor),
         texture = ifelse(tool %in% tools_texture, "yes", "no"),
         tools_db = factor(tools_db[tool], levels = tools_db_factor))


## unigenes identified by tool

unigenes <- readRDS(file = "code_R_analysis/output_abundance_diversity_resistome/unigenes_per_tool.rds") %>% 
  # convert MFS to efflux pump
  mutate(new_level = ifelse(new_level == "MFS efflux pump", "efflux pump", new_level)) %>%  
  mutate(tool = factor(tool, levels = tools_levels)) %>%
  mutate(tools_labels = factor(tools_labels[tool], levels = tools_labels_factor),
         texture = ifelse(tool %in% tools_texture, "yes", "no"),
         tools_db = factor(tools_db[tool], levels = tools_db_factor))


# Overlap between tools by gene class
# per class and tool
recall_fnr <- create_class_overlaps(unigenes)

# overlap by pipeline (without considering classes)
JI_all <- return_overlap_tools(unigenes)


# sort the gene classes by abundance and richness to decide which classes to
# show in sup_abundance plots
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


################################################################################################
# PLot 1 

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
                                            "ABRicate/\nResFinder\n(n = 37,772)", "AMRFinder-\nPlus (aa)\n(n = 3,141)"))) %>% 
  mutate(data_number = factor(data_type, levels = c( "nucleotide", "amino acid"))) %>% 
  ggplot(aes(x = id , color = tool_3, fill = tool_3)) +
  stat_ecdf(geom = "step", linewidth = 1) + 
  scale_color_manual(values = pal_7[c(1,5,3, 6)], 
                     guide = guide_legend(nrow = 1)) +
  ggtitle("") + 
  ylab("Empirical Cumulative \nDistribution Function") + 
  xlab("% Identity") +
  labs(color = "", fill = "") +
  scale_y_continuous(limits = c(0,1), expand = c(0.01, 0.01)) + 
  scale_x_continuous(limits = c(0,100), expand = c(0, 0.5)) + 
  theme_minimal() + 
  theme1 +
  theme(legend.position = "bottom", panel.border = element_blank())


# prepare the data for the proportion of gene class in each tool

df_plot <- unigenes %>% 
  filter(tool %in% basic_tools) %>% 
  group_by(tool, tools_labels, tools_db, new_level) %>% 
  summarise(n = n()) %>% 
  mutate(p = n / sum(n)) %>% 
  ungroup() %>% 
  group_by(new_level) %>%
  filter(max(p, na.rm = TRUE) >= 0.03) %>%  
  ungroup() %>%
  # change the name of gene classes to fit in the plot
  mutate(new_level = gsub(" beta-lactamase","", new_level)) %>%
  mutate(new_level = gsub("rifampin inactivation enzyme","RIF-inact. enz.", new_level)) %>%
  mutate(new_level = gsub("MFS efflux pump","MFS efflux", new_level)) %>%
  mutate(new_level = gsub("efflux pump","efflux", new_level)) %>%
  mutate(new_level = gsub("beta-lactam modulation resistance","beta-lactam\nmod.", new_level)) %>%
  mutate(new_level = gsub("target-modifying enzyme","target-modif.\nenzyme", new_level)) %>%
  mutate(new_level = gsub("self-resistance","self-resistance", new_level)) 

# plot the heatmap of the proportion of each gene class by tool
df_plot1 <-  df_plot %>% 
  filter(tool %in% basic_tools) %>% 
  ggplot(aes(x = tools_labels, y = new_level, fill = p)) + 
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
        plot.margin = margin(0, 0, 0, 20, unit = "pt"),
        strip.text.x = element_text(size = general_size, angle = 0, vjust = 0, hjust = 0.5))


# prepare the data for the Jaccard index plot 

mat <- JI_all %>% 
  filter(tool_ref %in% basic_tools, tool_comp %in% basic_tools) %>% 
  select(tool_ref, tool_comp, jaccard) %>%
  pivot_wider(names_from = tool_comp, values_from = jaccard)

mat_matrix <- as.matrix(mat[,-1])   # remove tool_ref column
rownames(mat_matrix) <- mat$tool_ref
dist_mat <- as.dist(1 - mat_matrix)

# hierarchical clustering to order the matrix
hc <- hclust(dist_mat)  # or "ward.D2"
plot(hc)
ordered_tools <- hc$labels[hc$order]

# filter the tools and rearrange the data

JI_all_plot <- JI_all %>% 
  filter(tool_ref %in% basic_tools, tool_comp %in% basic_tools) %>% 
  mutate(tool_lab_ref = factor(tools_labels[tool_ref], levels = tools_labels[ordered_tools]),
         tool_lab_comp = factor(tools_labels[tool_comp], levels = tools_labels[ordered_tools]))


# plot the jaccard index 

jaccard_plot <-  JI_all_plot %>% filter(as.numeric(tool_lab_ref) < as.numeric(tool_lab_comp)) %>% 
  filter(tool_ref %in% basic_tools, tool_comp %in% basic_tools) %>% 
  ggplot(aes(x = tool_lab_ref, y = tool_lab_comp, fill = jaccard)) + 
  geom_tile(color = "grey") + 
  scale_fill_gradientn( colors = brewer.pal(9, "YlOrBr"),
                        labels = percent_format(accuracy = 1),
                        breaks = c(0, 0.2, 0.4, 0.6, 0.8)) + 
  geom_text(
    data = JI_all_plot %>% 
      filter(tool_ref %in% basic_tools, tool_comp %in% basic_tools, 
             as.numeric(tool_lab_ref) < as.numeric(tool_lab_comp), 
             jaccard >= 0.05),
    aes(label = scales::percent(jaccard, accuracy = 1),
        color = jaccard >= 0.50), 
    size = general_size / .pt, show.legend = F) +
  scale_color_manual(values = c( "black", "white")) +
  scale_y_discrete(labels = function(x) {
    x <- gsub("ABRicate-\n", "ABRicate-", x)
    x <- gsub("-\nPlus", "Plus", x)
    x}, expand = 0) +
  labs(fill = "") + 
  ylab("") + 
  xlab("") + 
  theme1 +
  theme(legend.position = "bottom", panel.grid = element_blank(),
        panel.border =  element_blank(),
        plot.margin = margin(0, 0, 0, 10, unit = "pt"),
        strip.text.x = element_text(size = general_size, angle = 0, vjust = 0, hjust = 0.5))

jaccard_plot

# merge all panels of figure 1 together 


p2_1 <- ((p2a + ggtitle("a")) /
           (jaccard_plot + ggtitle("c")) /
           (id_plot + ggtitle("d")) +
           patchwork::plot_layout(heights = c(1,1,1)) |
           (df_plot1 + ggtitle("b"))) + 
  patchwork::plot_layout(widths = c(1,1))


ggsave("code_R_analysis/output_plots/fig1.svg", p2_1, width = 180, height = 180, unit = "mm")


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
                         width = 1.3, pattern_color = "black", pattern_fill = "black", pattern_density = 0.000000001,
                         pattern_spacing = 0.2,
                         pattern_size =  0.3, color = "black", outliers = FALSE, outlier.shape = NA, linewidth = 0.15) +
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
          plot.margin = margin(5.5, 5.5, 5.5, 5.5, unit = "pt")) 
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

a0 <- df_a0 %>% 
  filter(tool %in% basic_tools) %>% 
  filter(habitat %in% c("human gut")) %>% 
  ungroup() %>% group_by(habitat) %>%  
  mutate(N_samples = length(unique(sample))) %>% 
  ungroup() %>%
  mutate(habitat_label = paste0(habitat, "\n(n = ", scales::comma(N_samples[1]), ")")) %>%
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
                          paste0("'", x, "'"))
                 )
                 return(out)
               }, default = label_parsed)
             )) + 
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
    plot.margin = margin(0, 0, 0, 5.5, unit = "pt")) +
  guides(fill = guide_legend(
    override.aes = list(
      pattern = c(rep("none", 5),"stripe","none","stripe","none","stripe"),
      fill  = pal_10_q[1:10])), pattern = "none")

#a0

# merge the figure 2 panel a with the appropriate legend
left_block <- (left1/patchwork::wrap_elements(full = g_legend(a0))) + 
  patchwork::plot_layout(heights = c(15,1))

# merge all panels in figure 2 

p30 <- (left_block | (a0 + theme(legend.position = "none"))) + patchwork::plot_layout(widths = c(4,1))

ggsave("code_R_analysis/output_plots/fig2.svg", p30, width = 180, height = 130, unit = "mm")


# create the supplementary figure sup_abundance.svg abundance per class and habitat

sup_abundance_left <- 
  (abu_plots[[2]] + theme(axis.text.x = element_blank(), strip.text.x = element_blank()) + xlab("") + ylab("")) / 
  (abu_plots[[3]] + theme(axis.text.x = element_blank(), strip.text.x = element_blank()) + xlab("") + ylab("Relative abundance\n(aligned reads per million)")) /
  (abu_plots[[4]] + theme(axis.text.x = element_blank(), strip.text.x = element_blank()) + xlab("") + ylab("")) /
  (abu_plots[[5]] + theme(strip.text.x = element_blank()) + xlab("") + ylab(""))

sup_abundance_right <- 
  (abu_plots[[6]] + theme(axis.text.x = element_blank(), strip.text.x = element_blank()) + xlab("") + ylab("")) / 
  (abu_plots[[7]] + theme(axis.text.x = element_blank(), strip.text.x = element_blank()) + xlab("") + ylab("")) /
  (abu_plots[[8]] + theme(strip.text.x = element_blank()) + xlab("") + ylab("")) / 
  patchwork::wrap_elements(full = g_legend(a0 + 
                                             guides(fill = guide_legend(
                                               override.aes = list(
                                                 pattern = c(rep("none", 5),"stripe","none","stripe","none","stripe"),
                                                 fill  = pal_10_q[1:10]),
                                               nrow = 4), pattern = "none")))

sup_abundance <- sup_abundance_left | sup_abundance_right

ggsave("code_R_analysis/output_plots/sup_abundance.svg", sup_abundance, width = 180, height = 180, unit = "mm")

################################################################################################

# create vector for shaping the geom_points in the core-resistome and pan-resistome plots

shape_tools <- rep(21, length(tools_labels[tools_levels %in% basic_tools]))
shape_tools[tools_levels[tools_levels %in% basic_tools] %in% tools_texture] <- 24

# for the pan and core, define only one metric name for both and add the type of data (core or pan)

pan_core_df_plot <- pan_core %>% 
  select(!c(md, sd)) %>% 
  pivot_longer(cols = c(mn, core), names_to = "metric", values_to = "value") %>% 
  mutate(metric = ifelse(metric %in% "mn", "Pan-resistome", metric)) %>%
  mutate(metric = ifelse(metric %in% "core", "Core-resistome", metric)) %>%
  mutate(metric = factor(metric, levels = c("Pan-resistome", "Core-resistome"))) %>%
  filter(value > 0) 

# plot the pan resistome

p4a1 <- pan_core_df_plot %>% filter(metric %in% "Pan-resistome") %>% 
  filter(habitat %in% c("human gut", "pig gut", "wastewater", "marine", "freshwater", "soil")) %>% 
  filter(tool %in% basic_tools) %>%
  ggplot(aes(y = fct_rev(tool), x =  value)) +
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
  geom_segment(aes(x = 0, xend = value, y = tool, yend = tool, color = tool), linewidth = 0.2, show.legend = F) +
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

p4a2 <- pan_core_df_plot %>% filter(!metric %in% "Pan-resistome") %>% 
  filter(habitat %in% c("human gut", "pig gut", "wastewater", "marine", "freshwater", "soil")) %>% 
  filter(tool %in% basic_tools) %>%
  ggplot(aes(y = fct_rev(tool), x =  value)) +
  geom_point(aes(fill = tool, shape = texture),  color = "black", stroke = 0.3, size = 1.7) + 
  facet_grid(habitat ~ metric, scales = "free") +
  scale_fill_manual(values = pal_10_q) +
  scale_shape_manual(values = c("no" = 21, "yes" = 24)) +
  geom_segment(aes(x = 0, xend = value, y = tool, yend = tool, color = tool), linewidth = 0.2, show.legend = F) +
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

ggsave("code_R_analysis/output_plots/fig3.svg", p4, width = 120, height = 140, unit = "mm")

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
  mutate(val = recall, d = "csc") %>% 
  mutate(tools_labels_comp = factor(tools_labels[tool_comp], levels = tools_labels_factor),
         texture_comp = ifelse(tool_comp %in% tools_texture, "yes", "no"),
         tools_db_comp = factor(tools_db[tool_comp], levels = tools_db_factor)) %>% 
  mutate(tools_labels_comp = factor(gsub("-\n","",tools_labels_comp), levels = gsub("-\n","", levels(tools_labels_comp)))) %>%
  mutate(tools_labels_ref = factor(tools_labels[tool_ref], levels = tools_labels_factor))


# change the name of some classes for plotting

d1 <- d1 %>% 
  filter(tool_ref %in% basic_tools, tool_comp %in% basic_tools) %>%
  mutate(new_level = gsub(" beta-lactamase","", new_level)) %>%
  mutate(new_level = gsub("MFS efflux pump","MFS efflux", new_level)) %>%
  mutate(new_level = gsub("efflux pump","efflux", new_level)) %>% 
  group_by(tool_ref,tools_labels_ref, new_level) %>% mutate(n_obs = n()) %>% 
  mutate(n_obs = paste0('n = ',n_obs))


# plot the csc of top_cso classes

cs11 <- d1 %>% 
  filter(new_level %in% top_cso) %>% 
  ggplot(aes(x = recall*100, y = 0)) + 
  geom_boxplot_pattern(aes(fill = tool_ref, pattern = texture2),
                       position = position_dodge2(preserve = "single", width = 0.3, padding = 0), 
                       width = 1.3, pattern_color = "black", pattern_fill = "black", pattern_density = 0.000000001,
                       pattern_spacing = 0.2,
                       pattern_size =  0.3, color = "black", outliers = FALSE, outlier.shape = NA, linewidth = 0.15) +
  scale_pattern_manual(values = c('no' = 'none', 'yes' = 'stripe')) +
  geom_text(aes(x = 50, y = 0.9, label = n_obs), size = general_size / .pt) + 
  facet_grid(new_level ~ tools_labels_ref, scales = "free_y", space = "free",
             labeller = labeller(
               new_level = as_labeller(function(x) {
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


ggsave("code_R_analysis/output_plots/fig4.svg", cs11, width = 180, height = 100, unit = "mm")


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

rgi_query <- data.frame(gene = unique(c(lst$rgi.blast$query, lst$rgi.diamond$query, lst$rgi.diamond.prot$query)))
rgi_query <- rgi_query %>% 
  mutate(dataset = 
           ifelse(gene %in% intersect(intersect(lst$rgi.blast$query, lst$rgi.diamond$query), lst$rgi.diamond.prot$query), "all setups",
                  ifelse(gene %in% intersect(lst$rgi.blast$query, lst$rgi.diamond$query), "BLAST nt and DIAMOND (nt or aa)",
                         ifelse(gene %in% intersect(lst$rgi.diamond.prot$query, lst$rgi.diamond$query), "only\nDIAMOND\nnt and/or aa",
                                ifelse(gene %in% intersect(lst$rgi.diamond.prot$query, lst$rgi.blast$query), "BLAST nt and DIAMOND (nt or aa)",
                                       ifelse(gene %in% lst$rgi.diamond.prot$query, "only\nDIAMOND\nnt and/or aa",
                                              ifelse(gene %in% lst$rgi.diamond$query, "only\nDIAMOND\nnt and/or aa",
                                                     ifelse(gene %in% lst$rgi.blast$query, "only BLAST\nnt", NA)))))))) %>%
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

deeparg_over <- plot_overlaps(lst$deeparg.norm$query, lst$deeparg.norm.prot$query, "nt and aa", "only nt", "only aa", nam = "DeepARG")

# overlap for fargene

fargene_over <- plot_overlaps(lst$fargene$query, lst$fargene.prot$query, "nt and aa", "only nt", "only aa", nam = "fARGene")

# overlap for amrfinder

amrfinder_over <- plot_overlaps(lst$amrfinder.norm$query, lst$amrfinder.norm.prot$query, "nt and aa", "only nt", "only aa", nam = "AMRFinder-\nPlus")


#  merge rgi, fargene, deeparg, amrfinder

plot_overlaps_same_tool <- 
  ((rgi_over + xlab("") + ggtitle("a")) / 
     (deeparg_over + theme(legend.position = "none") + ggtitle("c")) /
     (fargene_over + theme(legend.position = "none")) /
     (amrfinder_over + xlab("ARGs"))) + 
  patchwork::plot_layout(heights = c(1, 1, 1, 1)) & 
  theme(plot.margin = margin(5, 5, 5, 5, unit = "pt"))


#  plot card, ncbi, resfinder

card_over <- plot_overlaps(lst$rgi.diamond.prot$query, lst$abricate.card.norm$query, "RGI and\nABRicate", "only RGI", "only ABRicate", nam = "CARD")
ncbi_over <- plot_overlaps(lst$amrfinder.norm.prot$query, lst$abricate.ncbi.norm$query, "AMRFinderPlus\nand ABRicate", "only AMRFinderPlus", "only ABRicate", nam = "NCBI")
resfinder_over <- plot_overlaps(lst$resfinder.norm$query, lst$abricate.resfinder.norm$query, "ResFinder\nand ABRicate", "only ResFinder", "only ABRicate", nam = "ResFinder")

#  merge

plot_db <- 
  (card_over + ggtitle("b")) / 
  ncbi_over / 
  resfinder_over + xlab("ARGs") & 
  theme(plot.margin = margin(5, 5, 5, 5, unit = "pt"))


# plot for megares vs other abricate 

other_abricate <- unique(c(lst$abricate.argannot.norm$query, lst$abricate.card.norm$query, 
                           lst$abricate.ncbi.norm$query, lst$abricate.resfinder.norm$query))

abricate_query <- data.frame(rbind(lst$abricate.megares.norm[,c("tool", "query")], 
                                   lst$abricate.argannot.norm[,c("tool", "query")], 
                                   lst$abricate.card.norm[,c("tool", "query")], 
                                   lst$abricate.ncbi.norm[,c("tool", "query")], 
                                   lst$abricate.resfinder.norm[,c("tool", "query")]))

ntools_abricate <- abricate_query %>% 
  group_by(query) %>% 
  summarise(n_tools = n_distinct(tool))

# modify the labels for intercepts with megares 

ntools_abricate <- ntools_abricate %>% mutate(s = ifelse(query %in% lst$abricate.megares.norm$query, "MEGARes","other"))
ntools_abricate <- ntools_abricate %>% mutate(s = ifelse(n_tools == 2 & s == "MEGARes" & query %in% lst$abricate.argannot.norm$query, "MEGARes and\nARGANNOT", s))
ntools_abricate <- ntools_abricate %>% mutate(s = ifelse(n_tools == 2 & s == "MEGARes" & query %in% lst$abricate.ncbi.norm$query, "MEGARes and\nNCBI", s))
ntools_abricate <- ntools_abricate %>% mutate(s = ifelse(n_tools == 2 & s == "MEGARes" & query %in% lst$abricate.resfinder.norm$query, "MEGARes and\nResFinder", s))
ntools_abricate <- ntools_abricate %>% mutate(s = ifelse(n_tools == 2 & s == "MEGARes" & query %in% lst$abricate.card.norm$query, "MEGARes and\nCARD", s))
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


# plot the supplementary figures of the pan and core resistome
# pan
pan_core_supp <- pan_core_df_plot %>% filter(metric %in% "Pan-resistome") %>% 
  filter(!habitat %in% c("human gut", "pig gut", "wastewater", "marine", "freshwater", "soil")) %>% 
  filter(tool %in% basic_tools) %>%
  ggplot(aes(y = fct_rev(tool), x =  value)) +
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
  geom_segment(aes(x = 0, xend = value, y = tool, yend = tool, color = tool), linewidth = 0.2, show.legend = F) +
  scale_color_manual(values = pal_10_q) +
  guides(fill = guide_legend(
    override.aes = list(
      shape = shape_tools,
      color  = pal_10_q[1:10])), shape = "none") + 
  theme(axis.text.y = element_blank(), 
        legend.position = "bottom",
        strip.text.x = element_text(size = general_size, vjust = 0, hjust = 0.5),
        strip.text.y = element_text(size = general_size, vjust = 0.5, hjust = 0.5, angle = 0),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.border = element_blank(),
        plot.margin = margin(5, 0, 0, 0, unit = "pt")) +
  scale_x_reverse(labels = label_comma()) 

# core
pan_core_supp2 <- pan_core_df_plot %>% filter(!metric %in% "Pan-resistome") %>% 
  filter(!habitat %in% c("human gut", "pig gut", "wastewater", "marine", "freshwater", "soil")) %>% 
  filter(tool %in% basic_tools) %>%
  ggplot(aes(y = fct_rev(tool), x =  value)) +
  geom_point(aes(fill = tool, shape = texture),  color = "black", stroke = 0.3, size = 1.7) + 
  facet_grid(habitat ~ metric, scales = "free") +
  scale_fill_manual(values = pal_10_q) +
  scale_shape_manual(values = c("no" = 21, "yes" = 24)) +
  geom_segment(aes(x = 0, xend = value, y = tool, yend = tool, color = tool), linewidth = 0.2, show.legend = F) +
  scale_color_manual(values = pal_10_q) +
  theme_minimal() +
  ylab("") +
  xlab("Number of ARGs") + 
  ggtitle("b") + 
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
  theme1 +
  theme(axis.text.y = element_blank(), 
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.text.y = element_blank(),
        panel.border = element_blank(),
        plot.margin = margin(5, 0, 0, 0, unit = "pt")) +
  scale_x_continuous(labels = label_comma()) 

# merge 
pan_core_supp | pan_core_supp2

# merge labels 

panel1_s <- (pan_core_supp + theme(legend.position = "none") | pan_core_supp2) + patchwork::plot_layout(widths = c(1,1))
p4_S <- (panel1_s / (g_legend(p4a1) )) + patchwork::plot_layout(heights = c(12,1))

ggsave("code_R_analysis/output_plots/sup_pan_core.svg", p4_S, width = 120, height = 180, unit = "mm")


####


# create the plots of richness

rich_plots <- list()
set.seed(2026)

for(j in 1:length(EN)){
  
  df_plot <- abundance_tool_sample %>%
    filter(tool %in% basic_tools) %>% 
    filter(habitat %in% EN[j]) %>% 
    group_by(habitat) %>% 
    mutate(N = n_distinct(sample),
           richness = richness)
  
  df_jitter <- df_plot %>% 
    ungroup() %>%
    group_by(habitat, tool) %>% 
    filter(richness < quantile(richness, 0.75) + 1.5*IQR(richness))
  
  lims_richness_max = lims_richness %>%
    filter(habitat %in% EN[j]) %>%
    filter(tool %in% basic_tools) %>% 
    ungroup() %>% summarise(max(w2)) %>% pull()
  
  rich_plots[[j]] <- lims_richness %>%
    mutate(habitat_label = paste0(habitat, "\n(n = ", scales::comma(df_plot$N[1]), ")")) %>%
    filter(habitat %in% EN[j]) %>%
    filter(tool %in% basic_tools) %>% 
    ggplot(aes(x = tools_labels, y=median, fill = tools_db))  + 
    geom_boxplot_pattern(aes(ymin = w1, lower = q25, middle = median, upper = q75, ymax = w2, pattern = texture), stat = "identity",
                         position = position_dodge2(preserve = "single", width = 0.3, padding = 0), 
                         width = 1.3, pattern_color = "black", pattern_fill = "black", pattern_density = 0.000000001,
                         pattern_spacing = 0.2,
                         pattern_size =  0.3, color = "black", outliers = FALSE, outlier.shape = NA, linewidth = 0.15) +
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
  
  all_top_abundance <- unique(
    levels_abundance_div %>% 
      filter(habitat %in% EN[j]) %>% 
      select(gene) %>% distinct() %>% slice_head(n=10) %>%
      ungroup() %>% pull())
  
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
        pattern = c(rep("none", 5),"stripe","none","stripe","none","stripe"),
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
       full = g_legend(a0 + 
                         guides(fill = guide_legend(
                           override.aes = list(
                             pattern = c(rep("none", 5),"stripe","none","stripe","none","stripe"),
                             fill  = pal_10_q[1:10]),
                           nrow = 2), pattern = "none"))))) + 
  patchwork::plot_layout(heights = c(10,1))

# other hosts 

sup_class_2 <- (((abu_class_plots[[6]] +  xlab("")) |
                   (abu_class_plots[[7]] +  xlab("") + ylab("")) |
                   (abu_class_plots[[8]] +   ylab("")) |  
                   (abu_class_plots[[9]] +  xlab("") + ylab(""))) / 
                  (patchwork::wrap_elements(
                    full = g_legend(a0 + 
                                      guides(fill = guide_legend(
                                        override.aes = list(
                                          pattern = c(rep("none", 5),"stripe","none","stripe","none","stripe"),
                                          fill  = pal_10_q[1:10]),
                                        nrow = 2), pattern = "none"))))) + 
  patchwork::plot_layout(heights = c(10,1))


# external habitats 

sup_class_3 <- (((abu_class_plots[[10]] +  xlab("")) |
                   (abu_class_plots[[11]] +  xlab("") + ylab("")) |
                   (abu_class_plots[[12]] +   ylab("")) |  
                   (abu_class_plots[[13]] +  xlab("") + ylab(""))) / 
                  (patchwork::wrap_elements(
                    full = g_legend(a0 + 
                                      guides(fill = guide_legend(
                                        override.aes = list(
                                          pattern = c(rep("none", 5),"stripe","none","stripe","none","stripe"),
                                          fill  = pal_10_q[1:10]),
                                        nrow = 2), pattern = "none"))))) + 
  patchwork::plot_layout(heights = c(10,1))


ggsave("code_R_analysis/output_plots/sup_abundance_class_human.svg", sup_class_1, width = 180, height = 200, unit = "mm")
ggsave("code_R_analysis/output_plots/sup_abundance_class_other_animal.svg", sup_class_2, width = 180, height = 200, unit = "mm")
ggsave("code_R_analysis/output_plots/sup_abundance_class_external_env.svg", sup_class_3, width = 180, height = 200, unit = "mm")



# create the CSC for the supplement, gene classes not in the main results 

d2 <- d1 %>% mutate(texture = unigenes$texture[match(tool_ref, unigenes$tool)])

cs11_sup <- d1 %>% 
  filter(!new_level %in% top_cso) %>%
  filter(!new_level %in% c("efflux", "class A", "class B", "class C", "class D")) %>%
  ggplot(aes(x = recall*100, y = 0)) + 
  geom_boxplot_pattern(aes(fill = tool_ref, pattern = texture2),
                       position = position_dodge2(preserve = "single", width = 0.3, padding = 0), 
                       width = 1.3, pattern_color = "black", pattern_fill = "black", pattern_density = 0.000000001,
                       pattern_spacing = 0.2,
                       pattern_size =  0.3, color = "black", outliers = FALSE, outlier.shape = NA, linewidth = 0.15) +
  scale_pattern_manual(values = c('no' = 'none', 'yes' = 'stripe')) +
  geom_text(aes(x = 50, y = 0.8, label = n_obs), size = general_size / .pt) + 
  facet_grid(new_level ~ tools_labels_ref, scales = "free_y", space = "free",
             labeller = labeller(
               new_level = as_labeller(function(x) {
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

core_class_human <- core %>% filter(cut == 0.5, cnt > 450, tool %in% basic_tools, habitat %in% "human gut") %>%
  ungroup() %>%
  group_by(X) %>% mutate(n_tool = ifelse(n_distinct(tool)==1,"Exclusive to the pipeline", "Found also in other pipelines")) %>% 
  ungroup() %>% 
  group_by(tool, n_tool, new_level) %>% 
  summarise(n = n()) %>% 
  mutate(tools_labels = factor(tools_labels[tool], levels =tools_labels_factor)) %>% 
  mutate(r = ifelse(tool %in% basic_tools[c(2,3,4,6)], "1",
                    ifelse(tool %in% basic_tools[c(7,8,9,10)], "2", "3"))) %>% 
  mutate(new_level = gsub(" beta-lactamase","", new_level)) %>%
  mutate(new_level = gsub("rifampin inactivation enzyme","RIF-inact. enz.", new_level)) %>%
  #mutate(new_level = gsub("cell wall ","cell wall\n", new_level)) %>%
  mutate(new_level = gsub("MFS efflux pump","MFS efflux", new_level)) %>%
  mutate(new_level = gsub("efflux pump","efflux", new_level)) %>%
  mutate(new_level = gsub("beta-lactam modulation resistance","beta-lactam\nmod.", new_level)) %>%
  mutate(new_level = gsub("target-modifying enzyme","target-modif.\nenzyme", new_level)) %>%
  mutate(new_level = gsub("self-resistance","self-resistance", new_level)) 

pcore1 <- ggplot(core_class_human %>% filter(r %in% "1"), aes(x = new_level, y = n, fill=n_tool)) +
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


pcore2 <- ggplot(core_class_human %>% filter(r %in% "2"), aes(x = new_level, y = n, fill=n_tool)) +
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


pcore3 <- ggplot(core_class_human %>% filter(r %in% "3"), aes(x = new_level, y = n, fill=n_tool)) +
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


tools_core <- core %>%
  ungroup() %>% 
  filter(cut %in% 0.5 & cnt > 450) %>%
  filter(tool %in% basic_tools, habitat %in% "human gut") %>% 
  group_by(X, new_level) %>% 
  summarise(
    tools = paste(sort(unique(tool)), collapse = ", "),
    number_tools = n_distinct(tool),
    .groups = "drop"
  ) %>% 
  arrange(desc(number_tools), tools) %>% 
  rename(unigene = X, gene_class = new_level)


ARO <- ARO %>% rename(gene_class = new_level, ARO_Term_ID = Term_ID) %>% 
  mutate(abbreviation = gene_class) %>% 
  mutate(abbreviation = gsub(" beta-lactamase","", abbreviation)) %>%
  mutate(abbreviation = gsub("rifampin inactivation enzyme","RIF-inact. enz.", abbreviation)) %>%
  mutate(abbreviation = gsub("MFS efflux pump","efflux", abbreviation)) %>%
  mutate(abbreviation = gsub("efflux pump","efflux", abbreviation)) %>%
  mutate(abbreviation = gsub("beta-lactam modulation resistance","beta-lactam mod.", abbreviation)) %>%
  mutate(abbreviation = gsub("target-modifying enzyme","target-modif. enzyme", abbreviation)) %>%
  mutate(abbreviation = gsub("abcF","ABC-F", abbreviation)) 

write.csv(ARO, file = "code_R_analysis/output_plots/TableS2.csv") 
write.csv(tools_core, file = "code_R_analysis/output_plots/TableS3.csv") 
write.csv(metadata %>% select(sample_id, habitat) %>% rename(sample = sample_id), "code_R_analysis/output_plots/TableS1.csv")






