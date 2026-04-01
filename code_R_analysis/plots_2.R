library(dplyr)
library(ggplot2)
#library(plotly)
library(ggVennDiagram)
library(gridExtra)
library(tidyverse)
library(RColorBrewer)
library(ggpattern)
library(grid)
library(eulerr)
library(reactable)
library(Cairo)
library(ggalluvial)
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

basic_tools <- c(
  "DeepARG", "fARGene","ABRicate-ARGANNOT", "ABRicate-MEGARes",
  "RGI-DIAMOND", "ABRicate-CARD","AMRFinderPlus", "ABRicate-NCBI",
  "ResFinder", "ABRicate-ResFinder")

tools_levels <- c(
  "DeepARG", "fARGene","ABRicate-ARGANNOT", "ABRicate-MEGARes",
  "RGI-DIAMOND", "ABRicate-CARD","AMRFinderPlus", "ABRicate-NCBI",
  "ResFinder", "ABRicate-ResFinder","DeepARG70","DeepARG80","DeepARG90",
  "RGI-DIAMOND70","RGI-DIAMOND80","RGI-DIAMOND90",
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
  
pal_10_q <- c(pal_10_q, da, rgi, other)
rm(da, rgi, other)

tools_labels <- c(
  "DeepARG", "fARGene", "ABRicate", "ABRicate",
  "RGI", "ABRicate", "AMRFinder-\nPlus", "ABRicate",
  "ResFinder", "ABRicate","DeepARG-70%","DeepARG-80%","DeepARG-90%","RGI-70%","RGI-80%","RGI-90%",
  "DeepARG-aa", "RGI/nBLAST", "RGI-aa", "fARGene-aa", "AMRFinder-\nPlus-nt")

names(tools_labels) <- tools_levels

tools_labels_factor <- c(
  "DeepARG", "fARGene", "RGI","AMRFinder-\nPlus", 
  "ResFinder", "ABRicate","DeepARG-70%","DeepARG-80%","DeepARG-90%",
  "RGI-70%","RGI-80%","RGI-90%",
  "DeepARG-aa", "RGI/nBLAST", "RGI-aa", "fARGene-aa", "AMRFinder-\nPlus-nt")

# one space for DeepARG
# two spaces for fARGene

tools_db <- c(" ", "  ", "ARG-\nANNOT", "MEGA-\nRes", "CARD", "CARD","NCBI","NCBI", 
              "Res-\nFinder","Res-\nFinder"," "," "," ","CARD","CARD","CARD",
              " ", "CARD", "CARD", "  ", "NCBI")

# one space for DeepARG
# two spaces for fARGene
tools_db_factor <- c(" ", "  ", "ARG-\nANNOT", "MEGA-\nRes",
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
         tools_db = factor(tools_db[tool], levels = tools_db_factor))


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
         tools_db = factor(tools_db[tool], levels = tools_db_factor))

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
  mutate(tool = factor(tool, levels = tools_levels)) %>%
  mutate(tools_labels = factor(tools_labels[tool], levels = tools_labels_factor),
         texture = ifelse(tool %in% tools_texture, "yes", "no"),
         tools_db = factor(tools_db[tool], levels = tools_db_factor))


# Jaccard index, recall/class specific concordance, fnr/class specific non-overlap
# per class and tool
recall_fnr <- create_class_overlaps(unigenes)


# per tool
JI_all <- return_overlap_tools(unigenes)



levels_unigenes <- unigenes %>% 
  ungroup() %>% 
  group_by(query) %>% 
  slice_head(n = 1) %>% 
  ungroup() %>% 
  group_by(new_level) %>% 
  summarise( n = n()) %>% 
  arrange(desc(n)) %>%
  ungroup() %>% 
  pull(new_level)


#"cell wall charge", "rpoB",
levels_abundance_div <- abundance_class %>% 
  group_by(habitat, tool, gene) %>% 
  summarise(a = sum(abundance), d = sum(richness)) %>% 
  arrange(desc(a), desc(d)) %>% 
  ungroup()


# levels_abundance_div %>% filter(habitat %in% "human gut") %>% 
#  select(gene) %>% distinct() %>% print(n = 15)

# levels_abundance_div %>% filter(habitat %in% "pig gut") %>% 
#   select(gene) %>% distinct() %>% print(n = 15)



top_abundance <- c("efflux pump", "van" , "class A beta-lactamase", 
                   "tet RPG",  "cell wall charge", "MFS efflux pump", 
                   "rpoB", "erm")


#top20

top_cso <- c("van", "efflux pump",  "tet RPG", "class A beta-lactamase", 
             "class B beta-lactamase","class C beta-lactamase", "class D beta-lactamase",
             "aph", "MFS efflux pump", "erm", "aac")


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
tools_texture_code <- rep("none", length(tools_labels))
tools_texture_code[tools_levels %in% tools_texture] <- "stripe"

tools_main <- c("DeepARG", "fARGene", "ABRicate-MEGARes",
                "RGI-DIAMOND", "AMRFinderPlus", "ResFinder")

type_tools <- c("solid", "solid", "solid", "dotted")

g1 <- guides(fill = guide_legend(
    override.aes = list(
      pattern = tools_texture,
      fill  = pal_7)), pattern = "none")

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
        strip.text.x = element_text(size = general_size, angle = 90, vjust = 0.5, hjust = 0))

p2a

id_plot <- unigenes %>% 
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
  ggplot(aes(x = id , color = tool_2, fill = tool_2)) +
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

id_plot


fig2ab <- plot_grid(
  p2a + ggtitle("a"), id_plot + ggtitle("c"),
  ncol = 1, align = "hv", axis = "tblr"
)


df_plot <- unigenes %>% 
  filter(tool %in% basic_tools) %>% 
  group_by(tool, tools_labels, tools_db, new_level) %>% 
  summarise(n = n()) %>% 
  mutate(p = n / sum(n)) %>% 
  ungroup() %>% 
  group_by(new_level) %>%
  filter(max(p, na.rm = TRUE) >= 0.03) %>%  
  ungroup() %>%
  mutate(new_level = gsub(" beta-lactamase","", new_level)) %>%
  mutate(new_level = gsub("rifampin inactivation enzyme","RIF-inact. enz.", new_level)) %>%
  #mutate(new_level = gsub("cell wall ","cell wall\n", new_level)) %>%
  mutate(new_level = gsub("MFS efflux pump","MFS efflux", new_level)) %>%
  mutate(new_level = gsub("efflux pump","efflux", new_level)) %>%
  mutate(new_level = gsub("beta-lactam modulation resistance","beta-lactam\nmod.", new_level)) %>%
  mutate(new_level = gsub("target-modifying enzyme","target-\nmodif. enz.", new_level)) %>%
  mutate(new_level = gsub("self-resistance","self-resistance", new_level)) 


df_plot1 <-  df_plot %>% 
  filter(tool %in% basic_tools) %>% 
  ggplot(aes(x = tools_labels, y = new_level, fill = p)) + 
  geom_tile(color = "grey") + 
  scale_fill_gradientn( colors = brewer.pal(9, "YlOrBr"),
                        labels = percent_format(accuracy = 1),
                        breaks = c(0.0001, 0.2, 0.4, 0.6)) + 
  facet_grid(.~tools_db, scales = "free_x", space = "free") +
  geom_text(
    data = df_plot %>% filter(p >= 0.03),  # only these get labels
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
        plot.margin = margin(0, 0, 0, 10, unit = "pt"),
        strip.text.x = element_text(size = general_size, angle = 90, vjust = 0.5, hjust = 0))

df_plot1



p2_1 <- ((p2a + ggtitle("a")) /
  (id_plot + ggtitle("c")) +
  patchwork::plot_layout(heights = c(1,1)) |
  (df_plot1 + ggtitle("b"))) + 
  patchwork::plot_layout(widths = c(1,1))


ggsave("code_R_analysis/output_plots/fig1.svg", p2_1, width = 180, height = 110, unit = "mm")


# tools 



env_to_plot_main <- EN[c(1, 9, 10, 11, 12, 13)]
env_to_plot_supp <- EN[!EN %in% env_to_plot_main]

scale_f <- 10^0

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

  
abu_plots <- list()
set.seed(2026)

for(j in 1:length(EN)){
  
  df_plot <- abundance_tool_sample %>%
    filter(tool %in% basic_tools) %>% 
    filter(habitat %in% EN[j]) %>% 
    group_by(habitat) %>% 
    mutate(N = n_distinct(sample),
           abundance = abundance)
  
  abu_plots[[j]] <- df_plot %>%
    filter(tool %in% basic_tools) %>% 
    ggplot(aes(x = tools_labels, y = abundance, fill = tools_db, pattern = texture)) + 
    geom_boxplot_pattern(position = position_dodge2(preserve = "single", width = 0.3, padding = 0), 
                         width = 1.3, pattern_color = "black", pattern_fill = "black", pattern_density = 0.000000001,
                         pattern_spacing = 0.2,
                         pattern_size =  0.3, color = "black", outliers = FALSE, outlier.shape = NA,
                         linewidth = 0.15) +
    scale_x_discrete(expand = expansion(add = 1)) +#, 
                     #breaks = levels(df_plot$tool)[levels(df_plot$tool) %in% df_plot$tool],
                     #labels = gsub("RGI-DIAMOND", "RGI", 
                     #levels(df_plot$tool)[levels(df_plot$tool) %in% df_plot$tool])) +
    #geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
    geom_jitter(data = df_plot %>% 
                  ungroup() %>%
                  group_by(habitat, tool) %>% 
                  filter(abundance < quantile(abundance, 0.75) + 1.5*IQR(abundance)) %>%
                  slice_sample(n = min(100, df_plot$N[1])), width = 0.35, size = 0.4, alpha = 0.1) + 
    facet_grid(habitat ~ tools_db  , scales = "free_x", space = "free_x") +
    scale_y_continuous(labels = scales::comma) + 
    scale_fill_manual(values = pal_7) +
    scale_pattern_manual(values = c('no' = 'none', 'yes' = 'stripe')) +
    theme1 + 
    theme(panel.border = element_blank(),
          strip.text.x = element_text(size = general_size, angle = 90, vjust = 0.5, hjust = 0),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          plot.margin = margin(5.5, 5.5, 5.5, 0, unit = "pt")) 
}

main_abundance_left <- 
  (abu_plots[[1]] + theme(axis.text.x = element_blank(), strip.text.x = element_text(size = general_size, angle = 90, vjust = 0.5, hjust = 0)) + xlab("") + ylab("") + ggtitle("a")) / 
  (abu_plots[[10]] + theme(axis.text.x = element_blank(), strip.text.x = element_blank()) + xlab("") +  ylab("Relative abundance")) /
  (abu_plots[[11]] + theme(strip.text.x = element_blank()) + xlab("") + ylab("")) 
  
main_abundance_right <- 
  (abu_plots[[9]] + theme(axis.text.x = element_blank(), strip.text.x = element_text(size = general_size, angle = 90, vjust = 0.5, hjust = 0)) + xlab("") + ylab("")) / 
  (abu_plots[[12]] + theme(axis.text.x = element_blank(), strip.text.x = element_blank()) + xlab("") +  ylab("")) /
  (abu_plots[[13]] + theme(strip.text.x = element_blank()) + xlab("") + ylab("")) 

left1 <- main_abundance_left | main_abundance_right


gene_levels <- levels(df_abundance_class_human$gene)
gene_levels <- gsub(" beta-lactamase","", gene_levels)
gene_levels <- gsub("cell wall ","cell\nwall\n", gene_levels)
gene_levels <- gsub("MFS efflux pump","MFS\nefflux", gene_levels)
gene_levels <- gsub("efflux pump","efflux", gene_levels)


df_abundance_class_all <- abundance_class %>% 
  mutate(gene = ifelse(gene %in% top_abundance, gene, "other")) %>% 
  mutate(texture = ifelse(tool %in% tools_texture, "yes", "no")) %>%
  mutate(gene = factor(gene, levels = c(top_abundance, "other"))) %>%
  ungroup() %>% 
  group_by(habitat, sample, tool, gene) %>%
  summarise(abundance = sum(abundance, richness = sum(richness)), texture = texture[1]) %>%
  mutate(tools_labels = factor(tools_labels[tool], levels = tools_labels[!duplicated(tools_labels)]),
         texture = ifelse(tool %in% tools_texture, "yes", "no"),
         tools_db = factor(tools_db[tool], levels = tools_db[!duplicated(tools_db)]))

df_a0 <- df_abundance_class_all  %>% 
  mutate(gene = gsub(" beta-lactamase","", gene)) %>%
  mutate(gene = gsub("cell wall ","cell\nwall\n", gene)) %>%
  mutate(gene = gsub("MFS efflux pump","MFS\nefflux", gene)) %>%
  mutate(gene = gsub("efflux pump","efflux", gene)) %>%
  mutate(gene = factor(gene, levels = gene_levels)) 

pal_10_q_2 <- pal_10_q
names(pal_10_q_2) <- tools_levels



a0 <- df_a0 %>% 
  filter(tool %in% basic_tools) %>% 
  filter(habitat %in% c("human gut")) %>% 
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
  facet_grid(gene ~ habitat, scales = "free", space = "free") +
  xlab("Relative abundance") +
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

a0

left_block <- (left1/patchwork::wrap_elements(full = g_legend(a0))) + patchwork::plot_layout(heights = c(15,1))
p30 <- (left_block | (a0 + theme(legend.position = "none"))) + patchwork::plot_layout(widths = c(3,1))
ggsave("code_R_analysis/output_plots/fig2.svg", p30, width = 180, height = 130, unit = "mm")


#a0_new <- (a0_van / (a0_not_van + theme(legend.position = "none"))) + patchwork::plot_layout(heights = c(1,9))




sup_abundance_left <- 
  (abu_plots[[2]] + theme(axis.text.x = element_blank()) + xlab("") + ylab("")) / 
  (abu_plots[[3]] + theme(axis.text.x = element_blank(), strip.text.x = element_blank()) + xlab("") + ylab("Relative abundance")) /
  (abu_plots[[4]] + theme(axis.text.x = element_blank(), strip.text.x = element_blank()) + xlab("") + ylab("")) /
  (abu_plots[[5]] + theme(strip.text.x = element_blank()) + xlab("") + ylab(""))

sup_abundance_right <- 
  (abu_plots[[6]] + theme(axis.text.x = element_blank()) + xlab("") + ylab("")) / 
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

###################
###################
###################

shape_tools <- rep(21, length(tools_labels[tools_levels %in% basic_tools]))
shape_tools[tools_levels[tools_levels %in% basic_tools] %in% tools_texture] <- 24


pan_core_df_plot <- pan_core %>% 
  select(!c(md, sd)) %>% 
  pivot_longer(cols = c(mn, core), names_to = "metric", values_to = "value") %>% 
  mutate(metric = ifelse(metric %in% "mn", "Pan-resistome", metric)) %>%
  mutate(metric = ifelse(metric %in% "core", "Core-resistome", metric)) %>%
  mutate(metric = factor(metric, levels = c("Pan-resistome", "Core-resistome"))) %>%
  filter(value > 0) 


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
        legend.position = "bottom",
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.border = element_blank(),
        plot.margin = margin(5, 0, 0, 0, unit = "pt")) +
  scale_x_reverse(labels = label_comma()) 

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
        strip.text.x = element_text(size = general_size, vjust = 0, hjust = 0.5),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.text.y = element_blank(),
        panel.border = element_blank(),
        plot.margin = margin(5, 0, 0, 0, unit = "pt")) +
  scale_x_continuous(labels = label_comma()) 

p4a1 | p4a2

df_sumcore <- sumcore %>% filter(habitat %in% "human gut") %>% 
  group_by(tool) %>% 
  mutate(p = unigenes / sum(unigenes)) %>% 
  mutate(new_level = gsub(" beta-lactamase","", new_level)) %>%
  mutate(new_level = gsub("rifampin inactivation enzyme","RIF-inact. enz.", new_level)) %>%
  mutate(new_level = gsub("inactivation enzyme","\ninact. enz.", new_level)) %>%
  mutate(new_level = gsub("cell wall ","cell wall\n", new_level)) %>%
  mutate(new_level = gsub("MFS efflux pump","MFS efflux", new_level)) %>%
  mutate(new_level = gsub("efflux pump","efflux", new_level)) %>%
  mutate(new_level = gsub("beta-lactam modulation resistance","beta-lactam\nmod.", new_level)) %>%
  mutate(new_level = gsub("target-modifying enzyme","target-\nmodif. enz.", new_level)) %>%
  mutate(new_level = gsub("self-resistance","self-resistance", new_level)) %>% 
  mutate(new_level = gsub("variant or","variant or\n", new_level)) 


df_sumcore_levels <- df_sumcore %>% ungroup() %>% 
  filter(tool %in% basic_tools) %>%
  arrange(desc(p)) %>% 
  group_by(new_level) %>%
  filter(p == max(p)) %>% 
  pull(new_level) 



panel1 <- (p4a1 + theme(legend.position = "none") | p4a2) + patchwork::plot_layout(widths = c(1,1))
p4 <- (panel1 / (g_legend(p4a1) )) + patchwork::plot_layout(heights = c(12,1))

ggsave("code_R_analysis/output_plots/fig3.svg", p4, width = 120, height = 140, unit = "mm")

#########
#########
#########




theme5 <- theme(
  legend.position = "none",
  legend.text = element_text(size = general_size ),
  #panel.border = element_blank(),
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




cs11 <- ggplot(d1 %>% 
                 filter(tool_ref %in% basic_tools, tool_comp %in% basic_tools) %>%
                 filter(tool_ref %in% c("DeepARG","fARGene", "RGI-DIAMOND")) %>% 
                 filter(new_level %in% top_cso) %>%
                 mutate(new_level = gsub(" beta-lactamase","", new_level)) %>%
                 mutate(new_level = gsub("MFS efflux pump","MFS efflux", new_level)) %>%
                 mutate(new_level = gsub("efflux pump","efflux", new_level)), 
               aes(x = recall*100, y = new_level)) + 
  geom_boxplot(aes(fill = tool_ref),
               position = position_dodge2(preserve = "single"),
               color = "black", outliers = FALSE) +
  facet_grid(tool_ref ~ " ", scales = "free_y") +
  scale_fill_manual(values = pal_10_q[match(c("DeepARG","fARGene", "RGI-DIAMOND","AMRFinderPlus"), tools_levels)]) +
  scale_y_discrete(drop = FALSE) +
  xlab("%") +
  ylab("ARG class") + 
  theme_minimal() +
  theme5 +
  theme( panel.grid = element_blank(),
         strip.text.x = element_text(size = general_size, vjust = 0, hjust = 0.5),
         strip.text.y = element_text(size = general_size, vjust = 0, hjust = 0.5))

cs11


cs2 <- ggplot(d1 %>% 
                filter(tool_ref %in% basic_tools, tool_comp %in% basic_tools) %>%
                filter(tool_ref %in% c("RGI-DIAMOND", "fARGene","DeepARG")), #%>%
              aes(x = recall*100, y = tools_labels_comp)) + 
  geom_boxplot(aes(fill = tool_ref, pattern = texture),
               position = position_dodge2(preserve = "single", width = 0, padding = 0), 
               width = 0.5, color = "black", outliers = FALSE) +
  facet_grid(tools_db_comp ~ tools_labels_ref, scales = "free_y", space = "free") +
  scale_fill_manual(values = pal_10_q[match(c("DeepARG","fARGene", "RGI-DIAMOND"), tools_levels)]) +
  scale_y_discrete(drop = T) +
  xlab("%") +
  ylab("Pipeline covered") +
  theme_minimal() +
  theme5 +
  theme(panel.grid = element_blank(), 
        plot.margin = margin(0, 10, 0, 0, unit = "pt"),
        strip.text.x = element_text(size = general_size, vjust = 0, hjust = 0.5),
        strip.text.y = element_text(size = general_size, vjust = 0, hjust = 0))

cs2

p5 <- (cs2 + ggtitle("a")| cs11 + ggtitle("b")) + patchwork::plot_layout(widths = c(2, 1))


ggsave("code_R_analysis/output_plots/fig4.svg", p5, width = 180, height = 110, unit = "mm")



########################################################################################################################
########################################################################################################################


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
        #legend.title = element_text(size = general_size, angle = 90, vjust = 0.5),
        #legend.text = element_text(size = general_size, angle = 90, vjust = 0.5),
        legend.text = element_text(size = general_size))


rgi_query <- data.frame(gene = unique(c(lst$rgi.blast$query, lst$rgi.diamond$query, lst$rgi.diamond.prot$query)))
rgi_query <- rgi_query %>% mutate(dataset = ifelse(gene %in% intersect(intersect(lst$rgi.blast$query, lst$rgi.diamond$query), lst$rgi.diamond.prot$query), "all setups",
                                                   ifelse(gene %in% intersect(lst$rgi.blast$query, lst$rgi.diamond$query), "BLAST nt and DIAMOND (nt or aa)",
                                                          ifelse(gene %in% intersect(lst$rgi.diamond.prot$query, lst$rgi.diamond$query), "only\nDIAMOND\nnt and/or aa",
                                                                 ifelse(gene %in% intersect(lst$rgi.diamond.prot$query, lst$rgi.blast$query), "BLAST nt and DIAMOND (nt or aa)",
                                                                        ifelse(gene %in% lst$rgi.diamond.prot$query, "only\nDIAMOND\nnt and/or aa",
                                                                               ifelse(gene %in% lst$rgi.diamond$query, "only\nDIAMOND\nnt and/or aa",
                                                                                      ifelse(gene %in% lst$rgi.blast$query, "only BLAST\nnt", NA)))))))) %>%
  mutate(dataset = factor(dataset, levels = c("all setups", "only\nDIAMOND\nnt and/or aa", 
                                              "BLAST nt and DIAMOND (nt or aa)",  "only BLAST\nnt")))

rgi_query2 <- rgi_query %>% 
  group_by(dataset) %>% summarise(n = n()) %>% mutate(p = n/sum(n)) %>% 
  ungroup() %>%
  arrange(max(as.numeric(dataset)) - as.numeric(dataset)) %>%
  mutate(N = cumsum(n)) %>%
  mutate(diff_by1 = (N - lag(N, n = 1, default = 0))/2)

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
  theme(plot.margin = margin(0, -10, 0, 0, unit = "pt"))


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
  
  #print(db_query2)
  
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

deeparg_query <- data.frame(gene = unique(c(lst$deeparg.norm$query, lst$deeparg.norm.prot$query)))
deeparg_query <- deeparg_query %>% mutate(dataset = ifelse(gene %in% intersect(lst$deeparg.norm$query, lst$deeparg.norm.prot$query), "nt and aa",
                                                           ifelse(gene %in% lst$deeparg.norm$query, "nt", "aa"))) %>% 
  mutate(dataset = factor(dataset, levels = c("nt and aa", "nt", "aa")))


deeparg_over <- plot_overlaps(lst$deeparg.norm$query, lst$deeparg.norm.prot$query, "nt and aa", "only nt", "only aa", nam = "DeepARG")
fargene_over <- plot_overlaps(lst$fargene$query, lst$fargene.prot$query, "nt and aa", "only nt", "only aa", nam = "fARGene")
amrfinder_over <- plot_overlaps(lst$amrfinder.norm$query, lst$amrfinder.norm.prot$query, "nt and aa", "only nt", "only aa", nam = "AMRFinder-\nPlus")

plot_overlaps_same_tool <- 
  ((rgi_over + xlab("") + ggtitle("a")) / 
  (deeparg_over + theme(legend.position = "none")) /
  (fargene_over + theme(legend.position = "none")) /
  (amrfinder_over + xlab("ARGs"))) + 
  patchwork::plot_layout(heights = c(1, 1, 1, 1)) & 
    theme(plot.margin = margin(10, 10, 10, 10, unit = "pt"))


card_over <- plot_overlaps(lst$rgi.diamond.prot$query, lst$abricate.card.norm$query, "RGI and\nABRicate", "only RGI", "only ABRicate", nam = "CARD")
ncbi_over <- plot_overlaps(lst$amrfinder.norm.prot$query, lst$abricate.ncbi.norm$query, "AMRFinderPlus\nand ABRicate", "only AMRFinderPlus", "only ABRicate", nam = "NCBI")
resfinder_over <- plot_overlaps(lst$resfinder.norm$query, lst$abricate.resfinder.norm$query, "ResFinder\nand ABRicate", "only ResFinder", "only ABRicate", nam = "ResFinder")

plot_db <- 
  (card_over + ggtitle("b")) / 
  ncbi_over / 
  resfinder_over + xlab("ARGs") & 
  theme(plot.margin = margin(10, 10, 10, 10, unit = "pt"))


other_abricate <- unique(c(lst$abricate.argannot.norm$query, lst$abricate.card.norm$query, 
lst$abricate.ncbi.norm$query, lst$abricate.resfinder.norm$query))

abricate_query <- data.frame(rbind(lst$abricate.megares.norm[,c("tool", "query")], 
           lst$abricate.argannot.norm[,c("tool", "query")], 
           lst$abricate.card.norm[,c("tool", "query")], 
           lst$abricate.ncbi.norm[,c("tool", "query")], 
           lst$abricate.resfinder.norm[,c("tool", "query")]))

ntools_abricate0 <- abricate_query %>% 
  group_by(query) %>% 
  summarise(n_tools = n_distinct(tool))


ntools_abricate <- abricate_query %>% 
  group_by(query) %>% 
  summarise(n_tools = n_distinct(tool))
  
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
  theme(plot.margin = margin(10, 10, 10, 10, unit = "pt"))


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
  theme(plot.margin = margin(10, 10, 10, 10, unit = "pt"))

plot_megares <- (abricate_plot + ggtitle("c")) / megares_plot

#abricate_meg1 <- plot_overlaps(lst$abricate.resfinder.norm$query, lst$abricate.megares.norm$query, "Both", "Only\nResFinder", "Only\nMEGARes", nam = "ResFinder\nMEGARes")
#abricate_meg2 <- plot_overlaps(lst$abricate.argannot.norm$query, lst$abricate.megares.norm$query, "Both", "Only\nARGANNOT", "Only\nMEGARes", nam = "ARGANNOT\nMEGARes")
#abricate_meg3 <- plot_overlaps(lst$abricate.card.norm$query, lst$abricate.megares.norm$query, "Both", "Only\nCARD", "Only\nMEGARes", nam = "CARD\nMEGARes")
#abricate_meg4 <- plot_overlaps(lst$abricate.ncbi.norm$query, lst$abricate.megares.norm$query, "Both", "Only\nNCBI", "Only\nMEGARes", nam = "NCBI\nMEGARes")

#plot_megares <- 
#  abricate_meg2 /
#  abricate_meg3 /
#  abricate_meg4 /
#  abricate_meg1 + xlab("% of ARGs")

  
overlap_all_plots <- ((plot_overlaps_same_tool) | 
 (plot_db ) | 
(((abricate_plot + ggtitle("c")) / megares_plot) + patchwork::plot_layout(widths = c(1, 1)))) + 
  patchwork::plot_layout(widths = c(1, 1, 1))


ggsave("code_R_analysis/output_plots/overlaps_tool.svg", plot_overlaps_same_tool, width = 100, height = 90, unit = "mm")
ggsave("code_R_analysis/output_plots/overlaps_db.svg", plot_db, width = 100, height = 90, unit = "mm")
ggsave("code_R_analysis/output_plots/overlaps_megares.svg", plot_megares, width = 100, height = 90, unit = "mm")


ggsave("code_R_analysis/output_plots/overlaps_all_plots.svg", overlap_all_plots, width = 180, height = 150, unit = "mm")



# supplementary pan and core 

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
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.border = element_blank(),
        plot.margin = margin(5, 0, 0, 0, unit = "pt")) +
  scale_x_reverse(labels = label_comma()) 

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

pan_core_supp | pan_core_supp2

panel1_s <- (pan_core_supp + theme(legend.position = "none") | pan_core_supp2) + patchwork::plot_layout(widths = c(1,1))
p4_S <- (panel1 / (g_legend(p4a1) )) + patchwork::plot_layout(heights = c(12,1))

ggsave("code_R_analysis/output_plots/sup_pan_core.svg", p4_S, width = 120, height = 180, unit = "mm")


####



rich_plots <- list()
set.seed(2026)

for(j in 1:length(EN)){
  
  df_plot <- abundance_tool_sample %>%
    filter(tool %in% basic_tools) %>% 
    filter(habitat %in% EN[j]) %>% 
    group_by(habitat) %>% 
    mutate(N = n_distinct(sample),
           richness = richness)
  
  limits_plot <- df_plot %>% 
    ungroup() %>%
    group_by(habitat, tool) %>% 
    filter(richness < quantile(richness, 0.75) + 1.5*IQR(richness))
  
  
  
  rich_plots[[j]] <- df_plot %>%
    filter(tool %in% basic_tools) %>% 
    ggplot(aes(x = tools_labels, y = richness, fill = tools_db, pattern = texture)) + 
    geom_boxplot_pattern(position = position_dodge2(preserve = "single", width = 0.3, padding = 0), 
                         width = 1.3, pattern_color = "black", pattern_fill = "black", pattern_density = 0.000000001,
                         pattern_spacing = 0.2,
                         pattern_size =  0.3, color = "black", outliers = FALSE, outlier.shape = NA,
                         linewidth = 0.15) +
    scale_x_discrete(expand = expansion(add = 1)) +#, 
    geom_jitter(data = limits_plot %>%
                  slice_sample(n = min(100, df_plot$N[1])), width = 0.35, size = 0.4, alpha = 0.1) + 
    facet_grid(habitat ~ tools_db  , scales = "free_x", space = "free_x") +
    scale_y_continuous(labels = scales::comma, limits = c(0, max(limits_plot$richness))) + 
    scale_fill_manual(values = pal_7) +
    scale_pattern_manual(values = c('no' = 'none', 'yes' = 'stripe')) +
    theme1 + 
    theme(panel.border = element_blank(),
          strip.text.x = element_text(size = general_size, angle = 90, vjust = 0.5, hjust = 1),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          plot.margin = margin(5.5, 5.5, 5.5, 0, unit = "pt"))  
}


sup_rich_left <- 
  (rich_plots[[1]] + theme(axis.text.x = element_blank()) + xlab("") + ylab("")) / 
  (rich_plots[[2]] + theme(axis.text.x = element_blank(), strip.text.x = element_blank()) + xlab("") + ylab("")) / 
  (rich_plots[[3]] + theme(axis.text.x = element_blank(), strip.text.x = element_blank()) + xlab("") + ylab("Richness")) /
  (rich_plots[[4]] + theme(axis.text.x = element_blank(), strip.text.x = element_blank()) + xlab("") + ylab("")) /
  (rich_plots[[5]] + theme(strip.text.x = element_blank()) + xlab("") + ylab(""))

sup_rich_mid <- 
  patchwork::plot_spacer() /
  (rich_plots[[6]] + theme(axis.text.x = element_blank()) + xlab("") + ylab("")) / 
  (rich_plots[[7]] + theme(axis.text.x = element_blank(), strip.text.x = element_blank()) + xlab("") + ylab("")) /
  (rich_plots[[8]] + theme(axis.text.x = element_blank(), strip.text.x = element_blank()) + xlab("") + ylab("")) / 
  (rich_plots[[9]] + theme(strip.text.x = element_blank()) + xlab("") + ylab(""))
   
  

sup_rich_right <- 
  patchwork::plot_spacer() /
  (rich_plots[[10]] + theme(axis.text.x = element_blank()) + xlab("") + ylab("")) / 
  (rich_plots[[11]] + theme(axis.text.x = element_blank(), strip.text.x = element_blank()) + xlab("") + ylab("")) /
  (rich_plots[[12]] + theme(axis.text.x = element_blank(), strip.text.x = element_blank()) + xlab("") + ylab("")) / 
  (rich_plots[[13]] + theme(strip.text.x = element_blank()) + xlab("") + ylab("")) 


sup_rich <- ((sup_rich_left | sup_rich_mid | sup_rich_right) / 
  patchwork::wrap_elements(full = g_legend(a0 + 
                                             guides(fill = guide_legend(
                                               override.aes = list(
                                                 pattern = c(rep("none", 5),"stripe","none","stripe","none","stripe"),
                                                 fill  = pal_10_q[1:10]),
                                               nrow = 2), pattern = "none")))) +
  patchwork::plot_layout(heights = c(8, 1))

ggsave("code_R_analysis/output_plots/sup_richness.svg", sup_rich, width = 180, height = 180, unit = "mm")



###
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
    summarise(abundance = sum(abundance, richness = sum(richness)), texture = texture[1]) %>%
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
                  quantile(abundance, 0.75) + 1.5*IQR(abundance))) 
  
  abu_class_plots[[j]] <- df_class_summary_plots[[j]] %>% 
    filter(tool %in% basic_tools) %>% 
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
    facet_grid(gene ~ habitat, scales = "free", space = "free") +
    xlab("Relative abundance") +
    ylab("") + 
    labs(fill = "") + 
    scale_x_continuous(labels = scales::comma, expand = c(0,0)) +
    ggtitle("") +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
    theme1 + 
    theme(#strip.text.y = element_blank(),
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




sup_class_1 <- (((abu_class_plots[[1]] +  xlab("")) |
(abu_class_plots[[2]] +  xlab("") + ylab("")) |
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
ggsave("code_R_analysis/output_plots/sup_abundance_class_other_animal.svg", sup_class_2, width = 200, height = 180, unit = "mm")
ggsave("code_R_analysis/output_plots/sup_abundance_class_external_env.svg", sup_class_3, width = 200, height = 180, unit = "mm")




d2 <- d1 %>% mutate(texture = unigenes$texture[match(tool_ref, unigenes$tool)])

cs11_sup <- ggplot(d2 %>% 
                 filter(tool_ref %in% basic_tools, tool_comp %in% basic_tools) %>%
                 filter(!tool_ref %in% c("DeepARG","fARGene", "RGI-DIAMOND")) %>% 
                 filter(new_level %in% top_cso) %>%
                 mutate(new_level = gsub(" beta-lactamase","", new_level)) %>%
                 mutate(new_level = gsub("MFS efflux pump","MFS efflux", new_level)) %>%
                 mutate(new_level = gsub("efflux pump","efflux", new_level)), 
               aes(x = recall*100, y = new_level)) + 
  geom_boxplot_pattern(aes(fill = tool_ref, pattern = texture),
               position = position_dodge2(preserve = "single"),
               color = "black", outliers = FALSE) +
  scale_pattern_manual(values = c('no' = 'none', 'yes' = 'stripe')) +
  facet_grid(. ~ tool_ref, scales = "free_y") +
  scale_fill_manual(values = pal_10_q[match(tools_levels[!tools_levels %in% c("DeepARG","fARGene", "RGI-DIAMOND")], tools_levels)]) +
  scale_y_discrete(drop = FALSE) +
  xlab("%") +
  ylab("ARG class") + 
  theme_minimal() +
  theme5 +
  theme( panel.grid = element_blank(),
         strip.text.x = element_text(size = general_size, vjust = 0, hjust = 0.5),
         strip.text.y = element_text(size = general_size, vjust = 0, hjust = 0.5))

cs11_sup

cs2_sup <- ggplot(d2 %>% 
                filter(tool_ref %in% basic_tools, tool_comp %in% basic_tools) %>%
                filter(!tool_ref %in% c("RGI-DIAMOND", "fARGene","DeepARG")), #%>%
              aes(x = recall*100, y = tool_comp)) + 
  geom_boxplot_pattern(aes(fill = tool_ref, pattern = texture),
               position = position_dodge2(preserve = "single", width = 0, padding = 0), 
               width = 0.5, color = "black", outliers = FALSE) +
  facet_grid(tools_db_comp ~ tool_ref, scales = "free_y", space = "free") +
  scale_fill_manual(values = pal_10_q[match(tools_levels[!tools_levels %in% c("DeepARG","fARGene", "RGI-DIAMOND")], tools_levels)]) +
  scale_y_discrete(drop = T) +
  xlab("%") +
  ylab("Pipeline covered") +
  theme_minimal() +
  theme5 +
  theme(panel.grid = element_blank(), 
        plot.margin = margin(0, 10, 0, 0, unit = "pt"),
        strip.text.x = element_text(size = general_size, vjust = 0, hjust = 0.5),
        strip.text.y = element_text(size = general_size, vjust = 0, hjust = 0))

cs2_sup


ggsave("code_R_analysis/output_plots/fig_csc_sup_class.svg", cs11_sup, width = 180, height = 110, unit = "mm")
ggsave("code_R_analysis/output_plots/fig_csc_sup_tool.svg", cs2_sup, width = 180, height = 110, unit = "mm")


