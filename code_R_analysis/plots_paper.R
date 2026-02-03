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

setwd("~/Documents/GitHub/arg_compare/code_R_analysis") 
options(dplyr.summarise.inform = FALSE)

source("helper.R")

# Sourced gene classes 
# gene_classes_list is a vector for labeling the plots 
# gene_classes is the list with the names from the prepared datasets and used for shiny

general_size <- 10
lab_fn <- function(x) {
  x <- gsub("-", "-\n", x)
  x <- gsub(" ", "\n", x)
  x <- gsub("/", "/\n", x)
  x
}

# FORMAT PLOTS OLD COLORS
pal_7 <- brewer.pal(8, "BrBG")
pal_7 <- pal_7[c(4,5,3,6,2,7,1)]
pal_10_q <- pal_7[c(1,2,3,4,5,5,6,6,7,7)]

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

tools_levels <- c("DeepARG", "fARGene",
                    "ABRicate-ARGANNOT", "ABRicate-MEGARes",
                    "RGI-DIAMOND", "ABRicate-CARD",
                    "AMRFinderPlus", "ABRicate-NCBI",
                    "ResFinder", "ABRicate-ResFinder")

tools_labels <- c("DeepARG", "fARGene",
                    "ABRicate-ARGANNOT", "ABRicate-MEGARes",
                    "RGI-CARD", "ABRicate-CARD",
                    "AMRFinderPlus-NCBI", "ABRicate-NCBI",
                    "ResFinder", "ABRicate-ResFinder")

# we repeat the color but we add texture to the plots
tools_texture <- c("ABRicate-ARGANNOT", "ABRicate-MEGARes", "ABRicate-CARD", "ABRicate-NCBI", "ABRicate-ResFinder")

# environments that we are not interested in
not_env <- c("amplicon", "isolate", "built-environment" )

EN2  <- EN[!EN %in% not_env]
h2 <- c("humans","mammals","wastewater", "freshwater", "soil", "marine")


################################################################################

## DATASETS 

# df2 <- readRDS(file = "data/conversion_ARO_parent_new_level.rds")
# results per individual tool 
lst <- readRDS("output_abundance_diversity_resistome/results_tools.rds")

# metagenomes' metadata
metadata <- read.delim("../data/metadata_GMGC10.sample.meta.tsv")

# abundance aggreatated per ARO and per gene class considering all the detected genes per tool
# abundance <- readRDS("output_abundance_diversity_resistome/abundance_diversity.rds")
abundance <- bind_rows(readRDS("output_abundance_diversity_resistome/abundance_diversity_part1.rds"),
                       readRDS("output_abundance_diversity_resistome/abundance_diversity_part2.rds"),
                       readRDS("output_abundance_diversity_resistome/abundance_diversity_part3.rds"))

# needs to be loaded and cocatenated 3 times if the file is split
# abundance1 <- abundance[1:5327637,]
# abundance2 <- abundance[(5327637+1):10655274,]
# abundance3 <- abundance[(10655274+1):15982911,]

# saveRDS(abundance1, "output_abundance_diversity_resistome/abundance_diversity_part1.rds", compress = T)
# saveRDS(abundance2, "output_abundance_diversity_resistome/abundance_diversity_part2.rds", compress = T)
# saveRDS(abundance3, "output_abundance_diversity_resistome/abundance_diversity_part3.rds", compress = T)

# the column distinct_unigenes_rarefied <- alpha diversity (number of different genes after rarefaction) 
# no need to complete information here

abundance <- abundance %>% 
  mutate(unigenes = distinct_unigenes_rarefied) %>%  # alpha diversity (number of different genes after rarefaction) 
  mutate(habitat = factor(habitat, levels = EN), # convert to factors for ordering the plots
         habitat2 = factor(SO[habitat], levels = h2), # convert to factors for ordering the plots
         tool = factor(tool, levels = tools_levels)) %>%  # convert to factors for ordering the plots
  mutate(location = ifelse(habitat2 %in% c("humans","mammals","wastewater","built-environment"), "human-related","external")) %>% 
  mutate(location = factor(location, levels = c("human-related","external"))) %>% 
  filter(tool %in% tools_levels & !habitat %in% not_env) %>%
  filter(aggregation %in% "new_level") %>%  # take only gene class aggregation
  mutate(habitat2 = factor(as.character(habitat2), h2))

# normed10m is abundance, unigenes is the diversity (number of different genes )
# we need to complete information, we need to make sure all samples appear in all tools 

abundance_tool_sample <- abundance %>%
  group_by(tool, sample, habitat, habitat2) %>%  
  summarise(normed10m = sum(normed10m), unigenes = sum(unigenes)) %>% # sum the abundance and diversity
  ungroup() %>% 
  complete(sample, tool) %>% # complete with NAs
  left_join(abundance %>% select(sample, habitat, habitat2) %>% 
              distinct(), by = "sample") %>% # get habitat and habitat2 
  mutate(habitat  = coalesce(habitat.x, habitat.y), 
         habitat2 = coalesce(habitat2.x, habitat2.y)) %>%
  select(-habitat.x, -habitat.y, -habitat2.x, -habitat2.y) %>% 
  mutate(normed10m = replace_na(normed10m, 0)) %>%  # change NAs to 0
  mutate(unigenes = replace_na(unigenes, 0)) %>% # change NAs to 0
  arrange(tool, sample)
  
# normed10m is abundance, unigenes is the diversity (number of different genes )
# we need to complete information, we need to make sure all samples have all gene classes and appear in all tools 

abundance_class <- abundance %>% 
  ungroup() %>% 
  filter(tool %in% tools_levels, 
         !habitat %in% not_env, aggregation %in% "new_level") %>% 
  mutate(gene = factor(gene), sample = factor(sample), 
         tool = factor(tool, levels = tools_levels)) %>%
  select(!c(raw, raw_unique, scaled, distinct_unigenes_rarefied, 
            distinct_unigenes_raw, distinct_unigenes_raw_unique)) %>% 
  complete( sample, gene, tool, 
            fill = list(normed10m = 0,  unigenes = 0)) %>% ## complete for all samples, tools, and gene classes
  mutate(aggregation = "new_level", 
         gene = as.character(gene),
         habitat = metadata$habitat[match(sample, metadata$sample_id)],
         habitat2 = metadata$habitat2[match(sample, metadata$sample_id)],
         location = metadata$location[match(sample, metadata$sample_id)])  # not so dplyr way to fetch habitat and environment


## CORE RESISTOME 

core <- readRDS(file = "output_abundance_diversity_resistome/core_resistome.rds")
core <- core %>% 
  rename(new_level = new_level_centroid, 
         X = centroid) %>% 
  filter(tool %in% tools_levels, 
         !habitat %in% not_env) %>% 
  mutate(habitat = factor(habitat, levels = EN2), 
         tool = factor(tool, levels =  tools_levels)) %>% 
  mutate(tool = factor(tool, levels = tools_levels))

sumcore <- sum_core_adjust(core, 450, 0.5)

## PAN RESISTOME

pan <- readRDS(file = "output_abundance_diversity_resistome/pan_resistome.rds")
pan <- pan %>% 
  filter(tool %in% tools_levels, 
         !habitat %in% not_env, 
         aggregation %in% "new_level_centroid") %>% 
  mutate(habitat = factor(habitat, levels = EN2), 
         tool = factor(tool, levels =  tools_levels)) %>% 
  mutate(tool = factor(tool, levels = tools_levels))

# summarise pan-resistomes by tool, gene class, 
# take the median, mean, and standard deviation across subsamples

sumpan2 <- pan %>% ungroup() %>% 
  group_by(tool, habitat, aggregation, epoch) %>% 
  summarise(s = sum(unigenes)) %>%
  ungroup() %>% 
  group_by(tool, habitat, aggregation) %>% 
  summarise(md = median(s), mn = mean(s), sd = sd(s))

pan_core <- sumpan2 %>% 
  left_join((sumcore %>% 
               ungroup() %>% 
               group_by(tool, habitat) %>% 
               summarise(core = sum(unigenes))), by = c("tool", "habitat")) %>%
  mutate(core = ifelse(is.na(core), 0, core)) %>% 
  mutate(prop = core / md) %>%
  ungroup() %>% group_by(tool, habitat) %>% 
  mutate(texture = ifelse(tool %in% tools_texture, "yes", "no"))


## unigenes identified by tool

unigenes <- as_tibble(
  do.call(rbind, 
  lapply(lst, function(x) 
    x[,c("query","tool", "ARO", "parent", 
         "parent_description", "new_level", "id")]))) %>%
  filter(tool %in% tools_levels)  %>%
  mutate(tool = factor(tool, levels = tools_levels))


# Jaccard index, recall/class specific concordance, fnr/class specific non-overlap
# per class and tool
recall_fnr <- create_class_overlaps(unigenes)


# per tool
JI_all <- return_overlap_tools(unigenes)

# Jaccard index, recall/class specific concordance, fnr/class specific non-overlap
# JI_core_class <- create_class_overlaps_core(core, 0.5, 450, "human gut")

#JI_core_class <- return_overlap_tools(core %>% 
#  filter(cut %in% 0.5 & cnt > 450, habitat %in% "human gut") %>% 
#  rename(query = X))



pan_core_env <- c("human gut", "human skin", "pig gut", "dog gut", "wastewater",  "marine", "freshwater", "soil")


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


top20 <- c("van", "efflux pump", "cell wall charge", "rpoB", "tet RPG", "class A beta-lactamase", 
           "aph", "MFS efflux pump", "erm", "target-modifying enzyme", "tet enzyme")


df_abundance_class_human <- abundance_class %>% 
  filter(habitat %in% "human gut") %>% 
  mutate(total = normed10m) %>%
  mutate(gene = ifelse(gene %in% top20, gene, "Other")) %>% 
  mutate(texture = ifelse(tool %in% tools_texture, "yes", "no")) %>%
  mutate(gene = factor(gene, levels = c(top20, "Other"))) %>%
  ungroup() %>% 
  group_by(sample, tool, gene) %>%
  summarise(total = sum(total), texture = texture[1])

max_abun_class_human <- df_abundance_class_human %>% group_by(tool, gene) %>% summarise(q75 = quantile(total, .75)) %>% ungroup() %>% summarise(m = max(q75)) %>% pull()

################################################################################################
# PLOTS

# In case we need legends
square_legends <- g_legend(
    plot_count_genes_tool_legends(
      unigenes = unigenes , tools_for_figure = tools_levels, 
      general_size = general_size, # font size
      pal_10_q = pal_10_q, # pallet
      tool_label = tools_labels,  # labels for the tools
      tools_levels = tools_levels) +  # factor names of the tools
    ggtitle("A") + theme( title = element_text(size = general_size + 2, face = "bold")))


### 

p2a <- plot_count_genes_tool(
  unigenes = unigenes, 
  tools_for_figure = tools_levels, 
  general_size = general_size, 
  pal_10_q = pal_10_q, 
  tool_label = tools_labels , 
  tools_levels = tools_levels, 
  texture = tools_texture) + # which tools have stripes
  theme(panel.background = element_rect(colour = "black", fill = NA))

p2a <- p2a + 
  ggtitle("A") + theme(legend.position = "none", title = element_text(size = general_size + 2, face = "bold"))
  



p3a <- plot_total_abundance_diversity_new_version(
  dataset = abundance_tool_sample, # 
  tools_labels = tools_labels,  #
  tools_to_plot = tools_levels,  #
  environments_plot = EN[1:13], # habitats to plot (aggregated humans and mammals)
  general_size = general_size, # font size
  pal_10_q = pal_10_q , # pallet
  metric = "abundance", # metric (abundance or diversity)
  sd = 2025, # seed to plot random samples in the distribution 
  obs = 200,  # number of samples to plot as dots per environment
  texture = tools_texture,
  tools_levels = tools_levels) # texture for repeated color 

ps4 <- plot_total_abundance_diversity_new_version(
  dataset = abundance_tool_sample, # 
  tools_labels = tools_labels,  #
  tools_to_plot = tools_levels,  #
  environments_plot = EN[1:13], # habitats to plot (aggregated humans and mammals)
  general_size = general_size, # font size
  pal_10_q = pal_10_q , # pallet
  metric = "diversity", # metric (abundance or diversity)
  sd = 2025, # seed to plot random samples in the distribution 
  obs = 200,  # number of samples to plot as dots per environment
  texture = tools_texture,
  tools_levels = tools_levels) # texture for repeated color 



################################################################################################
# Fig 2
  
shape_tools <- rep(21, length(tools_labels))
shape_tools[tools_labels %in% tools_texture] <- 24

# legends circles and triangels 
circle_triangle_legends <- g_legend(pan_core %>% 
                             ggplot(aes(x = md, y = core +1, fill = tool,  shape = texture)) +
                             geom_point(size = 3, color = "black") +
                             scale_fill_manual(values = pal_10_q, labels = lab_fn(tools_labels), name = NULL) + 
                             scale_shape_manual(values = c(21, 24)) +
                             guides(
                                 fill = guide_legend(
                                   override.aes = list(
                                     shape = shape_tools,
                                     fill  = pal_10_q)), shape = "none") +
                             theme_minimal() +
                             theme(legend.position = "bottom",
                                   legend.text = element_text(size = general_size)))





p4a <- pan_core %>% select(!c(md, sd)) %>% pivot_longer(cols = c(mn, core), names_to = "metric", values_to = "value") %>% 
  mutate(metric = ifelse(metric %in% "mn", "Pan-resistome", metric)) %>%
  mutate(metric = ifelse(metric %in% "core", "Core-resistome", metric)) %>%
  mutate(metric = factor(metric, levels = c("Pan-resistome", "Core-resistome"))) %>%
  filter(habitat %in% c("human gut", "pig gut", "wastewater", "marine", "freshwater", "soil")) %>% 
  ggplot(aes(x = habitat, y =  value)) +
  geom_jitter(aes(fill = tool, shape = texture),  color = "black", stroke = 0.3, size = 2.5, width = .5) + 
  facet_grid(metric ~ habitat, scales = "free") +
  scale_fill_manual(values = pal_10_q) +
  scale_shape_manual(values = c("no" = 21, "yes" = 24)) +
  theme_minimal() +
  xlab("") +
  ylab("ARGs") + 
  theme(
    legend.position = "none",
    strip.text.x   = element_text(size = general_size),
    #strip.text.y   = element_blank(),
    legend.text = element_text(size = general_size),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    plot.margin = margin(0, 0, 0, 0, unit = "pt"),
    legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
    legend.margin = margin(0, 0, 0, 0, unit = "pt"),
    panel.spacing = unit(0, "pt"),
    title = element_text(size = general_size + 2, face = "bold"),
    axis.title = element_text(size = general_size + 1, face = "bold"),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = general_size),
    panel.border = element_blank(),   
    panel.background = element_rect(colour = "black", fill = NA))

p4a


ps18 <- pan_core %>% select(!c(md, sd)) %>% pivot_longer(cols = c(mn, core), names_to = "metric", values_to = "value") %>% 
  mutate(metric = ifelse(metric %in% "mn", "Pan-resistome", metric)) %>%
  mutate(metric = ifelse(metric %in% "core", "Core-resistome", metric)) %>%
  mutate(metric = factor(metric, levels = c("Pan-resistome", "Core-resistome"))) %>%
  filter(!habitat %in% c("human gut", "pig gut", "wastewater", "marine", "freshwater", "soil", "amplicon", "isolate", "built-environment")) %>% 
  ggplot(aes(x = habitat, y =  value)) +
  geom_jitter(aes(fill = tool, shape = texture),  color = "black", stroke = 0.3, size = 2.5, width = .5) + 
  facet_grid(metric ~ habitat, scales = "free") +
  scale_fill_manual(values = pal_10_q) +
  scale_shape_manual(values = c("no" = 21, "yes" = 24)) +
  theme_minimal() +
  xlab("") +
  ylab("ARGs") + 
  theme(
    legend.position = "none",
    strip.text.x   = element_text(size = general_size),
    #strip.text.y   = element_blank(),
    legend.text = element_text(size = general_size),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    plot.margin = margin(0, 0, 0, 0, unit = "pt"),
    legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
    legend.margin = margin(0, 0, 0, 0, unit = "pt"),
    panel.spacing = unit(0, "pt"),
    title = element_text(size = general_size + 2, face = "bold"),
    axis.title = element_text(size = general_size + 1, face = "bold"),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = general_size),
    panel.border = element_blank(),   
    panel.background = element_rect(colour = "black", fill = NA))

ps18


p_alluvial  <- plot_alluvial_classes(unigenes, levels_unigenes, 0.99, 0.005, tools_levels, tools_labels, tools_levels, pal_10_q, general_size, gene_classes_list) +
  theme(panel.background = element_rect(colour = "black", fill = NA))


##########
##########
##########


p4b <- alluvial_pan_core_env(sumcore, pan, "human gut", 
                             c("DeepARG", "fARGene", "ABRicate-MEGARes", 
                               "RGI-DIAMOND","AMRFinderPlus", "ResFinder"),
                             pal_10_complete, general_size) + 
  theme(title = element_blank())



#########
#########
#########



# p42 <- plot_recall_fnr_texture(recall_fnr, tools_levels, unique(unigenes$new_level), tools_labels, pal_10_q, general_size, tools_levels, tools_texture)
# 
# p5b1 <- one_tool_recall_fnr_id(unigenes, recall_fnr, "DeepARG", pal_10_q[3:5], general_size, tool_label, tools_levels_2) + 
#   ggtitle("DeepARG")
# 
# p5b10 <- one_tool_recall_fnr_id(unigenes, recall_fnr, "DeepARG", pal_10_q[3:5], general_size, tool_label, tools_levels_2) + 
#   ggtitle("DeepARG") + 
#   theme(legend.position = "none")
# 
# p5b2 <- one_tool_recall_fnr_id(unigenes, recall_fnr, "RGI-DIAMOND", pal_10_q[3:5], general_size, tool_label, tools_levels_2) +
#   ggtitle("RGI-CARD") + 
#   theme(legend.position = "none")
# 
# p5b3 <- one_tool_recall_fnr_id(unigenes, recall_fnr, "fARGene", pal_10_q[3:5], general_size, tool_label, tools_levels_2) +
#   ggtitle("fARGene") + 
#   theme(legend.position = "none")
# 
# p5b4 <- one_tool_recall_fnr_id(unigenes, recall_fnr, "AMRFinderPlus", pal_10_q[3:5], general_size, tool_label, tools_levels_2) +
#   ggtitle("AMRFinderPlus-NCBI") + 
#   theme(legend.position = "none")
# 
# p5 <- grid.arrange(p5b10, p5b2, g_legend(p5b1), 
#                    layout_matrix = rbind(c(1, 2), c(1, 2), c(1, 2), c(1, 2), 
#                                          c(1, 2), c(1, 2), c(1, 2), c(1, 2),
#                                          c(4,4)))



tool_second_order <- c("DeepARG", "fARGene", "ABRicate-MEGARes", "RGI-DIAMOND","AMRFinderPlus", "ResFinder", 
                       "ABRicate-CARD", "ABRicate-ResFinder", "ABRicate-ARGANNOT", "ABRicate-NCBI")
pal2 <- pal_10_q[match(tool_second_order, tools_levels)]



fig5a1 <- recall_fnr %>% 
  filter(new_level %in% top20) %>% 
  filter(tool_ref %in% c("DeepARG", "fARGene", "ABRicate-MEGARes", "RGI-DIAMOND")) %>%
  mutate(tool_ref = factor(as.character(tool_ref), levels = tool_second_order)) %>% 
  mutate(texture = ifelse(tool_comp %in% tools_texture, "yes", "no")) %>%
  mutate(new_level = factor(new_level, levels = top20)) %>%
  ggplot(aes(x = new_level, fill = tool_ref, y = recall)) +
    geom_violin() +
    #geom_boxplot(outlier.shape = NA, position = position_dodge2(preserve = "single")) + 
    geom_jitter(aes(fill = tool_comp, shape = texture), color = "black", 
                stroke = 0.3, size = 1.5, width = 0.4,  height = 0) + 
    scale_pattern_manual(values = c('no' = 'none', 'yes' = 'stripe')) +
    scale_fill_manual(values = pal2, 
                      labels = tools_levels) +
    scale_color_manual(values = pal2, 
                       labels = tools_levels) +
    scale_x_discrete(labels = function(x) {
      x <- gsub("-", "-\n", x)
      x <- gsub(" ", "\n", x)
      x}) + 
    scale_shape_manual(values = c("no" = 21, "yes" = 24)) +
    facet_grid(tool_ref ~ new_level, scales = "free_x") +
    scale_y_continuous(limits = c(-0.2,1.2), 
                       breaks = seq(0, 1, length.out = 3),
                       labels = scales::label_number()) + 
    ylab("CSTC") + 
    xlab("") + 
    theme(
      legend.position = "none",
      legend.text = element_text(size = general_size ),
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
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
      strip.text = element_blank(),
      strip.background = element_blank())

fig5a1

fig5a2 <- recall_fnr %>% 
  filter(!is.na(recall)) %>% 
  filter(tool_ref %in% c("DeepARG", "fARGene", "ABRicate-MEGARes", "RGI-DIAMOND")) %>%
  mutate(texture = ifelse(tool_ref %in% tools_texture, "yes", "no")) %>%
  ungroup() %>% 
  group_by(tool_ref, new_level, texture) %>% 
  summarise(recall = median(recall)) %>%
  ggplot(aes(x = "All classes medians", fill = tool_ref, y = recall)) +
  geom_violin() +
  geom_jitter(aes(shape = texture), size = 1.5, width = 0.4, height = 0) + 
  scale_pattern_manual(values = c('no' = 'none', 'yes' = 'stripe')) +
  scale_fill_manual(values = pal_10_q[tools_levels %in% c("DeepARG", "fARGene", "ABRicate-MEGARes", "RGI-DIAMOND","AMRFinderPlus", "ResFinder")], 
                    labels = tools_levels[tools_levels %in% c("DeepARG", "fARGene", "ABRicate-MEGARes", "RGI-DIAMOND","AMRFinderPlus", "ResFinder")]) +
  scale_color_manual(values = pal_10_q[tools_levels %in% c("DeepARG", "fARGene", "ABRicate-MEGARes", "RGI-DIAMOND","AMRFinderPlus", "ResFinder")], 
                     labels = tools_levels[tools_levels %in% c("DeepARG", "fARGene", "ABRicate-MEGARes", "RGI-DIAMOND","AMRFinderPlus", "ResFinder")]) +
  facet_grid(tool_ref ~ ., scales = "free_x") +
  scale_y_continuous(limits = c(-0.2,1.2), 
                     breaks = seq(0, 1, length.out = 3),
                     labels = scales::label_number()) +
  scale_x_discrete(labels = function(x) {
    x <- gsub("-", "-\n", x)
    x <- gsub(" ", "\n", x)
    x}) + 
  ylab("") + 
  xlab("") + 
  theme(
    legend.position = "none",
    legend.text = element_text(size = general_size ),
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    plot.margin = margin(0, 0, 0, 0, unit = "pt"),
    legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
    legend.margin = margin(0, 0, 0, 0, unit = "pt"),
    panel.spacing = unit(0, "pt"),
    title = element_text(size = general_size + 2, face = "bold"),
    axis.title = element_text(size = general_size + 1, face = "bold"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = general_size),
    axis.text.y = element_blank(),
    panel.background = element_rect(colour = "black", fill = NA),
    strip.text = element_text(size = general_size, face = "bold"),
    strip.background = element_blank())


fig5a2


fig5b1 <- recall_fnr %>% 
  filter(new_level %in% top20) %>% 
  filter(tool_ref %in% c("DeepARG", "fARGene", "ABRicate-MEGARes", "RGI-DIAMOND")) %>%
  mutate(tool_ref = factor(as.character(tool_ref), levels = tool_second_order)) %>% 
  mutate(texture = ifelse(tool_comp %in% tools_texture, "yes", "no")) %>%
  mutate(new_level = factor(new_level, levels = top20)) %>%
  ggplot(aes(x = new_level, fill = tool_ref, y = fnr)) +
  geom_violin() +
  geom_jitter(aes(fill = tool_comp, shape = texture), color = "black", stroke = 0.3, size = 1.5, width = 0.4, height = 0) + 
  scale_pattern_manual(values = c('no' = 'none', 'yes' = 'stripe')) +
  scale_fill_manual(values = pal2, 
                    labels = tools_levels) +
  scale_color_manual(values = pal2, 
                     labels = tools_levels) +
  scale_x_discrete(labels = function(x) {
    x <- gsub("-", "-\n", x)
    x <- gsub(" ", "\n", x)
    x}) + 
  scale_shape_manual(values = c("no" = 21, "yes" = 24)) +
  facet_grid(tool_ref ~ new_level, scales = "free_x") +
  scale_y_continuous(limits = c(-0.2,1.2), 
                     breaks = seq(0, 1, length.out = 3),
                     labels = scales::label_number()) + 
  ylab("CSNO") + 
  xlab("") + 
  theme(
    legend.position = "none",
    legend.text = element_text(size = general_size ),
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
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
    strip.text = element_blank(),
    strip.background = element_blank())

fig5b1

fig5b2 <- recall_fnr %>% 
  filter(!is.na(recall)) %>% 
  filter(tool_ref %in% c("DeepARG", "fARGene", "ABRicate-MEGARes", "RGI-DIAMOND")) %>%
  mutate(texture = ifelse(tool_ref %in% tools_texture, "yes", "no")) %>%
  ungroup() %>% 
  group_by(tool_ref, new_level, texture) %>% 
  summarise(fnr = median(fnr)) %>%
  ggplot(aes(x = "Class medians", fill = tool_ref, y = fnr)) +
  geom_violin() +
  #geom_boxplot(outlier.shape = NA, position = position_dodge2(preserve = "single")) + 
  geom_jitter(aes(shape = texture), size = 1.5, width = 0.4, height = 0) + 
  scale_pattern_manual(values = c('no' = 'none', 'yes' = 'stripe')) +
  scale_fill_manual(values = pal_10_q[tools_levels %in% c("DeepARG", "fARGene", "ABRicate-MEGARes", "RGI-DIAMOND","AMRFinderPlus", "ResFinder")], 
                    labels = tools_levels[tools_levels %in% c("DeepARG", "fARGene", "ABRicate-MEGARes", "RGI-DIAMOND","AMRFinderPlus", "ResFinder")]) +
  scale_color_manual(values = pal_10_q[tools_levels %in% c("DeepARG", "fARGene", "ABRicate-MEGARes", "RGI-DIAMOND","AMRFinderPlus", "ResFinder")], 
                     labels = tools_levels[tools_levels %in% c("DeepARG", "fARGene", "ABRicate-MEGARes", "RGI-DIAMOND","AMRFinderPlus", "ResFinder")]) +
  facet_grid(tool_ref ~ ., scales = "free_x") +
  scale_y_continuous(limits = c(-0.2,1.2), 
                     breaks = seq(0, 1, length.out = 3),
                     labels = scales::label_number()) +
  scale_x_discrete(labels = function(x) {
    x <- gsub("-", "-\n", x)
    x <- gsub(" ", "\n", x)
    x}) + 
  ylab("") + 
  xlab("") + 
  theme(
    legend.position = "none",
    legend.text = element_text(size = general_size ),
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    plot.margin = margin(0, 0, 0, 0, unit = "pt"),
    legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
    legend.margin = margin(0, 0, 0, 0, unit = "pt"),
    panel.spacing = unit(0, "pt"),
    title = element_text(size = general_size + 2, face = "bold"),
    axis.title = element_text(size = general_size + 1, face = "bold"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = general_size),
    axis.text.y = element_blank(),
    panel.background = element_rect(colour = "black", fill = NA),
    strip.text = element_text(size = general_size, face = "bold"),
    strip.background = element_blank())


fig5b2



fig5c1 <- unigenes %>% 
  filter(tool %in% c("DeepARG",  "RGI-DIAMOND"), new_level %in% top20) %>%
  mutate(new_level = factor(new_level, levels = top20)) %>% 
  ggplot(aes(x = new_level, fill = tool, y = id )) +
  geom_violin() +
  #geom_jitter(size = 0.01, width = 0.3, alpha = 0.02, color = "grey") + 
  scale_pattern_manual(values = c('no' = 'none', 'yes' = 'stripe')) +
  scale_fill_manual(values = pal_10_q[tools_levels %in% c("DeepARG",  "RGI-DIAMOND")], 
                    labels = tools_levels[tools_levels %in% c("DeepARG",  "RGI-DIAMOND")]) +
  scale_color_manual(values = pal_10_q[tools_levels %in% c("DeepARG",  "RGI-DIAMOND")], 
                     labels = tools_levels[tools_levels %in% c("DeepARG",  "RGI-DIAMOND")]) +
  facet_grid(tool ~ new_level, scales = "free_x") +
  scale_y_continuous(limits = c(20, 101), 
                     breaks = seq(20, 100, length.out = 3),
                     labels = scales::label_number()) +
  scale_x_discrete(labels = function(x) {
    x <- gsub("-", "-\n", x)
    x <- gsub(" ", "\n", x)
    x}) + 
  ylab("% Identity") + 
  xlab("") + 
  theme(
    legend.position = "none",
    legend.text = element_text(size = general_size ),
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    plot.margin = margin(0, 0, 0, 0, unit = "pt"),
    legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
    legend.margin = margin(0, 0, 0, 0, unit = "pt"),
    panel.spacing = unit(0, "pt"),
    title = element_text(size = general_size + 2, face = "bold"),
    axis.title = element_text(size = general_size + 1, face = "bold"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = general_size),
    axis.text.y = element_text(size = general_size ),
    panel.background = element_rect(colour = "black", fill = NA),
    strip.text = element_blank(),
    strip.background = element_blank())


fig5c1


fig5c2 <- unigenes %>% 
  filter(tool %in% c("DeepARG",  "RGI-DIAMOND")) %>%
  ungroup() %>% 
  group_by(tool, new_level) %>% 
  summarise(id = median(id)) %>% 
  ggplot(aes(x = "All class medians", fill = tool, y = id )) +
  geom_violin() +
  geom_jitter(size = 1.5, width = 0.4, height = 0) +
  scale_pattern_manual(values = c('no' = 'none', 'yes' = 'stripe')) +
  scale_fill_manual(values = pal_10_q[tools_levels %in% c("DeepARG",  "RGI-DIAMOND")], 
                    labels = tools_levels[tools_levels %in% c("DeepARG",  "RGI-DIAMOND")]) +
  scale_color_manual(values = pal_10_q[tools_levels %in% c("DeepARG",  "RGI-DIAMOND")], 
                     labels = tools_levels[tools_levels %in% c("DeepARG",  "RGI-DIAMOND")]) +
  facet_grid(tool ~ ., scales = "free_x") +
  scale_y_continuous(limits = c(20,101), 
                     breaks = seq(20, 100, length.out = 3),
                     labels = scales::label_number()) +
  scale_x_discrete(labels = function(x) {
    x <- gsub("-", "-\n", x)
    x <- gsub(" ", "\n", x)
    x}) + 
  ylab("") + 
  xlab("") + 
  theme(
    legend.position = "none",
    legend.text = element_text(size = general_size ),
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    plot.margin = margin(0, 0, 0, 0, unit = "pt"),
    legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
    legend.margin = margin(0, 0, 0, 0, unit = "pt"),
    panel.spacing = unit(0, "pt"),
    title = element_text(size = general_size + 2, face = "bold"),
    axis.title = element_text(size = general_size + 1, face = "bold"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = general_size),
    axis.text.y = element_blank(),
    panel.background = element_rect(colour = "black", fill = NA),
    strip.text = element_text(size = general_size, face = "bold"),
    strip.background = element_blank())

fig5c2

plot5a <- plot_grid( fig5a1 + theme(axis.text.x = element_blank()), 
                     fig5a2 + theme(axis.text.x = element_blank()), 
                     ncol = 2,
                     align = "h",
                     axis = "tb",
                     rel_widths = c(3, 1))

plot5b <- plot_grid( fig5b1 + theme(axis.text.x = element_blank()), 
                     fig5b2 + theme(axis.text.x = element_blank()), 
                     ncol = 2,
                     align = "h",
                     axis = "tb",
                     rel_widths = c(3, 1))



plot5c <- plot_grid( fig5c1, 
                     fig5c2, 
                     ncol = 2,
                     align = "h",
                     axis = "tb",
                     rel_widths = c(3, 1))

p5_legend <- plot_grid( circle_triangle_legends, NULL,
                     ncol = 2,
                     align = "h",
                     axis = "tb",
                     rel_widths = c(6, 1))

p5 <- plot_grid( plot5a, plot5b, plot5c, p5_legend, ncol = 1,
                     align = "v",
                     axis = "lr",
                     rel_heights = c(14, 14, 9, 3))

ggsave("output_plots/fig5.svg", p5, width = 180, height = 210, unit = "mm")



abu_class_legend0 <- df_abundance_class_human %>%
  ggplot(aes(x = gene, y = total, fill = tool, pattern = texture)) + 
  stat_summary(fun.data = calc_boxplot_stat, geom="boxplot", 
               position = position_dodge2(preserve = "single"), 
               color = "black", linewidth = 0.2, 
               outlier.shape = NA, outlier.size = 0) +
  scale_fill_manual(values = pal_10_q, 
                    labels = lab_fn(tools_levels)) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, floor(max_abun_class_human/1000) * 1000, length.out = 4),
                     labels = scales::label_number()) + 
  theme_minimal() +
  facet_grid(. ~ gene, scales = "free_x", space = "free_x") + 
  scale_x_discrete(labels = function(x) {
    x <- gsub("-", "-\n", x)
    x <- gsub(" ", "\n", x)
    x})  + 
  ylab("Relative abundance") +
  xlab("") +
  labs(fill = "") +
  ggtitle("B") +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = general_size ),
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
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
    strip.text = element_blank())


abu_class_human <- df_abundance_class_human %>%
  ggplot(aes(x = gene, y = total, fill = tool, pattern = texture)) + 
  stat_summary(fun.data = calc_boxplot_stat, geom="boxplot_pattern", 
               position = position_dodge2(preserve = "single"), 
               color = "black", linewidth = 0.2, 
               pattern_color = "white",
               pattern_density = pattern_density, 
               pattern_spacing = pattern_spacing, 
               pattern_fill = pattern_fill,
               pattern_size = pattern_size,
               pattern_key_scale_factor = 1.2, outlier.shape = NA, outlier.size = 0) +
  scale_pattern_manual(values = c('no' = 'none', 'yes' = 'stripe')) +
  scale_fill_manual(values = pal_10_q, 
                    labels = lab_fn(tools_levels)) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, floor(max_abun_class_human/1000) * 1000, length.out = 4),
                     labels = scales::label_number()) + 
  theme_minimal() +
  facet_grid(. ~ gene, scales = "free_x", space = "free_x") + 
  scale_x_discrete(labels = function(x) {
    x <- gsub("-", "-\n", x)
    x <- gsub(" ", "\n", x)
    x})  + 
  ylab("Relative abundance") +
  xlab("") +
  labs(fill = "") +
  ggtitle("B") +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = general_size ),
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
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
    strip.text = element_blank())




########################################################################################
type_tools <- c("solid", "solid", "solid", "dotted")



id_plot <- unigenes %>% 
  ungroup() %>%
  mutate(tool_2 = ifelse(tool %in% c("ABRicate-ARGANNOT", "ABRicate-MEGARes", "ABRicate-CARD",
                                     "ABRicate-ResFinder", "ABRicate-NCBI", "ResFinder"), "ABRicate/ResFinder", as.character(tool))) %>% 
  filter(!((tool_2 %in% "AMRFinderPlus") & (query %in% lst$amrfinder.norm.prot$query[lst$amrfinder.norm.prot$Method %in% "HMM"]))) %>% 
  filter(!tool_2 %in% "fARGene") %>% 
  mutate(tool_2 = factor(tool_2, levels = c("DeepARG", "RGI-DIAMOND", "AMRFinderPlus", "ABRicate/ResFinder"))) %>% 
  mutate(data_type = ifelse(tool_2 %in% c("DeepARG", "RGI-DIAMOND", "AMRFinderPlus"), "amino acid", "nucleotide")) %>% 
  ggplot(aes(x = id , color = tool_2, fill = tool_2, linetype = data_type)) +
  stat_ecdf(geom = "step", linewidth = 1) + 
  scale_linetype_manual(values = c("solid","dotted")) +
  theme_minimal() + 
  scale_color_manual(values = pal_10_q[c(1,5,7,9)], 
                     labels = lab_fn(c(tools_labels[c(1,5,7)],"ABRicate/ResFinder")), 
                     guide = guide_legend(nrow = 2)) +
  ggtitle("B") + 
  guides(
    color = guide_legend(
      override.aes = list(
        linetype = type_tools,
        color  = pal_10_q[c(1,5,7,9)]),
        nrow = 2), linetype = "none") +
  ylab("ECDF") + 
  xlab("% Identity") +
  labs(color = "", fill = "", linetype = "") +
  scale_y_continuous(limits = c(0,1), expand = c(0.01, 0.01)) + 
  theme(legend.position = "bottom",
        panel.border = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = general_size),
        axis.text.y = element_text(size = general_size),
        strip.text = element_text(size = general_size, face = "bold"),
        plot.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.margin = margin(0, 0, 0, 0, unit = "pt"),
        panel.spacing = unit(0, "pt"),
        legend.text = element_text(size = general_size),
        title = element_text(size = general_size + 2, face = "bold"),
        axis.title = element_text(size = general_size + 1, face = "bold"),
        panel.background = element_rect(colour = "black", fill = NA))

id_plot


fig2ab <- plot_grid(
  p2a, id_plot,
  ncol = 2, align = "hv", axis = "tblr",
  rel_widths = c(7,4)
)


p2_1 <- plot_grid(fig2ab,
                  p_alluvial + ggtitle("C"),
                  ncol = 1,
                  rel_heights = c(1,1.4))

ggsave("output_plots/fig2.svg", p2_1, width = 180, height = 210, unit = "mm")
ggsave("output_plots/fig2ab.svg", fig2ab, width = 180, height = 70, unit = "mm")
ggsave("output_plots/fig2c.svg", p_alluvial + ggtitle("C"), width = 180, height = 140, unit = "mm")


fig_tot_abund_human_abund <- plot_grid(
  p3a + theme(legend.position = " none") + ggtitle("A"), 
  abu_class_human + theme(legend.position = " none") + ggtitle("B"),
  g_legend(abu_class_legend0),
  nrow = 3,
  align = "v",
  axis = "lr",
  rel_heights = c(14, 8, 2))


ggsave("output_plots/fig3.svg", fig_tot_abund_human_abund, width = 180, height = 210, unit = "mm")
ggsave("output_plots/fig3_alternative_legend.svg", g_legend(abu_class_human), width = 180, height = 30, unit = "mm")

ggsave("output_plots/fig_diversity.svg", p2diversity, width = 180, height = 140, unit = "mm")

p4 <- plot_grid(p4a + ggtitle("A"),
                  circle_triangle_legends,
                  p4b + ggtitle("B") + xlab("Pan- and core-resistomes"),
                  ncol = 1,
                  align = "hv",
                  axis = "tblr",
                  rel_heights = c(7,2,15))

ggsave("output_plots/fig4.svg", p4, width = 180, height = 210, unit = "mm")

# OVERLAPS 

theme_overlap <- theme_minimal() + 
  theme(legend.position = "none",
        panel.border = element_blank(),
        axis.text.x = element_text(size = general_size),
        axis.text.y = element_text( size = general_size),
        strip.text = element_text(size = general_size, face = "bold"),
        plot.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.margin = margin(0, 0, 0, 0, unit = "pt"),
        panel.spacing = unit(0, "pt"),
        legend.text = element_text(size = general_size),
        title = element_text(size = general_size + 2, face = "bold"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title = element_text(size = general_size + 1, face = "bold"),
        panel.background = element_rect(colour = "black", fill = NA))


rgi_query <- data.frame(gene = unique(c(lst$rgi.blast$query, lst$rgi.diamond$query, lst$rgi.diamond.prot$query)))
rgi_query <- rgi_query %>% mutate(dataset = ifelse(gene %in% intersect(intersect(lst$rgi.blast$query, lst$rgi.diamond$query), lst$rgi.diamond.prot$query), "All",
                                            ifelse(gene %in% intersect(lst$rgi.blast$query, lst$rgi.diamond$query), "B-nt and D-nt",
                                            ifelse(gene %in% intersect(lst$rgi.diamond.prot$query, lst$rgi.diamond$query), "D-nt and D-aa",
                                            ifelse(gene %in% intersect(lst$rgi.diamond.prot$query, lst$rgi.blast$query), "B-nt and D-aa",
                                            ifelse(gene %in% lst$rgi.diamond.prot$query, "D-aa",
                                            ifelse(gene %in% lst$rgi.diamond$query, "D-nt",
                                            ifelse(gene %in% lst$rgi.blast$query, "B-nt", NA)))))))) %>%
  mutate(dataset = factor(dataset, levels = c("All", "B-nt and D-nt", "D-nt and D-aa", "B-nt and D-aa", "D-aa", "D-nt", "B-nt")))
                                            
rgi_over <- rgi_query %>% 
  group_by(dataset) %>% summarise(n = n()) %>% mutate(p = n/sum(n)) %>%
  ggplot(aes(x = dataset, y = n, fill = dataset)) + 
  geom_col() + 
  scale_fill_manual(values = pal_10_q[3:10])  +
  geom_text(aes(label = scales::percent(p, accuracy = 0.1)), vjust = -0.3, size = general_size / .pt) +
  ylab("ARGs") + 
  xlab("") +
  scale_x_discrete(labels = function(x) {
    x <- gsub(" ", "\n", x)
    x}) +
  theme_overlap


deeparg_query <- data.frame(gene = unique(c(lst$deeparg.norm$query, lst$deeparg.norm.prot$query)))
deeparg_query <- deeparg_query %>% mutate(dataset = ifelse(gene %in% intersect(lst$deeparg.norm$query, lst$deeparg.norm.prot$query), "Both DeepARG",
                                                   ifelse(gene %in% lst$deeparg.norm$query, "nt", "aa"))) %>% 
  mutate(dataset = factor(dataset, levels = c("Both DeepARG", "nt", "aa")))


deeparg_over <- plot_overlaps(lst$deeparg.norm$query, lst$deeparg.norm.prot$query, "Both", "nt", "aa")
fargene_over <- plot_overlaps(lst$fargene$query, lst$fargene.prot$query, "Both", "nt", "aa")
amrfinder_over <- plot_overlaps(lst$amrfinder.norm$query, lst$amrfinder.norm.prot$query, "Both", "nt", "aa")

plot_overlaps_same_tool <- plot_grid(
  rgi_over + ggtitle("A"), 
  deeparg_over + ggtitle("B"), 
  fargene_over + ggtitle("C"), 
  amrfinder_over + ggtitle("D"),
  ncol = 2,
  align = "hv"
)


card_over <- plot_overlaps(lst$rgi.diamond.prot$query, lst$abricate.card.norm$query, "Both", "RGI", "ABRicate")
ncbi_over <- plot_overlaps(lst$amrfinder.norm.prot$query, lst$abricate.ncbi.norm$query, "Both", "AMRFinderPlus", "ABRicate")
resfinder_over <- plot_overlaps(lst$resfinder.norm$query, lst$abricate.resfinder.norm$query, "Both", "ResFinder", "ABRicate")

plot_db <- plot_grid(
  card_over  + ggtitle("A"), 
  ncbi_over + ggtitle("B"), 
  resfinder_over + ggtitle("C"),
  ncol = 2,
  align = "hv"
)


abricate_meg1 <- plot_overlaps(lst$abricate.resfinder.norm$query, lst$abricate.megares.norm$query, "Both", "ResFinder", "MEGARes")
abricate_meg2 <- plot_overlaps(lst$abricate.argannot.norm$query, lst$abricate.megares.norm$query, "Both", "ARGANNOT", "MEGARes")
abricate_meg3 <- plot_overlaps(lst$abricate.card.norm$query, lst$abricate.megares.norm$query, "Both", "CARD", "MEGARes")
abricate_meg4 <- plot_overlaps(lst$abricate.ncbi.norm$query, lst$abricate.megares.norm$query, "Both", "NCBI", "MEGARes")

plot_megares <- plot_grid(
  abricate_meg2 + ggtitle("A"), 
  abricate_meg3 + ggtitle("B"), 
  abricate_meg4 + ggtitle("C"),
  abricate_meg1  + ggtitle("D"), 
  ncol = 2,
  align = "hv"
)


ggsave("output_plots/overlaps_same_tool.svg", plot_overlaps_same_tool, width = 180, height = 160, unit = "mm")
ggsave("output_plots/overlaps_db.svg", plot_db, width = 180, height = 160, unit = "mm")
ggsave("output_plots/overlaps_megares.svg", plot_megares, width = 160, height = 210, unit = "mm")



alluvial_all <- list()
abundance_class_all <- list()

for(e in EN[-c(14:16)]){
  alluvial_all[[e]] <- alluvial_pan_core_env(sumcore, pan, e, tools_levels, pal_10_complete, general_size)
  abundance_class_all[[e]] <- plot_abundance_class_environment(abundance_class, e, general_size , pal_10_q)  
}

for(e in EN[-c(14:16)]){
  ggsave(gsub(" ", "", paste0("output_plots/s_", e, ".svg")), alluvial_all[[e]], width = 180, height = 180, unit = "mm")
  ggsave(gsub(" ", "", paste0("output_plots/s2_", e, ".svg")), abundance_class_all[[e]], width = 180, height = 100, unit = "mm")
}


###






habitat %in% c("wastewater", "marine", "freshwater", "soil")

sumpan2 %>% ungroup() %>% group_by(habitat) %>% mutate(deep = mn[tool == "AMRFinderPlus"][1]) %>% 
  ungroup() %>%
  filter(tool %in% c("ABRicate-CARD")) %>%
  mutate(over_deepARG = deep / mn) %>%
  summarise(mn = mean(over_deepARG))

sumpan2 %>% ungroup() %>% group_by(habitat) %>% mutate(deep = mn[tool == "RGI-DIAMOND"][1]) %>% 
  ungroup() %>%
  filter(tool %in% c("ABRicate-CARD")) %>%
  mutate(over_deepARG = deep / mn) %>%
  summarise(mn = mean(over_deepARG))


sumpan2 %>% ungroup() %>% group_by(habitat) %>% mutate(deep = mn[tool == "ResFinder"][1]) %>% 
  ungroup() %>%
  filter(tool %in% c("ABRicate-ResFinder")) %>%
  mutate(over_deepARG = deep / mn) %>%
  summarise(mn = mean(over_deepARG))



sumcore <- sum_core_adjust(core, 450, 0.5)
sumcore %>% 
  ungroup() %>%
  group_by(tool, habitat) %>% 
  summarise(unigenes = sum(unigenes)) %>%
  group_by(habitat) %>% mutate(deep = unigenes[tool == "DeepARG"][1], rgi = unigenes[tool == "RGI-DIAMOND"][1]) %>% 
  ungroup() %>%
  filter(!tool %in% c("DeepARG", "RGI-DIAMOND")) %>%
  filter(!habitat %in% c("wastewater", "marine", "freshwater", "soil", "human oral", "human skin", "human vagina", "human nose", "human oral")) %>% 
  mutate(over_deepARG = deep / unigenes, over_rgi = rgi / unigenes) %>%
  group_by(habitat) %>%
  summarise(min_deep = min(over_deepARG, na.rm = T), max_deep = max(over_deepARG, na.rm = T), 
            min_rgi = min(over_rgi, na.rm = T), max_rgi = max(over_rgi, na.rm = T)) %>% 
  ungroup() %>%
  summarise(min_deep = min(min_deep, na.rm = T), max_deep = max(max_deep, na.rm = T), 
            min_rgi = min(min_rgi, na.rm = T), max_rgi = max(max_rgi, na.rm = T))


sumcore %>% 
  ungroup() %>%
  group_by(tool, habitat) %>% 
  summarise(unigenes = sum(unigenes)) %>%
  group_by(habitat) %>% mutate(deep = unigenes[tool == "DeepARG"][1], rgi = unigenes[tool == "RGI-DIAMOND"][1]) %>% 
  ungroup() %>%
  filter(!tool %in% c("DeepARG", "RGI-DIAMOND")) %>%
  filter(habitat %in% c("human gut")) %>% 
  mutate(over_deepARG = deep / unigenes, over_rgi = rgi / unigenes) 


sumcore %>% 
  ungroup() %>%
  group_by(habitat) %>% mutate(deep = unigenes[tool == "DeepARG"][1], rgi = unigenes[tool == "RGI-DIAMOND"][1]) %>% 
  ungroup() %>%
  filter(habitat %in% c("pig gut"), new_level %in% "van")

sumcore %>% 
  ungroup() %>%
  group_by(tool, habitat) %>% 
  summarise(unigenes = sum(unigenes)) %>%
  group_by(habitat) %>% mutate(deep = unigenes[tool == "DeepARG"][1], rgi = unigenes[tool == "RGI-DIAMOND"][1]) %>% 
  ungroup() %>%
  filter(habitat %in% c("pig gut")) %>% 
  filter(!tool %in% c("DeepARG", "RGI-DIAMOND")) %>%
  mutate(over_deepARG = deep / unigenes, over_rgi = rgi / unigenes) 

sumcore %>% 
  ungroup() %>%
  group_by(tool, habitat) %>% 
  summarise(unigenes = sum(unigenes)) %>% 
  filter(habitat %in% "marine")


sumpan2 %>% ungroup() %>% group_by(habitat) %>% mutate(deep = mn[tool == "DeepARG"][1], rgi = mn[tool == "RGI-DIAMOND"][1]) %>% 
  ungroup() %>%
  filter(!tool %in% c("DeepARG", "RGI-DIAMOND")) %>%
  mutate(over_deepARG = deep / mn, over_rgi = rgi / mn) %>%
  group_by(habitat) %>%
  summarise(min_deep = min(over_deepARG), max_deep = max(over_deepARG), 
            min_rgi = min(over_rgi), max_rgi = max(over_rgi)) %>% 
  ungroup() %>%
  summarise(min_deep = min(min_deep), max_deep = max(max_deep), 
            min_rgi = min(min_rgi), max_rgi = max(max_rgi))



sumcore %>% 
  ungroup() %>%
  group_by(tool, habitat) %>% 
  summarise(unigenes = sum(unigenes)) %>%
  group_by(habitat) %>% mutate(deep = unigenes[tool == "DeepARG"][1], rgi = unigenes[tool == "RGI-DIAMOND"][1]) %>% 
  ungroup() %>%
  filter(!tool %in% c("DeepARG", "RGI-DIAMOND")) %>%
  filter(!habitat %in% c("wastewater", "marine", "freshwater", "soil", "human oral", "human skin", "human vagina", "human nose", "human oral")) %>% 
  mutate(over_deepARG = deep / unigenes, over_rgi = rgi / unigenes) %>%
  group_by(habitat) %>%
  summarise(min_deep = min(over_deepARG, na.rm = T), max_deep = max(over_deepARG, na.rm = T), 
            min_rgi = min(over_rgi, na.rm = T), max_rgi = max(over_rgi, na.rm = T)) %>% 
  ungroup() %>%
  summarise(min_deep = min(min_deep, na.rm = T), max_deep = max(max_deep, na.rm = T), 
            min_rgi = min(min_rgi, na.rm = T), max_rgi = max(max_rgi, na.rm = T))








core %>% filter(habitat %in% "human gut", cut %in% 0.5, cnt > 450) %>% 
  ungroup() %>% 
  group_by(tool, new_level) %>%
  summarise(n_tool = n_distinct(X)) %>%
  mutate(p = n_tool/sum(n_tool)) %>%
  arrange(desc(n_tool))%>%
  filter(new_level %in% "tet RPG")

28 / 12
core %>% filter(habitat %in% "pig gut", cut %in% 0.5, cnt > 450) %>% 
  ungroup() %>% 
  group_by(tool, new_level) %>%
  summarise(n_tool = n_distinct(X)) %>%
  mutate(p = n_tool/sum(n_tool)) %>%
  arrange(desc(n_tool))

core %>% filter(habitat %in% "pig gut", cut %in% 0.5, cnt > 450) %>% 
  ungroup() %>% 
  group_by(tool) %>%
  summarise(n_tool = n_distinct(X)) 
  arrange(desc(n_tool))

core %>% filter(habitat %in% "human gut", cut %in% 0.5, cnt > 450) %>% 
  ungroup() %>% 
  group_by(X, new_level) %>%
  summarise(n_tool = n_distinct(tool)) %>% 
  mutate(n_tool50 = ifelse(n_tool>=5,"yes","no")) %>% 
  ungroup() %>% 
  group_by( new_level) %>% 
  mutate(N = n_distinct(X)) %>%
  ungroup() %>%
  group_by(new_level, n_tool) %>% 
  summarise(N = N[1], np = n()/N) %>%
  arrange(desc(N)) %>%
  print(n = 50) %>% filter(n_tool == 1)



core_1 <- core %>% filter(habitat %in% "human gut", cut %in% 0.5, cnt > 450) %>% 
  ungroup() %>% 
  filter(!tool %in% c("DeepARG", "RGI-DIAMOND")) %>%
  group_by(X) %>%
  summarise(n_tool = n_distinct(tool)) 


#167 unique human gut genes 

core %>% filter(habitat %in% "human gut", cut %in% 0.5, cnt > 450) %>% 
  ungroup() %>% 
  filter(tool %in% c("DeepARG", "RGI-DIAMOND"), !X %in% core_1$X) %>%
  ungroup() %>% 
  summarise(n_distinct(X))
  


only_rgi_deep_arg <- core %>% filter(habitat %in% "human gut", cut %in% 0.5, cnt > 450) %>% 
  ungroup() %>% 
  filter(tool %in% c("DeepARG", "RGI-DIAMOND"), !X %in% core_1$X) %>%
  ungroup() 

core %>% filter(habitat %in% "human gut", cut %in% 0.5, cnt > 450) %>% 
  ungroup() %>% 
  filter(tool %in% c("DeepARG", "RGI-DIAMOND"), !X %in% core_1$X) %>%
  ungroup() %>% 
  group_by(new_level) %>%
  summarise(n = n_distinct(X)) %>%
  arrange(desc(n))





fig_numbe_genes_core <- core %>% filter(habitat %in% "human gut", cut %in% 0.5, cnt > 450) %>% 
  ungroup() %>% 
  group_by(X, new_level) %>%
  summarise(n_tool = n_distinct(tool)) %>% 
  mutate(tools = ifelse(X %in% only_rgi_deep_arg$X, "Exlusive RGI/DeepARG", "Other")) %>%
  ggplot(aes (x = new_level, y = n_tool, color = tools)) +
  geom_jitter(height = 0.1) +
  ylab("Number of tools") +
  xlab("") + 
  labs(color = "") +
  scale_y_continuous(expand = c(.01, .01), 
                     limits = c(0, 10), 
                     breaks = seq(from = 0, to = 10, by = 2)) +
  scale_color_manual(values = pal_10_q[3:4]) +
  scale_x_discrete(labels = function(x) {
    x <- gsub("-", "-\n", x)
    x <- gsub(" ", "\n", x)
    x}) +
  facet_grid(. ~ new_level, scales = "free_x", space = "free_x") +
  theme_minimal() + 
  theme(legend.position = "bottom",
        panel.border = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = general_size),
        axis.text.y = element_text( size = general_size),
        strip.text = element_blank(),
        plot.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.margin = margin(0, 0, 0, 0, unit = "pt"),
        panel.spacing = unit(0, "pt"),
        legend.text = element_text(size = general_size),
        title = element_text(size = general_size + 2, face = "bold"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title = element_text(size = general_size + 1, face = "bold"),
        panel.background = element_rect(colour = "black", fill = NA))




#######






