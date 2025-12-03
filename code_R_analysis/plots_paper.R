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

setwd("~/Documents/GitHub/arg_compare/code_R_analysis") 
options(dplyr.summarise.inform = FALSE)

source("helper.R")

# Sourced gene classes 
# gene_classes_list is a vector for labeling the plots 
# gene_classes is the list with the names from the prepared datasets and used for shiny


# FORMAT PLOTS OLD COLORS
pal_7 <- brewer.pal(8, "BrBG")
pal_7 <- pal_7[c(4,5,3,6,2,7,1)]
pal_10_q <- pal_7[c(1,2,3,4,5,5,6,6,7,7)]

pal_10_complete <- brewer.pal(10, "BrBG")
pal_10_complete <- pal_10_complete[c(-1,-10)]

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
tools_texture <- c("ABRicate-CARD", "ABRicate-NCBI", "ABRicate-ResFinder")

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
         location = metadata$location[match(sample, metadata$sample_id)]) %>% # not so dplyr way to fetch habitat and environment


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

JI_core_class <- return_overlap_tools(core %>% 
  filter(cut %in% 0.5 & cnt > 450, habitat %in% "human gut") %>% 
  rename(query = X))


top20 <- c("van", "efflux pump", "cell wall charge", "rpoB", "tet RPG", "class A beta-lactamase", "class B beta-lactamase", "class C beta-lactamase", "class D beta-lactamase", 
           "aac", "aph", "MFS efflux pump", "erm", "target-modifying enzyme", "tet enzyme", "mph", "fos", "cat", "abcF")


pan_core_env <- c("human gut", "human skin", "pig gut", "dog gut", "wastewater",  "marine", "freshwater", "soil")

list_plot_core_class <- list()
list_plot_pan_class <- list()

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



################################################################################################
# Fig 1

square_legends <- g_legend(
    plot_count_genes_tool_legends(
      unigenes = unigenes , tools_for_figure = tools_levels, 
      general_size = general_size, # font size
      pal_10_q = pal_10_q, # pallet
      tool_label = tools_labels,  # labels for the tools
      tools_levels = tools_levels) +  # factor names of the tools
    ggtitle("A") + theme( title = element_text(size = general_size + 2, face = "bold")))

p_counts <- plot_count_genes_tool(
  unigenes = unigenes, 
  tools_for_figure = tools_levels, 
  general_size = general_size, 
  pal_10_q = pal_10_q, 
  tool_label = tools_labels , 
  tools_levels = tools_levels, 
  texture = tools_texture) # which tools have stripes

square_texture_legends <- g_legend( p_counts + guides(pattern = "none"))

p1a <- p_counts + 
  ggtitle("A") + theme(legend.position = "none", title = element_text(size = general_size + 2, face = "bold"))

p1b <- return_heatmap_overalp(
  JI_all = JI_all, # overlap between tools
  tools_selected = as.character(unique(JI_all$tool_ref)), # tools selected to be plotted 
  general_size = general_size, # font size
  tool_label = tools_labels, # tool names in the plots
  tools_levels = unique(as.character(JI_all$tool_ref))) + # factor names 
  ggtitle("B") + theme(title = element_text(size = general_size + 2, face = "bold"))


p2a <- plot_total_abundance_diversity_new_version(
  dataset = abundance_tool_sample, # 
  tools_labels = tools_labels,  #
  tools_to_plot = tools_levels,  #
  environments_plot = h2, # habitats to plot (aggregated humans and mammals)
  general_size = general_size, # font size
  pal_10_q = pal_10_q , # pallet
  metric = "abundance", # metric (abundance or diversity)
  sd = 2025, # seed to plot random samples in the distribution 
  obs = 200,  # number of samples to plot as dots per environment
  texture = tools_texture) # texture for repeated color 

p2b <- plot_total_abundance_diversity_new_version(
  abundance_tool_sample, 
  tools_labels, 
  tools_levels, 
  h2, 
  general_size, 
  pal_10_q , 
  "diversity", 
  sd = 2025, 
  obs = 200, 
  tools_texture)

################################################################################################
# Fig 2
  

circle_legends <- g_legend(pan_core %>% 
  ggplot(aes(x = md, y = core +1, fill = tool)) +
    geom_point(size = 3, color = "black", shape = 21) +
    scale_fill_manual(values = pal_10_q, labels = tools_labels) + 
    labs(fill = "") +
    theme(legend.position = "bottom",
      legend.text = element_text(size = general_size)))

triangle_legends <- g_legend(pan_core %>% 
  ggplot(aes(x = md, y = core +1, fill = tool)) +
    geom_point(size = 3, color = "black", shape = 24) +
    scale_fill_manual(values = pal_10_q, labels = tools_labels) + 
    labs(fill = "") +
    theme(legend.position = "bottom",
      legend.text = element_text(size = general_size)))


p_pan_core <- dot_plot_pan_core(
    pan_core, 
    pan_core_env, 
    tools_labels, 
    tools_levels, 
    tools_levels, 
    tools_texture, 
    general_size, 
    pal_10_q) + 
  ggtitle("A") + theme(legend.position = "none")

p_pan_core_supp <- dot_plot_pan_core(
  pan_core, 
  EN2[!EN2 %in% pan_core_env], 
  tools_labels, 
  tools_levels, 
  tools_levels, 
  tools_texture, 
  general_size, 
  pal_10_q) + 
  ggtitle("A") + theme(legend.position = "none")



# PLOTS


for( e in EN2){
  plot_list <- plot_pan_core_class(
    pan = pan, # pan resistome data
    sumcore = sumcore,  # core resistome data
    top20 = top20, # gene classes to plot, the rest go to Other
    tools_labels = tools_labels, #
    tools_to_plot = tools_levels, #
    tools_levels = tools_levels, #
    environments_plot = e,  #
    general_size = general_size, 
    pal_10_q = pal_10_q, 
    tools_texture = tools_texture)  
  
  list_plot_core_class[[e]] <- plot_list[[1]]
  list_plot_pan_class[[e]] <- plot_list[[2]]
}

p_alluvial  <- plot_alluvial_classes(unigenes, levels_unigenes, 0.95, 0.01, tools_levels, tools_labels, tools_levels, pal_10_q, general_size, gene_classes_list)


p42 <- plot_recall_fnr_texture(recall_fnr, tools_levels, unique(unigenes$new_level), tools_labels, pal_10_q, general_size, tools_levels, tools_texture)

p5b10 <- one_tool_recall_fnr_id(unigenes, recall_fnr, "DeepARG", pal_10_q[3:5], general_size, tool_label, tools_levels_2) + 
  ggtitle("DeepARG") + 
  theme(legend.position = "none")

p5b2 <- one_tool_recall_fnr_id(unigenes, recall_fnr, "RGI-DIAMOND", pal_10_q[3:5], general_size, tool_label, tools_levels_2) +
  ggtitle("RGI-CARD") + 
  theme(legend.position = "none")

p5b3 <- one_tool_recall_fnr_id(unigenes, recall_fnr, "fARGene", pal_10_q[3:5], general_size, tool_label, tools_levels_2) +
  ggtitle("fARGene") + 
  theme(legend.position = "none")

p5b4 <- one_tool_recall_fnr_id(unigenes, recall_fnr, "AMRFinderPlus", pal_10_q[3:5], general_size, tool_label, tools_levels_2) +
  ggtitle("AMRFinderPlus-NCBI") + 
  theme(legend.position = "none")

p5 <- grid.arrange(p5a, p5b1, p5b2, g_legend(p5b10), layout_matrix = rbind(c(1, 1), c(1, 1), c(1, 1),
                                                                           c(2, 3), c(2, 3), c(2, 3), c(2, 3), c(2, 3), c(2, 3), c(2, 3), c(2, 3),
                                                                           c(4,4)))


ggsave("output_plots/fig1c.svg", p1a, width = 60, height = 75, unit = "mm")
ggsave("output_plots/fig1c_legends_1.svg", p1a_legends_1, width = 120, height = 75, unit = "mm")
ggsave("output_plots/fig1c_legends.svg", p1a_legends, width = 120, height = 75, unit = "mm")
ggsave("output_plots/figS19.svg", p1b + ggtitle(""), width = 90, height = 75, unit = "mm")
ggsave("output_plots/fig2.svg", p2a, width = 180, height = 120, unit = "mm")
ggsave("output_plots/figS1.svg", p2b, width = 180, height = 120, unit = "mm")
ggsave("output_plots/fig3.svg", p2 + theme(legend.position = "none"), width = 180, height = 120, unit = "mm")
ggsave("output_plots/fig3_supp.svg", p2_other + theme(legend.position = "none"), width = 180, height = 120, unit = "mm")

for( e in EN2){
  ggsave(paste0("output_plots/pan_core_class_", e,".svg"), 
         grid.arrange(list_plot_core_class[[e]], list_plot_pan_class[[e]], nrow = 2) , 
         width = 180, height = 180, unit = "mm")
}

ggsave("output_plots/alluvial.svg", p_alluvial , width = 180, height = 180, unit = "mm")
ggsave("output_plots/fig5.svg", p42, width = 180, height = 80, unit = "mm")


ggsave("output_plots/fig6.svg", p5, width = 180, height = 220, unit = "mm")

################################################################################################   
  # Fig 3

# abundance_env_class <- abundance_medians_env_class(abundance_class, "human gut", tools_levels, unique(unigenes$new_level))


## 

# abundance_env_class <- abundance_medians_env_class(abundance_class, "human gut", tools_levels, top20)
# unigenes_class <- get_unigenes_class(unigenes, tools_levels,  top20) %>% 
#   ungroup() %>% 
#   group_by(tool) %>% 
#   mutate(proportion = n / sum(n)) %>% 
#   arrange(desc(proportion)) %>% 
#   mutate(cum_p = cumsum(proportion)) %>% 
#   ungroup() 
# 
# p3a1 <- plot_unigenes_env_tool_class(unigenes_class,  general_size, pal_10_q, tools_levels, tools_levels, tools_levels, 0, 15000, 5, T) + 
#   ggtitle("") + theme(title = element_text(size = general_size + 2, face = "bold"))
# 
# p3a2 <- plot_unigenes_env_tool_class(unigenes_class,  general_size, pal_10_q, tools_levels, tools_levels, tools_levels, 15000, max(unigenes_class$n), 2, T) + 
#   ggtitle("A") + ylab("") + theme(title = element_text(size = general_size + 2, face = "bold"), axis.text.x  = element_blank())
# 
# p3a <- grid.arrange(p3a2 + theme(plot.margin = margin(0, 0, -10, 0)), 
#                     p3a1 + theme(plot.margin = margin(-10, 0, 0, 0)), 
#                     layout_matrix = rbind(c(1, 1), c(2, 2), c(2, 2), c(2, 2)))
#   
# p3b0 <- plot_abundance_env_tool_class(abundance_env_class, tools_levels, "human gut", top20, general_size, pal_10_q, tools_levels, tools_levels) +
#   ggtitle("B") + theme(title = element_text(size = general_size + 2, face = "bold"), strip.text = element_blank())
# 
# 
# p3b <- p3b0 + theme(legend.position = "none")
# p3c <- plot_core_env_tool_class(sumcore, tools_levels, "human gut", top20, general_size, pal_10_q, tools_labels, tools_levels) +
#   ggtitle("C") + theme(legend.position = "none", title = element_text(size = general_size + 2, face = "bold"), strip.text = element_blank()) +
#   ylab("Number of ARGs")
# 
# p3 <- grid.arrange(p3a, p3b, p3c, g_legend(p3b0), 
#                    layout_matrix = rbind(c(1, 1), c(1, 1), c(1, 1), c(1, 1), c(1, 1), 
#                                          c(2, 2), c(2, 2), c(2, 2), c(2, 2), c(2, 2), 
#                                          c(3, 3), c(3, 3), c(3, 3), c(3, 3), c(3, 3), c(4,4)))

# ggsave("~/Documents/plots_project6/fig4.svg", p3, width = 180, height = 220, unit = "mm")


# for(microbiome in 1:length(microbiomes)){
#   abundance_env_class <- abundance_medians_env_class(abundance_class, microbiomes[microbiome], tools_levels_2, top20)
#   
#   p3b0 <- plot_abundance_env_tool_class(abundance_env_class, tools_levels_2, microbiomes[microbiome], top20, general_size, pal_10_q, tools_levels_2, tools_levels_2) +
#     ggtitle(paste0(microbiomes[microbiome], "\nA")) + theme(title = element_text(size = general_size + 2, face = "bold"), strip.text = element_blank())
#   
#   p3b <- p3b0 + theme(legend.position = "none")
#   p3c <- plot_core_env_tool_class(sumcore, tools_levels_2, microbiomes[microbiome], top20, general_size, pal_10_q, tools_levels_2, tools_levels_2) +
#     ggtitle("B") + theme(legend.position = "none", title = element_text(size = general_size + 2, face = "bold"), strip.text = element_blank()) +
#     ylab("Number of ARGs")
#   
#   p3 <- grid.arrange(p3b, p3c, g_legend(p3b0), layout_matrix = rbind(c(1, 1), c(1, 1), c(1, 1), c(1, 1), c(1, 1), 
#                                                                      c(2, 2), c(2, 2), c(2, 2), c(2, 2), c(2, 2),
#                                                                      c(3,3)))
#   ggsave(paste0(paste0("~/Documents/plots_project6/figS", 6 + microbiome),".svg"), p3, width = 180, height = 150, unit = "mm")
# }




  ################################################################################################
  # Fig 4
  


########################################################################################




df0 <- unigenes %>% ungroup() %>%
  mutate(tool_2 = ifelse(tool %in% c("ABRicate-ARGANNOT", "ABRicate-MEGARes", "ABRicate-CARD",
                                     "ABRicate-ResFinder", "ABRicate-NCBI", "ResFinder"), "ABRicate and ResFinder", as.character(tool))) %>% 
  filter(!((tool %in% "AMRFinderPlus") & (query %in% lst$amrfinder.norm.prot$query[lst$amrfinder.norm.prot$Method %in% "HMM"]))) %>% 
  filter(!tool %in% "fARGene") %>% 
  mutate(bin = cut_width(round(id), width = 1, boundary = min(round(id), na.rm = T))) %>%
  group_by(tool_2, bin) %>%
  summarise(n = n(), .groups = "drop_last") %>%
  mutate(frac = n / sum(n)) %>%
  ungroup() %>%
  mutate(tool_2 = factor(tool_2, levels = c("DeepARG", "RGI-DIAMOND", "AMRFinderPlus", "ABRicate and ResFinder"))) %>%
  rename(tool = tool_2)


df_min <- unigenes %>% ungroup() %>%
  mutate(tool_2 = ifelse(tool %in% c("ABRicate-ARGANNOT", "ABRicate-MEGARes", "ABRicate-CARD",
                                     "ABRicate-ResFinder", "ABRicate-NCBI", "ResFinder"), "ABRicate and ResFinder", as.character(tool))) %>% 
  filter(!((tool %in% "AMRFinderPlus") & (query %in% lst$amrfinder.norm.prot$query[lst$amrfinder.norm.prot$Method %in% "HMM"]))) %>% 
  filter(!tool %in% "fARGene") %>% 
  select(-tool) %>% 
  rename(tool = tool_2) %>% 
  group_by(tool) %>% summarise(mi = min(id, na.rm = T), ma = max(id, na.rm = T))

df0 <- df0 %>% left_join(df_min) %>% 
  mutate(tool = factor(tool, levels = c("DeepARG", "RGI-DIAMOND", "AMRFinderPlus", "ABRicate and ResFinder"))) 


df0 <- df0 %>% mutate(texture = ifelse(tool %in% tools_texture, "yes", "no")) %>%
  ungroup() %>%
  mutate(id2 =  as.numeric(sub(".*,(.*)\\]", "\\1", bin))) %>%
  mutate(id1 =  as.numeric(sub("^[\\(\\[]([0-9]+),.*\\]$", "\\1", bin))) %>% 
  group_by(tool) %>%
  mutate(xmin = ifelse(id1 < mi, mi, id1)) %>%
  mutate(xmax = ifelse(id2 > ma, ma, id2))

id_plot <- df0 %>%
  ggplot(aes(x = id2 , color = tool, y = frac, fill = tool, pattern = texture)) +
  geom_rect_pattern(aes(xmin = xmin, xmax = xmax, 
                ymin = 0, ymax = frac), linewidth = 0,
                pattern_color = "white",
                pattern_density = 0.1, 
                pattern_spacing = 0.2, 
                pattern_fill = "white",
                pattern_key_scale_factor = 0.05) + 
  scale_pattern_manual(values = c('no' = 'none', 'yes' = 'stripe')) +
  theme_minimal() + 
  scale_fill_manual(values = pal_10_q[c(1,5,7,9)], labels = tools_labels[c(1,5,7,9)], guide = guide_legend(nrow = 2)) +
  scale_color_manual(values = pal_10_q[c(1,5,7,9)], labels = tools_labels[c(1,5,7,9)], guide = guide_legend(nrow = 2)) +
  #scale_x_continuous(limits = c(20, 100), breaks = seq(from = 0, to = 100, by = 20)) + 
  ggtitle("") + 
  ylab("Proportion of genes") + 
  xlab("% Identity level") + 
  labs(color = "", fill = "") +
  facet_grid(. ~ tool , scales = "free") +
  theme(legend.position = "none",
        panel.border = element_blank(),
        axis.text.x = element_text(angle = 0, hjust = 1, vjust = 0.5, size = general_size),
        axis.text.y = element_blank(),
        strip.text = element_text(size = general_size, face = "bold"),
        plot.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.margin = margin(0, 0, 0, 0, unit = "pt"),
        panel.spacing = unit(0, "pt"),
        legend.text = element_text(size = general_size),
        title = element_text(size = general_size + 2, face = "bold"),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = general_size + 1, face = "bold"),
        panel.background = element_rect(colour = "black", fill = NA))

id_plot



ggsave("output_plots/ids.svg", id_plot,  width = 180, height = 220, unit = "mm")


################################################################################################
# Fig S1
# 
# rgi.ven = ggVennDiagram(list("RGI-DIAMOND" = lst$rgi.diamond$query,
#                              "RGI-BLAST" = lst$rgi.blast$query,
#                              "RGI-DIAMOND (aa)" = lst$rgi.diamond.prot$query),
#                         color = 1, lwd = 0.7, label_size = 2, set_size = 2) + 
#   scale_fill_gradient(low = "#F6E8C3", high = pal_10_q[5]) +
#   theme(legend.position = "none") +
#   ggtitle("A") +
#   theme(
#     legend.text = element_text(size = general_size),
#     legend.title = element_text(size = general_size + 1),
#     text = element_text(size = general_size),
#     axis.text  = element_blank(),
#     title = element_text(size = general_size + 2, , face = "bold"))
# 
# rgi.ven
# 
# 
# deeparg.ven = ggVennDiagram(list("DeepARG" = lst$deeparg.norm$query,
#                                  "DeepARG (aa)" = lst$deeparg.norm.prot$query),
#                             label_percent_digit = 1,
#                             color = 1, lwd = 0.1, label_size = 2, set_size = 2) + 
#   scale_fill_gradient(low = "#F6E8C3", high = pal_10_q[5]) +
#   theme(legend.position = "none") +
#   scale_x_continuous(labels = ) +
#   ggtitle("B") +
#   theme(
#     legend.text = element_text(size = general_size),
#     legend.title = element_text(size = general_size + 1),
#     text = element_text(size = general_size),
#     axis.text  = element_blank(),
#     title = element_text(size = general_size + 2, , face = "bold"))
# 
# deeparg.ven
# 
# 
# fargene.ven = ggVennDiagram(list("fARGene" = lst$fargene$query,
#                                  "fARGene (aa)" = lst$fargene.prot$query),
#                             label_percent_digit = 2,
#                             color = 1, lwd = 0.1, label_size = 2, set_size = 2) + 
#   scale_fill_gradient(low = "#F6E8C3", high = pal_10_q[5]) +  
#   theme(legend.position = "none") +
#   ggtitle("C") +
#   theme(
#     legend.text = element_text(size = general_size),
#     legend.title = element_text(size = general_size + 1),
#     text = element_text(size = general_size),
#     axis.text  = element_blank(),
#     title = element_text(size = general_size + 2, , face = "bold"))
# 
# fargene.ven
# 
# amrfinder.ven = ggVennDiagram(list("AMRFinderPlus" = lst$amrfinder.norm$query,
#                                    "AMRFinderPlus (aa)" = lst$amrfinder.norm.prot$query),
#                               label_percent_digit = 1,
#                               color = 1, lwd = 0.1, label_size = 2, set_size = 2) + 
#   scale_fill_gradient(low = "#F6E8C3", high = pal_10_q[5]) +  
#   theme(legend.position = "none") +
#   ggtitle("D") +
#   theme(
#     legend.text = element_text(size = general_size),
#     legend.title = element_text(size = general_size + 1),
#     text = element_text(size = general_size),
#     axis.text  = element_blank(),
#     title = element_text(size = general_size + 2, , face = "bold"))
# 
# amrfinder.ven

theme_overlap <- theme_minimal() + 
  theme(legend.position = "none",
        panel.border = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = general_size),
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

id_plot

rgi_query <- data.frame(gene = unique(c(lst$rgi.blast$query, lst$rgi.diamond$query, lst$rgi.diamond.prot$query)))
rgi_query <- rgi_query %>% mutate(dataset = ifelse(gene %in% intersect(intersect(lst$rgi.blast$query, lst$rgi.diamond$query), lst$rgi.diamond.prot$query), "All",
                                            ifelse(gene %in% intersect(lst$rgi.blast$query, lst$rgi.diamond$query), "BLAST nt and DIAMOND nt",
                                            ifelse(gene %in% intersect(lst$rgi.diamond.prot$query, lst$rgi.diamond$query), "DIAMOND nt and DIAMOND aa",
                                            ifelse(gene %in% intersect(lst$rgi.diamond.prot$query, lst$rgi.blast$query), "BLAST nt and DIAMOND aa",
                                            ifelse(gene %in% lst$rgi.diamond.prot$query, "DIAMOND aa",
                                            ifelse(gene %in% lst$rgi.diamond$query, "DIAMOND nt",
                                            ifelse(gene %in% lst$rgi.blast$query, "BLAST nt", NA)))))))) %>%
  mutate(dataset = factor(dataset, levels = c("All", "BLAST nt and DIAMOND nt", "DIAMOND nt and DIAMOND aa", "BLAST nt and DIAMOND aa", "DIAMOND aa", "DIAMOND nt", "BLAST nt")))
                                            
rgi_over <- rgi_query %>% 
  group_by(dataset) %>% summarise(n = n()) %>% mutate(p = n/sum(n)) %>%
  ggplot(aes(x = dataset, y = n, fill = dataset)) + 
  geom_col() + 
  scale_fill_manual(values = pal_10_q)  +
  ylab("ARGs") + 
  theme_overlap


deeparg_query <- data.frame(gene = unique(c(lst$deeparg.norm$query, lst$deeparg.norm.prot$query)))
deeparg_query <- deeparg_query %>% mutate(dataset = ifelse(gene %in% intersect(lst$deeparg.norm$query, lst$deeparg.norm.prot$query), "Both",
                                                   ifelse(gene %in% lst$deeparg.norm$query, "DeepARG nt", "DeepARG aa"))) %>% 
  mutate(dataset = factor(dataset, levels = c("Both", "DeepARG nt", "DeepARG aa")))


deeparg_over <- deeparg_query %>% 
  group_by(dataset) %>% summarise(n = n()) %>% mutate(p = n/sum(n)) %>%
  ggplot(aes(x = dataset, y = n, fill = dataset)) + 
  geom_col() + 
  scale_fill_manual(values = pal_10_q)  +
  ylab("ARGs") + 
  theme_overlap


fargene_query <- data.frame(gene = unique(c(lst$fargene$query, lst$fargene.prot$query)))
fargene_query <- fargene_query %>% mutate(dataset = ifelse(gene %in% intersect(lst$fargene$query, lst$fargene.prot$query), "Both",
                                                           ifelse(gene %in% lst$fargene$query, "fARGene nt", "fARGene aa"))) %>% 
  mutate(dataset = factor(dataset, levels = c("Both", "fARGene nt", "fARGene aa")))


fargene_over <- fargene_query %>% 
  group_by(dataset) %>% summarise(n = n()) %>% mutate(p = n/sum(n)) %>%
  ggplot(aes(x = dataset, y = n, fill = dataset)) + 
  geom_col() + 
  scale_fill_manual(values = pal_10_q)  +
  ylab("ARGs") + 
  theme_overlap



amrfinder_query <- data.frame(gene = unique(c(lst$amrfinder.norm$query, lst$amrfinder.norm.prot$query)))
amrfinder_query <- amrfinder_query %>% mutate(dataset = ifelse(gene %in% intersect(lst$amrfinder.norm$query, lst$amrfinder.norm.prot$query), "Both",
                                                           ifelse(gene %in% lst$amrfinder.norm$query, "AMRFinderPlus nt", "AMRFinderPlus aa"))) %>% 
  mutate(dataset = factor(dataset, levels = c("Both", "AMRFinderPlus nt", "AMRFinderPlus aa")))


amrfinder_over <-amrfinder_query %>% 
  group_by(dataset) %>% summarise(n = n()) %>% mutate(p = n/sum(n)) %>%
  ggplot(aes(x = dataset, y = n, fill = dataset)) + 
  geom_col() + 
  scale_fill_manual(values = pal_10_q)  +
  ylab("ARGs") + 
  theme_overlap





################################################################################################
# Fig S2




# ggsave("~/Documents/plots_project6/figS3.svg", ps2, width = 180, height = 150, unit = "mm")
# ggsave("~/Documents/plots_project6/figS4.svg", ps2, width = 180, height = 150, unit = "mm")
# ggsave("~/Documents/plots_project6/figS5.svg", ps2, width = 180, height = 150, unit = "mm")


################################################################################################
# Fig S3


# class overlap for core resistome 


# data.frame(JI_core_class %>% ungroup() %>% group_by(new_level) %>% summarise(m = sum(jaccard), n = n()) %>% mutate(m = ifelse(is.na(m), 0, m / 90)) %>% arrange(desc(m)))
# length(unique(paste(JI_core_class$tool_ref, JI_core_class$tool_comp)))
# 
# core_jaccard_human <- JI_core_class %>% ungroup() %>% 
#   filter(as.numeric(tool_comp) < as.numeric(tool_ref)) %>% #complete(tool, new_level, fill = list(unigenes = 0, total = 0)) %>%
#   ungroup() %>% 
#   filter(as.numeric(tool_ref) != as.numeric(tool_comp)) %>% 
#   ggplot(aes(x = tool_ref, y = tool_comp, fill = jaccard)) +
#   geom_tile(na.rm = T) +
#   scale_x_discrete(labels = tools_levels) +
#   scale_y_discrete(labels = tools_levels) +
#   facet_wrap(new_level ~ .) +
#   scale_fill_gradientn(colors = RColorBrewer::brewer.pal(9, "YlOrBr"), na.value = "white") +
#   theme_minimal() +
#   labs(fill = "") +
#   xlab("") +
#   ylab("") +
#   ggtitle("") + 
#   theme(
#     legend.position = "bottom",
#     legend.text = element_text(size = general_size ),
#     panel.border = element_rect(fill = "transparent", color = "black", linewidth = 1),
#     #panel.grid.major.x = element_blank(),
#     panel.grid.minor.x = element_blank(),
#     plot.margin = margin(0, 0, 0, 0, unit = "pt"),
#     legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
#     legend.margin = margin(0, 0, 0, 0, unit = "pt"),
#     panel.spacing = unit(0, "pt"),
#     strip.text = element_text(size = general_size + 1, face = "bold"),
#     axis.title = element_text(size = general_size + 1, face = "bold"),
#     axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = general_size),
#     axis.text.y = element_text(size = general_size))
# 
# core_jaccard_human
# 
# ggsave("~/Documents/plots_project6/figS6.svg", core_jaccard_human, width = 300, height = 300, unit = "mm") 


################################################################################################
# Fig S4

##  Tables 

# sum(table(abundance$habitat[!duplicated(abundance$sample)]))
# table(abundance$habitat[!duplicated(abundance$sample)])

# write.csv(abundance %>% ungroup() %>% select(sample, habitat) %>% filter(!habitat %in% c("amplicon","built-environment","isolate")) %>%
#  rename(accession = sample) %>% group_by(accession) %>% slice_head(n = 1), "~/Documents/plots_project3/accession_metagenomes.csv", row.names = FALSE)

# write.csv(abundance %>% ungroup() %>% select(sample, habitat) %>% filter(!habitat %in% c("amplicon","built-environment","isolate")) %>%
#             rename(accession = sample) %>% group_by(accession) %>% slice_head(n = 1) %>% 
#            ungroup() %>% group_by(habitat) %>% summarise(samples = n()), "~/Documents/plots_project3/summary_metagenomes.csv", row.names = FALSE)




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


