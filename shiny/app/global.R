library(shiny)
library(bslib)
library(tidyverse)
library(shinyWidgets)
library(gridExtra)
library(RColorBrewer)
library(ggpattern)
library(grid)
library(reactable)
library(Cairo)
library(ggalluvial)
library(cowplot)
library(shinycssloaders)
library(ggrastr)
library(scales)
library(ragg)
library(patchwork)


APP_DIR  <- normalizePath(getwd(), mustWork = TRUE)                
ROOT_DIR <- normalizePath(file.path(APP_DIR, "..", ".."), mustWork = TRUE)

options(dplyr.summarise.inform = FALSE)
options(shiny.useragg = TRUE)


general_size <- 12
pan_core_size <- 16

pal_7 <- brewer.pal(8, "Dark2")
pal_7 <- pal_7[-7]
pal_7 <- pal_7[c(1,2,3,4,6,5,7)]
pal_10_q <- pal_7[c(1,2,3,4,5,5,6,6,7,7)]


pattern_density <- 0.001 
pattern_spacing <- 0.025
pattern_fill <- "black"
pattern_size <- 0.12

EN <- c("human gut", "human oral",  "human skin", 
        "human nose", "human vagina", 
        "dog gut", "cat gut", "mouse gut", 
        "pig gut", "wastewater", "marine", 
        "freshwater", "soil" )

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

pal_10_q <- c(
  pal_10_q,
  setNames(rep(pal_10_q["DeepARG"],   3), c("DeepARG70","DeepARG80","DeepARG90")),
  setNames(rep(pal_10_q["RGI-DIAMOND"],3), c("RGI-DIAMOND70","RGI-DIAMOND80","RGI-DIAMOND90")),
  setNames(c(pal_10_q["DeepARG"], pal_10_q["RGI-DIAMOND"], pal_10_q["RGI-DIAMOND"],
             pal_10_q["fARGene"], pal_10_q["AMRFinderPlus"]),
           c("DeepARG-aa","RGI-BLAST","RGI-DIAMOND-aa","fARGene-aa","AMRFinderPlus-nt"))
)

add_texture <- function(df) {
  df %>% mutate(texture = case_when(
    tool %in% c("DeepARG70",   "RGI-DIAMOND70") ~ "y70",
    tool %in% c("DeepARG80",   "RGI-DIAMOND80") ~ "y80",
    tool %in% c("DeepARG90",   "RGI-DIAMOND90") ~ "y90",
    TRUE ~ texture
  ))
}

scale_pattern_shared <- scale_pattern_manual(
  values = c('no'='none','yes'='stripe','y70'='crosshatch','y80'='crosshatch','y90'='crosshatch')
)

tools_labels <- c(
  "DeepARG", "fARGene", "ABRicate-\nARGANNOT", "ABRicate-\nMEGARes",
  "RGI", "ABRicate-\nCARD", "AMRFinder-\nPlus", "ABRicate-\nNCBI",
  "ResFinder", "ABRicate-\nResFinder",
  "DeepARG-70%","DeepARG-80%","DeepARG-90%","RGI-70%","RGI-80%","RGI-90%",
  "DeepARG-aa", "RGI/nBLAST", "RGI-aa", "fARGene-aa", "AMRFinder-\nPlus-nt")

names(tools_labels) <- tools_levels
tools_labels_lookup <- tools_labels

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

tools_db_factor <- c(" ", "  ", "   ", "    ",
                     "CARD", "NCBI", "ResFinder")

tools_texture <- c("ABRicate-CARD", "ABRicate-NCBI", "ABRicate-ResFinder")


lst_results <- readRDS(file.path(ROOT_DIR, "shiny", "app", "data.rds"))

abundance_tool_sample <- lst_results$abundance_tool_sample
core <- lst_results$core
sumpan2 <- lst_results$sumpan2
unigenes <- lst_results$unigenes
csc_fnr <- lst_results$csc_fnr
abundance_class_summary <- lst_results$abundance_class_summary
sumcore <- lst_results$sumcore
JI_all_plot <- lst_results$JI_all_plot
rm(lst_results)

ARO <- read.csv(file.path(ROOT_DIR, "code_R_analysis/output_abundance_diversity_resistome/conversion_ARO_parent_new_level.csv"))

table_s1 <- abundance_tool_sample %>%
  select(sample, habitat) %>%
  distinct() %>%
  rename(Sample = sample, Habitat = habitat)

table_s2 <- ARO %>%
  rename(gene_class = new_level, ARO_Term_ID = Term_ID) %>%
  mutate(abbreviation = gene_class) %>%
  mutate(abbreviation = gsub(" beta-lactamase", "", abbreviation)) %>%
  mutate(abbreviation = gsub("rifampin inactivation enzyme", "RIF-inact. enz.", abbreviation)) %>%
  mutate(abbreviation = gsub("MFS efflux pump", "efflux", abbreviation)) %>%
  mutate(abbreviation = gsub("efflux pump", "efflux", abbreviation)) %>%
  mutate(abbreviation = gsub("beta-lactam modulation resistance", "beta-lactam mod.", abbreviation)) %>%
  mutate(abbreviation = gsub("target-modifying enzyme", "target-modif. enzyme", abbreviation)) %>%
  mutate(abbreviation = gsub("abcF", "ABC-F", abbreviation))

top_abundance <- c("efflux pump", "van" , "class A beta-lactamase", 
                   "tet RPG",  "cell wall charge",
                   "rpoB", "erm")

top_cso <- c("van", "efflux pump",  "tet RPG", "class A beta-lactamase", 
             "class B beta-lactamase","class C beta-lactamase", "class D beta-lactamase",
             "aph", "erm", "aac")


tool_choices <- c(basic_tools,"DeepARG70","DeepARG80","DeepARG90",
                  "RGI-DIAMOND70", "RGI-DIAMOND80", "RGI-DIAMOND90")

tool_lab1 <- tools_db[match(tool_choices, names(tools_labels))]
tool_lab2 <- as.vector(tools_labels[tool_choices])
tool_lab1[!grepl("ABRicate", tool_lab2)] <- ""
tool_lab1 <- gsub("-\n","",tool_lab1)
tool_choices_label <- paste0(tool_lab2, tool_lab1)
tool_choices <- as.list(setNames(tool_choices, tool_choices_label))

tool_choices_single <- setNames(
  as.list(unlist(tool_choices)),
  gsub("\n", " ", names(tool_choices))
)

gene_levels <- c(top_abundance, "other")
gene_levels <- gsub(" beta-lactamase", "", gene_levels)
gene_levels <- gsub("cell wall ", "cell\nwall\n", gene_levels)
gene_levels <- gsub("MFS efflux pump", "MFS\nefflux", gene_levels)
gene_levels <- gsub("efflux pump", "efflux", gene_levels)

gene_classes <- unigenes %>% 
  ungroup() %>% 
  group_by(query) %>% 
  slice_head(n = 1) %>% 
  ungroup() %>% 
  group_by(new_level) %>% 
  summarise( n = n()) %>% 
  arrange(desc(n)) %>%
  ungroup() %>% 
  pull(new_level)

unigenes_proportion <- unigenes %>% 
  group_by(tool, tools_labels, tools_db, new_level) %>% 
  summarise(n = n_distinct(query)) %>%      
  group_by(tool, tools_labels, tools_db) %>% 
  mutate(p = n / sum(n)) %>%
  ungroup()

gene_classes_default <- unigenes_proportion %>%
  group_by(new_level) %>%
  filter(max(p) >= 0.05) %>%
  pull(new_level) %>%
  unique()

sum_core_adjust <- function(core, cnt_subset = 900, threshold_samples = 0.5){
  return(core %>% 
           filter(cut %in% threshold_samples & cnt > cnt_subset) %>% 
           ungroup() %>% 
           group_by(new_level, tool, habitat) %>% 
           summarise(unigenes = n_distinct(X)))
}

shape_tools <- rep(21, length(tools_labels))
shape_tools[tools_levels %in% tools_texture] <- 24
shape_tools[tools_levels %in% c("DeepARG70", "RGI-DIAMOND70")] <- 22
shape_tools[tools_levels %in% c("DeepARG80", "RGI-DIAMOND80")] <- 23
shape_tools[tools_levels %in% c("DeepARG90", "RGI-DIAMOND90")] <- 25

names(shape_tools) <- tools_levels
names(pal_7) <- tools_db_factor


abundance_tool_sample <- add_texture(abundance_tool_sample)

unigenes_proportion <- 
  unigenes_proportion %>% 
  mutate(texture = abundance_tool_sample$texture[match(tool, abundance_tool_sample$tool)])

abundance_class_summary <- add_texture(abundance_class_summary)

csc_fnr <- 
  csc_fnr %>% 
  mutate(texture = abundance_tool_sample$texture[match(tool_ref, abundance_tool_sample$tool)]) %>% 
  mutate(tools_db_comp =  abundance_tool_sample$tools_db[match(tool_comp, abundance_tool_sample$tool)],
         tools_db_ref =  abundance_tool_sample$tools_db[match(tool_ref, abundance_tool_sample$tool)],
         tools_labels_ref =  abundance_tool_sample$tools_labels[match(tool_ref, abundance_tool_sample$tool)],
         tools_labels_comp =  abundance_tool_sample$tools_labels[match(tool_comp, abundance_tool_sample$tool)]) %>% 
  mutate(tool_comp2 = factor(tools_labels[tool_comp], levels = tools_labels_factor),
         tool_ref2 = factor(tools_labels[tool_ref], levels = tools_labels_factor))


unigenes <- add_texture(unigenes)

sumcore <- add_texture(sumcore)
sumcore2 <- sumcore %>%
  group_by(tool, habitat, tools_labels, tools_db, texture, cnt, cut) %>%
  summarise(core = sum(unigenes)) %>%
  mutate(tool2 = factor(
    as.vector(tools_labels_lookup)[match(as.character(tool), names(tools_labels_lookup))],
    levels = tools_labels_factor))


sumpan2 <- add_texture(sumpan2) %>%
  mutate(tool2 = factor(
    as.vector(tools_labels_lookup)[match(as.character(tool), names(tools_labels_lookup))],
    levels = tools_labels_factor))

pal_10_q_2 <- pal_10_q
names(pal_10_q_2) <- as.vector(tools_labels[names(pal_10_q)])


shape_tools_2 <- shape_tools
names(shape_tools_2) <- as.vector(tools_labels[names(shape_tools_2)])

unigenes <- unigenes %>% 
  mutate(tool2 = factor(tools_labels[tool], levels = tools_labels_factor))

abundance_tool_sample <- abundance_tool_sample %>% 
  mutate(tool2 = factor(tools_labels[tool], levels = tools_labels_factor))

habitat_n_samples <- abundance_tool_sample %>%
  group_by(habitat) %>%
  summarise(N_samples = n_distinct(sample))

lims_richness <- abundance_tool_sample %>%
  group_by(tool, habitat, tools_db, tools_labels, texture) %>%
  summarise(
    median = ifelse(quantile(richness, 0.5)  < 0, 0, quantile(richness, 0.5)),
    q25    = ifelse(quantile(richness, 0.25) < 0, 0, quantile(richness, 0.25)),
    q75    = ifelse(quantile(richness, 0.75) < 0, 0, quantile(richness, 0.75)),
    w1     = ifelse(quantile(richness, 0.25) - 1.5*IQR(richness) < 0, 0,
                    quantile(richness, 0.25) - 1.5*IQR(richness)),
    w2     = ifelse(quantile(richness, 0.75) + 1.5*IQR(richness) < 0, 0,
                    quantile(richness, 0.75) + 1.5*IQR(richness))
  ) %>%
  ungroup()

abundance_class_summary <- abundance_class_summary %>% 
  mutate(tool2 = factor(tools_labels[tool], levels = tools_labels_factor))

unigenes_proportion <- unigenes_proportion %>% 
  mutate(tool2 = factor(tools_labels[tool], levels = tools_labels_factor))

table_s3_base <- unigenes %>%
  filter(tool %in% basic_tools) %>%
  ungroup() %>%
  select(query, new_level, tool)

g_legend <- function(a.gplot){
  tmp <- ggplotGrob(a.gplot)
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

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

theme_pan_core <- theme_minimal() +
  theme(
    text          = element_text(size = pan_core_size, color = "black"),
    title         = element_text(size = pan_core_size, face = "bold"),
    axis.title    = element_text(size = pan_core_size, face = "bold"),
    axis.text.x   = element_text(size = pan_core_size, angle = 90, vjust = 0.5, hjust = 1),
    axis.text.y   = element_blank(),
    strip.text.x  = element_text(size = pan_core_size, angle = 0, vjust = 0, hjust = 0.5),
    panel.background = element_blank(),
    panel.border     = element_blank(),
    panel.spacing    = unit(0, "pt"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    plot.margin      = margin(5, 5, 5, 5, unit = "pt"),
    legend.position  = "none",
    legend.text      = element_text(size = general_size),
    legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
    legend.margin    = margin(0, 0, 0, 0, unit = "pt")
  )