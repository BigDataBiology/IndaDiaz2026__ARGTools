library(shiny)
library(bslib)
library(tidyverse)
library(shinyWidgets)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(RColorBrewer)
library(ggpattern)
library(grid)
library(reactable)
library(Cairo)
library(ggalluvial)
library(cowplot)
library(scales)
library(tidyr)
library(magrittr)
library(shinycssloaders)
library(ggrastr)
library(ragg)
library(patchwork)

# Defining the root and app directory 
APP_DIR  <- normalizePath(getwd(), mustWork = TRUE)                
ROOT_DIR <- normalizePath(file.path(APP_DIR, "..", ".."), mustWork = TRUE)

options(dplyr.summarise.inform = FALSE)
options(shiny.useragg = TRUE)


general_size <- 12
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
pal_10_q <- pal_7[c(1,2,3,4,5,5,6,6,7,7)]


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
  "DeepARG", "fARGene", "ABRicate-\nARGANNOT", "ABRicate-\nMEGARes",
  "RGI", "ABRicate-\nCARD", "AMRFinder-\nPlus", "ABRicate-\nNCBI",
  "ResFinder", "ABRicate-\nResFinder",
  "DeepARG-70%","DeepARG-80%","DeepARG-90%","RGI-70%","RGI-80%","RGI-90%",
  "DeepARG-aa", "RGI/nBLAST", "RGI-aa", "fARGene-aa", "AMRFinder-\nPlus-nt")


tools_labels2 <- c(
  "DeepARG", "fARGene", "ABRicate-\nARGANNOT", "ABRicate-\nMEGARes",
  "RGI", "ABRicate-\nCARD", "AMRFinder-\nPlus", "ABRicate-\nNCBI",
  "ResFinder", "ABRicate-\nResFinder",
  "DeepARG-70%","DeepARG-80%","DeepARG-90%","RGI-70%","RGI-80%","RGI-90%",
  "DeepARG-aa", "RGI/nBLAST", "RGI-aa", "fARGene-aa", "AMRFinder-\nPlus-nt")

names(tools_labels) <- tools_levels
names(tools_labels2) <- tools_levels

tools_labels_factor <- c(
  "DeepARG", "fARGene", "ABRicate-\nARGANNOT", "ABRicate-\nMEGARes", 
  "RGI", "ABRicate-\nCARD", "AMRFinder-\nPlus", "ABRicate-\nNCBI",
  "ResFinder", "ABRicate-\nResFinder", "DeepARG-70%","DeepARG-80%","DeepARG-90%",
  "RGI-70%","RGI-80%","RGI-90%",
  "DeepARG-aa", "RGI/nBLAST", "RGI-aa", "fARGene-aa", "AMRFinder-\nPlus-nt")

tools_labels2_factor <- c(
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

# one space for DeepARG
# two spaces for fARGene
tools_db_factor <- c(" ", "  ", "   ", "    ",
                     "CARD", "NCBI", "ResFinder")

tools_texture <- c("ABRicate-CARD", "ABRicate-NCBI", "ABRicate-ResFinder")
#lst_results <- readRDS("data.rds")

lst_results <- readRDS(file.path(ROOT_DIR, "shiny", "app", "data.rds"))

abundance_tool_sample <- lst_results$abundance_tool_sample
core <- lst_results$core
sumpan2 <- lst_results$sumpan2
unigenes <- lst_results$unigenes
recall_fnr <- lst_results$recall_fnr
abundance_class_summary <- lst_results$abundance_class_summary
sumcore <- lst_results$sumcore
rm(lst_results)

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
#tool_lab2 <- grep(" ", "", tool_lab2)


tool_lab1 <- gsub("-\n","",tool_lab1)
tool_choices_label <- paste0(tool_lab2, tool_lab1)

tool_choices <- as.list(setNames(tool_choices, tool_choices_label))

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

unigenes_propotion <- unigenes %>% 
  group_by(tool, tools_labels, tools_db, new_level) %>% 
  summarise(n = n_distinct(query)) %>%      
  group_by(tool, tools_labels, tools_db) %>% 
  mutate(p = n / sum(n))


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


abundance_tool_sample <- 
  abundance_tool_sample %>% 
  mutate(texture = ifelse(tool %in% c("DeepARG70", "RGI-DIAMOND70"), "y70",
                          ifelse(tool %in% c("DeepARG80", "RGI-DIAMOND80"), "y80",
                                 ifelse(tool %in% c("DeepARG90", "RGI-DIAMOND90"), "y90", 
                                        texture))))

unigenes_propotion <- 
  unigenes_propotion %>% 
  mutate(texture = abundance_tool_sample$texture[match(tool, abundance_tool_sample$tool)])

abundance_class_summary <- 
  abundance_class_summary %>% 
  mutate(texture = ifelse(tool %in% c("DeepARG70", "RGI-DIAMOND70"), "y70",
                          ifelse(tool %in% c("DeepARG80", "RGI-DIAMOND80"), "y80",
                                 ifelse(tool %in% c("DeepARG90", "RGI-DIAMOND90"), "y90", 
                                        texture))))

recall_fnr <- 
  recall_fnr %>% 
  mutate(texture = abundance_tool_sample$texture[match(tool_ref, abundance_tool_sample$tool)]) %>% 
  mutate(tools_db_comp =  abundance_tool_sample$tools_db[match(tool_comp, abundance_tool_sample$tool)],
         tools_db_ref =  abundance_tool_sample$tools_db[match(tool_ref, abundance_tool_sample$tool)],
         tools_labels_ref =  abundance_tool_sample$tools_labels[match(tool_ref, abundance_tool_sample$tool)],
         tools_labels_comp =  abundance_tool_sample$tools_labels[match(tool_comp, abundance_tool_sample$tool)]) %>% 
  mutate(tool_comp2 = factor(tools_labels2[tool_comp], levels = tools_labels2_factor),
         tool_ref2 = factor(tools_labels2[tool_ref], levels = tools_labels2_factor))


unigenes <- 
  unigenes %>% 
  mutate(texture = ifelse(tool %in% c("DeepARG70", "RGI-DIAMOND70"), "y70",
                          ifelse(tool %in% c("DeepARG80", "RGI-DIAMOND80"), "y80",
                                 ifelse(tool %in% c("DeepARG90", "RGI-DIAMOND90"), "y90", 
                                        texture))))

sumcore <- sumcore %>% 
  mutate(texture = ifelse(tool %in% c("DeepARG70", "RGI-DIAMOND70"), "y70",
                          ifelse(tool %in% c("DeepARG80", "RGI-DIAMOND80"), "y80",
                                 ifelse(tool %in% c("DeepARG90", "RGI-DIAMOND90"), "y90", 
                                        texture))))

sumcore2 <- sumcore %>% 
  group_by(tool, habitat, tools_labels, tools_db, texture, cnt, cut) %>% 
  summarise(core = sum(unigenes)) %>% 
  mutate(tool2 = factor(tools_labels2[tool], levels = tools_labels2_factor))

sumpan2 <- sumpan2 %>% 
  mutate(texture = ifelse(tool %in% c("DeepARG70", "RGI-DIAMOND70"), "y70",
                          ifelse(tool %in% c("DeepARG80", "RGI-DIAMOND80"), "y80",
                                 ifelse(tool %in% c("DeepARG90", "RGI-DIAMOND90"), "y90", 
                                        texture)))) %>% 
  mutate(tool2 = factor(tools_labels2[tool], levels = tools_labels2_factor))

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

pal_10_q
pal_10_q_2 <- pal_10_q
names(pal_10_q_2) <- as.vector(tools_labels2[names(pal_10_q)])


shape_tools_2 <- shape_tools
names(shape_tools_2) <- as.vector(tools_labels2[names(shape_tools_2)])

unigenes <- unigenes %>% 
  mutate(tool2 = factor(tools_labels2[tool], levels = tools_labels2_factor))

abundance_tool_sample <- abundance_tool_sample %>% 
  mutate(tool2 = factor(tools_labels2[tool], levels = tools_labels2_factor))

abundance_class_summary <- abundance_class_summary %>% 
  mutate(tool2 = factor(tools_labels2[tool], levels = tools_labels2_factor))

unigenes_propotion <- unigenes_propotion %>% 
  mutate(tool2 = factor(tools_labels2[tool], levels = tools_labels2_factor))