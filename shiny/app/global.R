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


lst_results <- readRDS("data.rds")

abundance_tool_sample <- lst_results$abundance_tool_sample
core <- lst_results$core
sumpan2 <- lst_results$sumpan2
unigenes <- lst_results$unigenes
recall_fnr <- lst_results$recall_fnr
abundance_class_summary <- lst_results$abundance_class_summary
rm(lst_results)

top_abundance <- c("efflux pump", "van" , "class A beta-lactamase", 
                   "tet RPG",  "cell wall charge", "MFS efflux pump", 
                   "rpoB", "erm")

top_cso <- c("van", "efflux pump",  "tet RPG", "class A beta-lactamase", 
             "class B beta-lactamase","class C beta-lactamase", "class D beta-lactamase",
             "aph", "MFS efflux pump", "erm", "aac")



basic_tools <- c(
  "DeepARG", "fARGene","ABRicate-ARGANNOT", "ABRicate-MEGARes",
  "RGI-DIAMOND", "ABRicate-CARD","AMRFinderPlus", "ABRicate-NCBI",
  "ResFinder", "ABRicate-ResFinder")

tool_choices <- c(basic_tools,"DeepARG70","DeepARG80","DeepARG90",
                  "RGI-DIAMOND70", "RGI-DIAMOND80", "RGI-DIAMOND90")


tool_lab1 <- tools_db[match(tool_choices, names(tools_labels))]
tool_lab2 <- as.vector(tools_labels[tool_choices])
tool_lab1[!grepl("ABRicate", tool_lab2)] <- ""
tool_lab1 <- gsub("-\n","",tool_lab1)

tool_choices_label <- paste(tool_lab2, tool_lab1)

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
  summarise(n = n()) %>% 
  mutate(p = n / sum(n)) %>%
  ungroup() %>%
  mutate(new_level = gsub(" beta-lactamase","", new_level)) %>%
  mutate(new_level = gsub("rifampin inactivation enzyme","RIF-inact. enz.", new_level)) %>%
  mutate(new_level = gsub("inactivation enzyme","\ninact. enz.", new_level)) %>%
  mutate(new_level = gsub("cell wall ","cell wall\n", new_level)) %>%
  mutate(new_level = gsub("MFS efflux pump","MFS efflux", new_level)) %>%
  mutate(new_level = gsub("efflux pump","efflux", new_level)) %>%
  mutate(new_level = gsub("beta-lactam modulation resistance","beta-lactam\nmod.", new_level)) %>%
  mutate(new_level = gsub("target-modifying enzyme","target-\nmodif. enz.", new_level)) %>%
  mutate(new_level = gsub("self-resistance","self-resistance\n", new_level)) %>% 
  mutate(new_level = gsub("host-dependent","host-dependent\n", new_level)) %>% 
  mutate(new_level = gsub("variant or","variant or\n", new_level)) 




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
  mutate(texture = abundance_tool_sample$texture[match(tool_ref, abundance_tool_sample$tool)])

unigenes <- 
  unigenes %>% 
  mutate(texture = ifelse(tool %in% c("DeepARG70", "RGI-DIAMOND70"), "y70",
                          ifelse(tool %in% c("DeepARG80", "RGI-DIAMOND80"), "y80",
                                 ifelse(tool %in% c("DeepARG90", "RGI-DIAMOND90"), "y90", 
                                        texture))))



g_legend <- function(a.gplot){
  tmp <- ggplotGrob(a.gplot)
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

