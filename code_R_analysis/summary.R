library(dplyr)
library(ggplot2)
library(ggVennDiagram)
library(gridExtra)
library(tidyverse)
library(RColorBrewer)

setwd("~/Documents/GitHub/arg_compare/")
df2 <- readRDS(file = "code_R_analysis/output_abundance_diversity_resistome/conversion_ARO_parent_new_level.rds")
lst <- readRDS("code_R_analysis/output_abundance_diversity_resistome/results_tools.rds")
# args_abundances <- read.delim("data/abundances/args_abundances.tsv") 
# total numbrer of unigenes captured with abundance (for the ok environments 188 441 genes 
# metadata <- read.delim("data/metadata_GMGC10.sample.meta.tsv")
# args_abundances <- args_abundances %>% mutate(habitat = metadata$habitat[match(sample, metadata$sample_id)])


lst$rgi.blast <- lst$rgi.blast %>% mutate(tool =  "RGI (BLAST - nt)")
lst$rgi.diamond <- lst$rgi.diamond %>% mutate(tool =  "RGI (DIAMOND - nt)")
lst$rgi.diamond.prot <- lst$rgi.diamond.prot %>% mutate(tool =  "RGI (DIAMOND - aa)")


abundance <- readRDS("code_R_analysis/output_abundance_diversity_resistome/abundance_diversity.rds")

# d is divergent / q is qualitative
pal_8_q <- brewer.pal(8, "Dark2")
pal_10_d <- brewer.pal(10, "BrBG")
pal_4_d <- brewer.pal(4, "BrBG")
pal_11_d <- brewer.pal(11, "BrBG")
pal_9_d <- brewer.pal(9, "BrBG")
# pallet for 15 tool + format repeated by tool
pal15_rep <- c(rep(pal_10_d[1], 2), rep(pal_10_d[2], 3), rep(pal_10_d[3], 2),
               pal_10_d[4], rep(pal_10_d[5], 2), pal_10_d[6:10])
# size for font in plots 
general_size = 10


# HABITATS
EN <- c("human gut", "human oral",  "human skin", "human nose", "human vagina", 
        "dog gut", "cat gut", "mouse gut", "pig gut", "wastewater", "marine", 
        "freshwater", "soil" , "amplicon", "isolate",  "built-environment" )

# SOURCE FOR EACH HABITAT
SO <- c(rep("Humans", 5), rep("Mammals", 4),  
        "Wastewater", "Marine", "Freshwater", 
        "Soil", rep("Other", 3))

names(SO) <- EN

tools_levels <- c("DeepARG (nt)", "DeepARG (aa)", "RGI (BLAST - nt)", "RGI (BLAST - aa)", 
                  "RGI (DIAMOND - nt)", "RGI (DIAMOND - aa)", "fARGene (nt)", 
                  "fARGene (aa)", "ResFinder (nt)", "AMRFinderPlus (nt)", "AMRFinderPlus (aa)", 
                  "ABRicate (ARG-ANNOT - nt)", "ABRicate (CARD - nt)", 
                  "ABRicate (MEGARES - nt)", "ABRicate (NCBI - nt)", "ABRicate (ResFinder - nt)")

tool_2 <- c("DeepARG (nt)", "RGI (DIAMOND - nt)", "fARGene (nt)", "AMRFinderPlus (aa)", "ResFinder (nt)",
            "ABRicate (ARG-ANNOT - nt)", "ABRicate (CARD - nt)", "ABRicate (MEGARES - nt)", "ABRicate (NCBI - nt)",
            "ABRicate (ResFinder - nt)")


abundance <- abundance %>% mutate( tool = 
                     ifelse(tool == "RGI (BLAST nt)", "RGI (BLAST - nt)",
                     ifelse(tool == "RGI (BLAST aa)", "RGI (BLAST - aa)",
                     ifelse(tool == "RGI (DIAMOND aa)", "RGI (DIAMOND - aa)",
                     ifelse(tool == "RGI (DIAMOND nt)", "RGI (DIAMOND - nt)", tool)))))

# changing habitats and tools to factor
abundance <- abundance %>% mutate(habitat = factor(habitat, levels = EN),
                                  habitat2 = factor(SO[habitat], levels = unique(SO)),
                                  tool = factor(tool, levels = tools_levels))

## factors for new_level, we take the highest abundance per ontology by tool

factor_aro <- abundance %>% ungroup() %>% filter(aggregation %in% "ARO") %>%
  group_by(aggregation, tool, gene ) %>% summarise(total = sum(normed10m)) %>%
  ungroup() %>% arrange(aggregation, tool, desc(total)) %>% 
  group_by(aggregation, tool, gene) %>%  ungroup() %>% select(gene) %>% distinct() %>% pull()
  
factor_parent <- abundance %>% ungroup() %>% filter(aggregation %in% "parent_description") %>%
  group_by(aggregation, tool, gene ) %>% summarise(total = sum(normed10m)) %>%
    ungroup() %>% arrange(aggregation, tool, desc(total)) %>% 
    group_by(aggregation, tool, gene) %>%  ungroup() %>% select(gene) %>% distinct() %>% pull()
  
factor_new_level <- abundance %>% ungroup() %>% filter(aggregation %in% "new_level") %>%
  group_by(aggregation, tool, gene ) %>% summarise(total = sum(normed10m)) %>%
    ungroup() %>% arrange(aggregation, tool, desc(total)) %>% 
    group_by(aggregation, tool, gene) %>%  ungroup() %>% select(gene) %>% distinct() %>% pull()

# this is mainly for the circular / radial plots
factor_new_level2 <- c(factor_new_level[seq(1, length(factor_new_level), by = 3)],
                       factor_new_level[seq(1, length(factor_new_level), by = 3) + 1],
                       factor_new_level[seq(1, length(factor_new_level), by = 3) + 2])
                       #factor_new_level[seq(1, length(factor_new_level), by = 3) + 3])

# factor for habitat

abundance$habitat <- factor(abundance$habitat, levels = EN)
abundance$habitat <- as.character(abundance$habitat)

# environments that we are not interested in
not_env <- c("amplicon", "isolate")

# SUMMARIES
# total numbrer of unigenes captured with abundance 188 441 genes 
# args_abundances  %>% filter(!habitat %in% not_env) %>%
#   select(X) %>% distinct() %>% summarise(n = n())

#double_level <- do.call(rbind, lapply(lst, function(x) x[,c("query","ARO","tool","new_level")])) %>%
#  ungroup() %>% group_by(query) %>% mutate(n = n_distinct(new_level)) %>% ungroup() %>% filter(n>1) %>%  arrange(desc(n), query, tool )

#double_level %>% ungroup() %>% select(query) %>% distinct() %>% pull()
#double_level %>% ungroup() %>% select(ARO) %>% distinct() %>% pull()
#sort(table(double_level %>% ungroup() %>% select(new_level)  %>% pull()))
#double_level %>% group_by(tool, new_level) %>% summarise(n = n()) %>% arrange(desc(n))
#print(double_level, n = 20)
#double_level %>% filter(query %in% (double_level %>% filter(is.na(new_level)) %>% select(query) %>% pull())) %>% print(n =100)
#double_level %>% filter(query %in% (double_level %>% filter(tool %in% "ABRicate (MEGARES - nt)", new_level %in% "Efflux p.") %>% select(query) %>% pull())) %>% group_by(tool, new_level) %>% summarise(n = n())
#double_level %>% filter(query %in% (double_level %>% filter(tool %in% "DeepARG (nt)") %>% select(query) %>% pull())) %>% group_by(tool, new_level) %>% summarise(n = n())
#double_level %>% filter(query %in% (double_level %>% filter(tool %in% "DeepARG (nt)") %>% select(query) %>% pull())) %>% group_by(tool) %>% summarise(n = n())
############## Pan and core resistome 

# unigenes captured with the tools
unigenes <- do.call(rbind, lapply(lst, function(x) x[,c("query","tool", "ARO", "parent", "parent_description", "new_level", "id")])) 
unigenes <- unigenes %>% mutate(tool = factor(tool, levels = tools_levels))

# total number of unique unigenes captured with the tools
unigenes %>% select(query) %>% distinct() %>% summarise(n = n())

# unigenes per tool
unigenes  %>% group_by(tool) %>% summarise(n = n_distinct(query)) %>% ungroup() %>% arrange(n)


plot_count_genes_tool <- unigenes  %>% 
  ggplot(aes( x = tool)) +
  geom_bar(aes(fill = tool), color = "black") +
  scale_fill_manual(values = pal15_rep) +
  theme_minimal() +
  ylab("Number of detected ARGs") +
  xlab("Tool") +
  labs(fill = "") +
  scale_y_continuous(breaks = c(0, 25000, 50000,75000,100000,125000), labels = c("0","25,000", "50,000", "75,000","100,000","125,000")) +
  theme(
    legend.position = "none",
    panel.border = element_rect(fill = "transparent", color = "black", linewidth = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title = element_text(size = general_size + 2, face = "bold"),
    axis.text.x = element_text(angle = 90, size = general_size),
    axis.text.y = element_text(size = general_size))

plot_count_genes_tool

dev.off()
ggsave("~/Documents/plots_project/count_genes_per_tool.svg", plot_count_genes_tool, width = 10, height = 6)
dev.off()  


## UNIGENES missing aro per tool NONE ANYMORE!!!
unigenes %>% 
  group_by(tool) %>% mutate(N = n()) %>% filter(ARO == "") %>% mutate(n = n()) %>% 
  summarise( p = n[1] / N[1], n = n[1])

unigenes %>% 
  group_by(tool) %>% mutate(N = n()) %>% filter(is.na(new_level)) %>% mutate(n = n()) %>% 
  summarise( p = n[1] / N[1])

## NUMBER OF PARENT CLASSES 
unigenes %>%  
  summarise(N = n_distinct(parent_description))

## NUMBER OF new levels
unigenes %>% 
  summarise(N = n_distinct(new_level))


### 

get_accumulated_count <- function(x) {
  y <- x %>% group_by(tool, new_level) %>% 
    summarise(n = n()) %>% 
    arrange(n) %>% 
    ungroup() %>% 
    mutate(prop = cumsum(n) / sum(n))
  return(y)
}

# NUMBER OF GENES PER CLASS PER TOOL  ORDERED BY PROPORTION ON THE TOOL
class_per_tool <- do.call(rbind, lapply(lst, function(x) get_accumulated_count(x)))

# classes contributing with at least 50%  of genes in each tool
data.frame(class_per_tool %>% group_by(tool) %>% 
  arrange(desc(prop), .by_group = TRUE) %>%  # or arrange by time/order column if you have one
  mutate(next_x = lead(prop)) %>%
  filter(prop > 0.5 | (lag(prop) > 0.5 & prop <= 0.5)) %>%
  select(-next_x))

# NUMBER OF GENES PER CLASS PER TOOL, 
# N total genes in all tools and classes, 
# M total genes per class in all tools
# P proportion of class in total (all tools)
# ntool total genes per tool
# n total genes per class in specific tool
# p proportion of class in specific tool

proportion_new_level_tool <- unigenes %>% ungroup() %>%  
  mutate(N = n_distinct(query)) %>% 
  group_by(new_level) %>% mutate(M = n_distinct(query)) %>% ungroup() %>%
  group_by(tool) %>% mutate(ntool = n_distinct(query)) %>% ungroup() %>%
  mutate(P = M / N) %>%
  group_by(new_level, tool) %>% summarise(N = N[1], M = M[1], P = P[1], ntool = ntool[1], n = n_distinct(query)) %>% 
  arrange(desc(n)) %>% mutate(p = n/ntool)


# sum of propotions for GPA AND EFFLUX
proportion_new_level_tool %>% filter(new_level %in% c("GPA","Efflux p.")) %>% group_by(tool) %>% summarise(N = N[1], M = sum(M), ntool = ntool[1], n = sum(n), p = sum(p))

# sum of propotions for BETA-LACTAM
proportion_new_level_tool %>% filter(new_level %in% c("Class A", "Class B", "Class D", "Class C")) %>% group_by(tool) %>% summarise(N = N[1], M = sum(M), ntool = ntool[1], n = sum(n), p = sum(p))


# proportion of class among all unigenes detected
top10class <- unigenes %>% 
  mutate(N = n_distinct(query)) %>% 
  group_by(new_level) %>% summarise(N = N[1], n = n_distinct(query)) %>% 
  arrange(desc(n)) %>% mutate(p = n/N) %>% 
  select(new_level) %>% slice_head(n=10) %>% pull()

top10class <- factor(c(top10class, "Other"), 
                     levels = c(top10class, "Other"))

top10class_plot <- bind_rows(proportion_new_level_tool %>% filter(new_level %in% top10class, tool %in% tool_2) %>% 
            select(new_level, tool, p) %>% ungroup() %>% group_by(tool) %>% mutate(P = sum(p)),
  
          proportion_new_level_tool %>% filter(new_level %in% top10class, tool %in% tool_2) %>% 
            select(new_level, tool, p) %>% ungroup() %>% group_by(tool) %>% mutate(P = sum(p)) %>% 
            summarise(p = 1 - max(P)) %>% mutate(new_level = "Other")) %>% 
  mutate(tool = factor(tool, levels = tool_2),
         new_level = factor(new_level, levels = top10class)) %>% 
  ggplot(aes( x = tool, y = p)) +
  geom_col(aes(fill = new_level), color = "black") +
  scale_fill_manual(values = pal_11_d) +
  theme_minimal() +
  #scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
  #              labels = scales::trans_format("log10", scales::math_format(10^.x)))  + 
  ylab("Number of detected ARGs") +
  xlab("Tool") +
  labs(fill = "") +
  theme(
    #legend.position = "none",
    panel.border = element_rect(fill = "transparent", color = "black", linewidth = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title = element_text(size = general_size + 2, face = "bold"),
    axis.text.x = element_text(angle = 90, size = general_size),
    axis.text.y = element_text(size = general_size))

top10class_plot

dev.off()
ggsave("~/Documents/plots_project/top10classes.svg", top10class_plot, width = 10, height = 6)
dev.off()  


##############################################################
##############################################################

abundance_plot <- abundance  %>% 
  group_by(habitat2, tool, sample, aggregation) %>% summarise(total = sum(normed10m)+1) %>%
  #filter(!habitat %in% not_env, tool %in% tool_selected, !sample %in% extreme_samples) %>% 
  filter(tool %in% tool_2, aggregation %in% "new_level") %>% 
  filter(!habitat2 %in% "Other") %>% 
  ggplot(aes( x = habitat2)) +
  geom_boxplot(aes(y = total, fill = tool), outlier.shape = NA) +
  scale_fill_manual(values = pal_10_d) +
  theme_minimal() +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)),
                limits = c(1, 20000)) + 
  ylab("Normalized abundance") +
  xlab("Source") +
  labs(fill = "") +
  theme(
    legend.position = "bottom",
    panel.border = element_rect(fill = "transparent", color = "black", linewidth = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title = element_text(size = general_size + 2, face = "bold"),
    axis.text.x = element_text(angle = 90, size = general_size),
    axis.text.y = element_text(size = general_size))

abundance_plot

dev.off()
ggsave("~/Documents/plots_project/abundance.svg", abundance_plot, width = 10, height = 6)
dev.off()  

##############################################################
##############################################################

diversity_plot <- abundance  %>% 
  group_by(habitat2, tool, sample, aggregation) %>% summarise(total = sum(unigenes) + 1) %>%
  #filter(!habitat %in% not_env, tool %in% tool_selected, !sample %in% extreme_samples) %>% 
  filter(tool %in% tool_2, aggregation %in% "new_level") %>% 
  filter(!habitat2 %in% "Other") %>% 
  ggplot(aes( x = habitat2)) +
  geom_boxplot(aes(y = total, fill = tool), outlier.shape = NA) +
  scale_fill_manual(values = pal_10_d) +
  theme_minimal() +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)),
                limits = c(1, 6000)) + 
  ylab("Diversity") +
  xlab("Source") +
  labs(fill = "") +
  theme(
    legend.position = "bottom",
    panel.border = element_rect(fill = "transparent", color = "black", linewidth = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title = element_text(size = general_size + 2, face = "bold"),
    axis.text.x = element_text(angle = 90, size = general_size),
    axis.text.y = element_text(size = general_size))

diversity_plot

dev.off()
ggsave("~/Documents/plots_project/diversity.svg", diversity_plot, width = 10, height = 6)
dev.off()  

##############################################################
##############################################################

abu_tool_habitat <- abundance  %>% 
  group_by(habitat, tool, sample, aggregation, gene) %>% 
  summarise(abundance = sum(normed10m) + 1,
            diversity = sum(unigenes) + 1) %>%
  filter(!habitat %in% not_env, tool %in% tool_2) 

# tool_2[c(1:3,5)] = "DeepARG (nt)" "RGI (DIAMOND - nt)" "fARGene (nt)" "ResFinder (nt)" 
human.genes <- abu_tool_habitat  %>% 
  filter(habitat %in% c("human gut"), aggregation %in% "new_level", tool %in% tool_2[c(1:3,5)]) %>%
  ungroup() %>% group_by( tool, gene) %>% summarise(n = median(abundance)) %>%
  ungroup() %>% arrange( tool, desc(n)) %>% group_by(tool)  %>% slice_head(n = 5) %>% 
  ungroup() %>% 
  select(gene) %>% 
  distinct() %>% 
  pull()

abundance_plot_habitat <- 
  abu_tool_habitat  %>% 
  filter(habitat %in% c("human gut"), aggregation %in% "new_level", tool %in% tool_2[c(1:3,5)], gene %in% human.genes) %>%
  mutate(gene = factor(gene, levels = factor_new_level)) %>% 
  ggplot(aes( x = gene)) +
  geom_boxplot(aes(y = abundance, fill = tool), outlier.shape = NA, position = position_dodge(preserve = "single")) +
  scale_fill_manual(values = pal_4_d) +
  theme_minimal() +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)),
                limits = c(1, 20000)) + 
  ylab("Normalized abundance") +
  xlab("Class") +
  labs(fill = "") +
  facet_grid(habitat ~ gene, scales = "free_x") +
  theme(
    legend.position = "bottom",
    panel.border = element_rect(fill = "transparent", color = "black", linewidth = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    strip.text = element_blank(),
    axis.title = element_text(size = general_size + 2, face = "bold"),
    axis.text.x = element_text(angle = 90, size = general_size),
    axis.text.y = element_text(size = general_size))

abundance_plot_habitat


diversity_plot_habitat <- 
  abu_tool_habitat  %>% 
  filter(habitat %in% c("human gut"), aggregation %in% "new_level", tool %in% tool_2[c(1:3,5)], gene %in% human.genes) %>%
  mutate(gene = factor(gene, levels = factor_new_level)) %>% 
  ggplot(aes( x = gene)) +
  geom_boxplot(aes(y = diversity, fill = tool), outlier.shape = NA, position = position_dodge(preserve = "single")) +
  scale_fill_manual(values = pal_4_d) +
  theme_minimal() +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)),
                limits = c(1, 1300)) + 
  ylab("Diversity") +
  xlab("Class") +
  labs(fill = "") +
  facet_grid(habitat ~ gene, scales = "free_x") +
  theme(
    legend.position = "bottom",
    panel.border = element_rect(fill = "transparent", color = "black", linewidth = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    strip.text = element_blank(),
    axis.title = element_text(size = general_size + 2, face = "bold"),
    axis.text.x = element_text(angle = 90, size = general_size),
    axis.text.y = element_text(size = general_size))

diversity_plot_habitat

dev.off()
ggsave("~/Documents/plots_project/abundance_class_tool.svg", abundance_plot_habitat, width = 10, height = 6)
dev.off()  

dev.off()
ggsave("~/Documents/plots_project/diversity_class_tool.svg", diversity_plot_habitat, width = 10, height = 6)
dev.off()  


################################
##############################################################
##############################################################

### conformity (overlap)

## genes per tool and ID level
rownames(unigenes) <- NULL


tools_per_unigene <- unigenes %>% ungroup()  %>% 
  arrange(query) %>% 
  group_by(query) %>% 
  mutate(n_tools = n_distinct(tool)) %>% 
  mutate(single = (n_tools ==1)) 


# overlap all genes between tools 
# overlap all genes between tools 
# overlap all genes between tools 
# overlap all genes between tools 

# create lists
sets <- tools_per_unigene %>%  
  group_by(tool) %>%
  summarise(query = list(query), .groups = "drop") # put every query in a list

# Create all pairwise combinations
pairwise <- expand_grid(tool_ref = sets$tool, tool_comp = sets$tool) 

# Compute Jaccard, recall, fnr, for each pair of tools using the set lists 
JI_all <- pairwise %>%
  left_join(sets, by = c("tool_ref" = "tool")) %>%
  rename(values1 = query) %>%
  left_join(sets, by = c("tool_comp" = "tool")) %>%
  rename(values2 = query) %>%
  mutate(jaccard = map2_dbl(values1, values2, ~ length(intersect(.x, .y)) / length(union(.x, .y)))) %>%
  mutate(recall = map2_dbl(values1, values2, ~ length(intersect(.x, .y)) / length( .y))) %>%
  mutate(fnr = map2_dbl(values1, values2, ~ length(setdiff(.x, .y)) / length( .x))) %>%
  select(tool_ref, tool_comp, jaccard, recall, fnr)  
data.frame(JI_all)

# overlap all genes between tools and classes
# create lists

sets <- tools_per_unigene %>%
  group_by(new_level, tool) %>% 
  summarise(query = list(query), .groups = "drop")

# Pairwise combinations within each d
pairwise <- sets %>%
  group_by(new_level) %>%
  summarise(pairs = list(expand_grid(tool_ref = tool, tool_comp = tool)), .groups = "drop") %>%
  unnest(pairs)

# Add the value lists to each pair
JI_class <- pairwise %>%
  left_join(sets, by = c("new_level", "tool_ref" = "tool")) %>%
  rename(values1 = query) %>%
  left_join(sets, by = c("new_level", "tool_comp" = "tool")) %>%
  rename(values2 = query) %>%
  mutate(jaccard = map2_dbl(values1, values2, ~ length(intersect(.x, .y)) / length(union(.x, .y)))) %>%
  mutate(recall = map2_dbl(values1, values2, ~ length(intersect(.x, .y)) / length( .y))) %>%
  mutate(fnr = map2_dbl(values1, values2, ~ length(setdiff(.x, .y)) / length( .x))) %>%
  select(new_level, tool_ref, tool_comp, jaccard, recall, fnr)


JI_class %>% filter(tool_ref %in% "RGI (DIAMOND - nt)")
JI_class %>% filter(tool_ref %in% "fARGene (nt)")
JI_class <- JI_class %>% filter(tool_ref != tool_comp) %>% ungroup()


# double new level
double_level <- do.call(rbind, lapply(lst, function(x) x[,c("query","tool","new_level")])) %>%
  ungroup() %>% 
  group_by(query) %>% 
  summarise(n = n_distinct(new_level)) %>% 
  filter(n>1) %>% select(query) %>% 
  pull()

# double new level for unigenes present in fargene 
unigenes %>% 
  filter(query %in% 
           (unigenes %>% filter(query %in% double_level, tool %in% c("fARGene (nt)","fARGene (aa)")) %>% 
              select(query) %>% pull())) %>% arrange(query, id, tool) %>% select(query, tool, new_level, id)

unigenes %>% filter(tool %in% c("fARGene (nt)","fARGene (aa)")) %>% group_by(query) %>% summarise(n = n_distinct(new_level)) %>% filter(n > 1)

# double_new_level_count = pairwise count for each new level
sets <- unigenes %>% 
  filter(query %in% double_level) %>% 
  ungroup() %>% group_by(query, new_level) %>% distinct() %>%
  group_by(new_level) %>%  
  summarise(query = list(query), .groups = "drop")

pairwise <- expand_grid(nl1 = sets$new_level, nl2 = sets$new_level) 

double_new_level_count <- pairwise %>%
  left_join(sets, by = c("nl1" = "new_level")) %>%
  rename(values1 = query) %>%
  left_join(sets, by = c("nl2" = "new_level")) %>%
  rename(values2 = query) %>%
  filter(nl1 != nl2) %>% 
  rowwise() %>%
  mutate( pair1 = min(c(nl1, nl2)), pair2 = max(c(nl1, nl2))) %>%
  ungroup() %>% 
  distinct(pair1, pair2, .keep_all = TRUE) %>% 
  mutate(n_genes = map2_dbl(values1, values2, ~ length(intersect(.x, .y)))) %>%
  filter(n_genes > 0 ) %>% 
  select(nl1, nl2, n_genes) %>% arrange(desc(n_genes))

data.frame(double_new_level_count)


# double_tool_count = pairwise count for tool

sets <- unigenes %>% 
  filter(query %in% double_level) %>% 
  group_by(tool) %>%  
  summarise(query = list(query),
            new_level = list(new_level), .groups = "drop")

pairwise <- expand_grid(nl1 = sets$tool, nl2 = sets$tool) %>% filter(nl1 != nl2)

library(purrr)


pull_unigenes <- function(x,y){
  g1 <- unigenes %>% ungroup() %>% 
    filter(query %in% double_level) %>% 
    group_by(tool) %>% filter(tool %in% x) %>% select(query, new_level)
  g2 <- unigenes %>% ungroup() %>% 
    filter(query %in% double_level) %>% 
    group_by(tool) %>% filter(tool %in% y) %>% select(query, new_level)
  g <-  bind_rows(g1, g2)
  g <- g %>%group_by(query) %>% summarise(n = n_distinct(new_level)) %>% filter(n>1) %>% ungroup()
  g <- g %>% select(query) %>% pull()
  return(g)
}


double_tool_count <- pairwise %>%
  mutate(q = map2(nl1, nl2, ~  pull_unigenes(as.character(.x), as.character(.y)))) %>% 
  filter(nl1 != nl2) %>% 
  rowwise() %>%
  mutate(n_genes = length(q)) %>%
  mutate( pair1 = min(c(as.character(nl1), as.character(nl2))), pair2 = max(c(as.character(nl1), as.character(nl2)))) %>% 
  ungroup() %>% 
  distinct(pair1, pair2, .keep_all = TRUE) %>% 
  filter(n_genes > 0 ) %>% 
  select(nl1, nl2, n_genes) %>% arrange(nl1, desc(n_genes))

data.frame(double_tool_count %>% filter(nl1 %in% tool_2, nl2 %in% tool_2))



# now, we have a problem, the recall will be biased if two tools report one unigene with different class

# overlap all genes between tools and classes
# 
#
# create lists

# sets0 <- tools_per_unigene %>%
#   group_by(tool) %>%
#   summarise(query = list(query), .groups = "drop")
# 
# sets1 <- tools_per_unigene %>%
#   group_by(new_level, tool) %>%
#   summarise(query = list(query), .groups = "drop")
# 
# # Pairwise combinations within each d
# pairwise <- sets %>%
#   group_by(new_level) %>%
#   summarise(pairs = list(expand_grid(tool_ref = tool, tool_comp = tool)), .groups = "drop") %>%
#   unnest(pairs)
# 
# # Add the value lists to each pair
# JI_class_vs_all <- pairwise %>%
#   left_join(sets1, by = c("new_level", "tool_ref" = "tool")) %>%
#   rename(values1 = query) %>%
#   left_join(sets0, by = c( "tool_comp" = "tool")) %>%
#   rename(values2 = query) %>%
#   mutate(jaccard = map2_dbl(values1, values2, ~ length(intersect(.x, .y)) / length(union(.x, .y)))) %>%
#   mutate(recall = map2_dbl(values1, values2, ~ length(intersect(.x, .y)) / length( .y))) %>%
#   mutate(fnr = map2_dbl(values1, values2, ~ length(setdiff(.x, .y)) / length( .x))) %>%
#   select(new_level, tool_ref, tool_comp, jaccard, recall, fnr)
# 
# 
# JI_class %>% filter(tool_ref %in% "RGI_DIAMOND_nt")
# JI_class %>% filter(tool_ref %in% "fARGene_nt")
# JI_class_vs_all %>% filter(tool_ref %in% "RGI_DIAMOND_nt")
# JI_class_vs_all %>% filter(tool_ref %in% "fARGene_nt")
# rm (sets0, sets1, sets, pairwise)
# rm (sets0, sets1, sets, pairwise)



JI_all %>% 
  ggplot(aes(x = tool_comp, y = tool_ref, fill = jaccard)) +
  geom_tile() +
  #scale_fill_gradient2(low = "white",  high = "red") +
  scale_fill_viridis_c() + 
  theme_minimal() +
  labs(fill = "Jaccaard index") +
  xlab("") +
  theme(axis.title.y = element_blank(),
        legend.position = "bottom",
        panel.border = element_rect(fill = "transparent", color = "black", linewidth = 1),
        axis.text.x = element_text(angle = 90, size = general_size),
        axis.text.y = element_text(size = general_size),
        strip.text = element_text(size = general_size, face = "bold"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title = element_text(size = general_size + 2, face = "bold"))


JI_all %>% 
  ggplot(aes(x = tool_comp, y = tool_ref, fill = fnr)) +
  geom_tile() +
  #scale_fill_gradient2(low = "white",  high = "red") +
  scale_fill_viridis_c() + 
  theme_minimal() +
  labs(fill = "Genes not found in opposite tool") +
  xlab("") +
  theme(axis.title.y = element_blank(),
        legend.position = "bottom",
        panel.border = element_rect(fill = "transparent", color = "black", linewidth = 1),
        axis.text.x = element_text(angle = 90, size = general_size),
        axis.text.y = element_text(size = general_size),
        strip.text = element_text(size = general_size, face = "bold"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title = element_text(size = general_size + 2, face = "bold"))


JI_all %>% 
  ggplot(aes(x = tool_comp, y = tool_ref, fill = recall)) +
  geom_tile() +
  scale_fill_viridis_c() + 
  theme_minimal() +
  labs(fill = "Genes found in opposite tool") +
  xlab("") +
  theme(axis.title.y = element_blank(),
        legend.position = "bottom",
        panel.border = element_rect(fill = "transparent", color = "black", linewidth = 1),
        axis.text.x = element_text(angle = 90, size = general_size),
        axis.text.y = element_text(size = general_size),
        strip.text = element_text(size = general_size, face = "bold"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title = element_text(size = general_size + 2, face = "bold"))


# summary 1 - calculation of recall and fnr per tool and gene class


plot_recall <- JI_class_filter %>% filter(tool_ref %in% c("RGI (DIAMOND - nt)", "DeepARG (nt)", "fARGene (nt)"),
                                          tool_comp %in% tool_2,
                                          new_level %in% c(human.genes, "Class C", "Class D", "MPH", "APH", "QNR")) %>% 
  ggplot(aes(x = new_level, y = recall)) +
  geom_boxplot(aes(fill = tool_ref)) +
  facet_grid( tool_ref ~ new_level , scales = "free_x") +
  scale_color_manual(values = pal_10_d) +
  scale_fill_manual(values = pal_10_d) +
  theme_minimal() +
  ylab("Robustness") +
  xlab("Class") +
  labs(fill = "") +
  theme(
    legend.position = "bottom",
    panel.border = element_rect(fill = "transparent", color = "black", linewidth = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    strip.text.y = element_text(size = general_size + 2, face = "bold"),
    strip.text.x = element_blank(),
    axis.title = element_text(size = general_size + 2, face = "bold"),
    axis.text.x = element_text(angle = 90, size = general_size),
    axis.text.y = element_text(size = general_size))

plot_recall

dev.off()
ggsave("~/Documents/plots_project/recall_box_1.svg", plot_recall, width = 18, height = 8)
dev.off()  




plot_fdr <- JI_class_filter %>% filter(tool_ref %in% c("RGI (DIAMOND - nt)", "DeepARG (nt)", "fARGene (nt)"),
                                       tool_comp %in% tool_2,
                                       new_level %in% c(human.genes, "Class C", "Class D", "MPH", "APH", "QNR")) %>% 
  ggplot(aes(x = new_level, y = fnr)) +
  geom_boxplot(aes(fill = tool_ref)) +
  facet_grid( tool_ref ~ new_level , scales = "free_x") +
  scale_color_manual(values = pal_10_d) +
  scale_fill_manual(values = pal_10_d) +
  theme_minimal() +
  ylab("Unique") +
  xlab("Class") +
  labs(fill = "") +
  theme(
    legend.position = "bottom",
    panel.border = element_rect(fill = "transparent", color = "black", linewidth = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    strip.text.y = element_text(size = general_size + 2, face = "bold"),
    strip.text.x = element_blank(),
    axis.title = element_text(size = general_size + 2, face = "bold"),
    axis.text.x = element_text(angle = 90, size = general_size),
    axis.text.y = element_text(size = general_size))

plot_fdr

dev.off()
ggsave("~/Documents/plots_project/fnr_box_1.svg", plot_fdr, width = 18, height = 8)
dev.off()  



data.frame(JI_class_filter %>% filter(tool_ref %in% "AMRFinderPlus (aa)", tool_comp %in% tool_2) %>% 
             group_by(tool_ref, new_level) %>% summarise(M_R = median(recall)))


data.frame(hm_recall_fnr_1 %>% filter(tool_ref %in% "AMRFinderPlus (aa)" & !tool_comp %in% "fARGene_nt" ) %>% 
             group_by(tool_ref, new_level) %>% summarise(M_R = median(recall)))


# summary 2 -  medians of recall and fnr per tool and gene class, accross all other tools

hm_recall_fnr_2 <- hm_recall_fnr_1 %>% 
  group_by(tool, new_level) %>%
  summarise(mn_recall = mean(recall, na.rm = T), med_recall = median(recall, na.rm = T),
            mn_fnr = mean(fnr, na.rm = T), med_fnr = median(fnr, na.rm = T),
            iqr_recall = IQR(recall, na.rm = T), iqr_fnr = IQR(fnr, na.rm = T),
            sum_row = sum(row_sum - value, na.rm = T))

# summary 2 -  medians of the medians of recall and fnr per tool accross gene class after accross tool

 hm_recall_fnr_2 %>% filter(!tool %in% c("fARGene (nt)", "AMRFinderPlus (aa)")) %>% 
  ungroup() %>%
  group_by(tool) %>% 
  summarise(MN_recall = mean(mn_recall, na.rm = T), MED_recall = median(med_recall, na.rm = T),
            MN_fnr = mean(mn_fnr, na.rm = T), MED_fnr = median(med_fnr, na.rm = T),
            iqr_med_recall = IQR(med_recall, na.rm = T), iqr_med_fnr = IQR(med_fnr, na.rm = T),
            SUM_row = sum(sum_row, na.rm = T))
 
 hm_recall_fnr_2 %>% filter(tool %in% c("fARGene (nt)", "AMRFinderPlus (aa)"), !is.na(med_recall)) 
 hm_recall_fnr_2 %>% filter(tool %in% c("fARGene (nt)"), !is.na(med_recall)) 
 hm_recall_fnr_2 %>% filter(tool %in% c("AMRFinderPlus (aa)"), !is.na(med_recall)) 

hm_recall_fnr_3 <- hm_recall_fnr_2 %>%
  ungroup() %>%
  group_by(tool) %>% 
  summarise(MN_recall = mean(mn_recall, na.rm = T), MED_recall = median(med_recall, na.rm = T),
            MN_fnr = mean(mn_fnr, na.rm = T), MED_fnr = median(med_fnr, na.rm = T),
            iqr_med_recall = IQR(med_recall, na.rm = T), iqr_med_fnr = IQR(med_fnr, na.rm = T),
            SUM_row = sum(sum_row, na.rm = T))

# unique genes per tool
tools_per_unigene %>% group_by(query) %>% filter(tool == "AMRFinderPlus_aa") %>% 
  ungroup() %>% group_by(n_tools) %>% summarise(n = n()) %>% mutate( p = n/sum(n)) %>% filter(n_tools == 1)

tools_per_unigene %>% group_by(query) %>% filter(tool == "ABRicate_ResFinder_nt") %>% 
  ungroup() %>% group_by(n_tools) %>% summarise(n = n()) %>% mutate( p = n/sum(n)) %>% filter(n_tools == 1)

tools_per_unigene %>% group_by(query) %>% filter(tool == "ABRicate_CARD_nt") %>% 
  ungroup() %>% group_by(n_tools) %>% summarise(n = n()) %>% mutate( p = n/sum(n)) %>% filter(n_tools == 1)

tools_per_unigene %>% group_by(query) %>% filter(tool == "ABRicate_ARGANNOT_nt") %>% 
  ungroup() %>% group_by(n_tools) %>% summarise(n = n()) %>% mutate( p = n/sum(n)) %>% filter(n_tools == 1)

tools_per_unigene %>% group_by(query) %>% filter(tool == "ABRicate_NCBI_nt") %>% 
  ungroup() %>% group_by(n_tools) %>% summarise(n = n()) %>% mutate( p = n/sum(n)) %>% filter(n_tools == 1)

tools_per_unigene %>% group_by(query) %>% filter(tool == "ABRicate_MEGARES_nt") %>% 
  ungroup() %>% group_by(n_tools) %>% summarise(n = n()) %>% mutate( p = n/sum(n)) %>% filter(n_tools == 1)

tools_per_unigene %>% group_by(query) %>% filter(tool == "RGI_DIAMOND_nt") %>% 
  ungroup() %>% group_by(n_tools) %>% summarise(n = n()) %>% mutate( p = n/sum(n)) %>% filter(n_tools == 1)

tools_per_unigene %>% group_by(query) %>% filter(tool == "DeepARG_nt") %>% 
  ungroup() %>% group_by(n_tools, new_level) %>% summarise(n = n())  %>% filter(n_tools == 1) %>% arrange(desc(n))



# ggplot(hm_recall_fnr_2, aes(x = new_level, y = tool, size = 1 - iqr_recall, color = exp(med_recall))) +
#   geom_point() +
#   scale_size_continuous(range = c(0.05, 10)) +
#   scale_color_gradient2(low = "white",  high = "red", midpoint = max(exp(1))/2) +
#   theme_minimal() +
#   theme(axis.title.y = element_blank(),
#         axis.ticks.y = element_blank())

alt_name_tools_rev <- names(alt_name_tools)
names(alt_name_tools_rev) <- alt_name_tools

 

#hm_recall_fnr_2 <- hm_recall_fnr_2 %>% mutate(tool = as.vector(alt_name_tools_rev[as.character(tool)]))
#hm_recall_fnr_2 <- hm_recall_fnr_2 %>% mutate(tool = factor(tool, levels = tools_levels))

general_size <- 10

plot_recall <- ggplot(hm_recall_fnr_2, aes(x = tool, y = med_recall, fill = tool)) +
  geom_boxplot() +
  scale_fill_manual(values = pal_10_d) +
  theme_minimal() +
  ylab("Recall") +
  xlab("Tool") +
  theme(panel.border = element_rect(fill = "transparent", color = "black", linewidth = 2),
        legend.position = "none",
        axis.text.x = element_text(angle = 90, size = general_size),
        axis.text.y = element_text(size = general_size),
        axis.title = element_text(size = general_size + 2, face = "bold"))

plot_recall

plot_fnr <- ggplot(hm_recall_fnr_2, aes(x = tool, y = med_fnr, fill = tool)) +
  geom_boxplot() +
  scale_fill_manual(values = pal_10_d) +
  theme_minimal() +
  ylab("False Negative Rate") +
  xlab("Tool") +
  #scale_x_continuous(breaks = c(0, 0.5, 1), labels = c(0, 0.5, 1)) +
  theme(panel.border = element_rect(fill = "transparent", color = "black", linewidth = 2),
        legend.position = "none",
        axis.text.x = element_text(angle = 90, size = general_size),
        axis.text.y = element_text(size = general_size),
        axis.title = element_text(size = general_size + 2, face = "bold"))

plot_fnr


#df$facet_var <- gsub(" ", "\n", df$facet_var)

hm_recall_fnr_2 <- hm_recall_fnr_2 %>% 
  filter(!is.na(new_level)) %>% 
  mutate(strip = gsub("\\(", "\n\\(", tool)) %>%
  mutate(strip = factor(strip, levels = gsub("\\(", "\n\\(", levels(tool))))

recall_heatmap <- ggplot(hm_recall_fnr_2, aes(x = med_recall, y = new_level, fill = 1-iqr_recall)) +
  geom_col() +
  #scale_fill_gradient2(low = "white",  high = "red") +
  scale_fill_viridis_c() + 
  facet_grid(. ~ strip) +
  theme_minimal() +
  labs(fill = "1 - IQR") +
  xlab("Recall") +
  scale_x_continuous(breaks = c(0, 0.5, 1), labels = c(0, 0.5, 1)) +
  theme(axis.title.y = element_blank(),
        legend.position = "bottom",
        panel.border = element_rect(fill = "transparent", color = "black", linewidth = 1),
        axis.text.x = element_text( size = general_size),
        axis.text.y = element_text(size = general_size),
        strip.text = element_text(size = general_size, face = "bold"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title = element_text(size = general_size + 2, face = "bold"))

recall_heatmap

fnr_heatmap <- ggplot(hm_recall_fnr_2 %>% filter(!is.na(new_level)), aes(x = med_fnr, y = new_level, fill = 1-iqr_fnr)) +
  geom_col() +
  #scale_fill_gradient2(low = "white",  high = "red") +
  scale_fill_viridis_c() + 
  facet_grid(. ~ strip) +
  theme_minimal() +
  labs(fill = "1 - IQR") +
  xlab("False Negative Rate") +
  scale_x_continuous(breaks = c(0, 0.5, 1), labels = c(0, 0.5, 1)) +
  theme(axis.title.y = element_blank(),
        legend.position = "bottom",
        panel.border = element_rect(fill = "transparent", color = "black", linewidth = 1),
        axis.text.x = element_text( size = general_size),
        axis.text.y = element_text(size = general_size),
        strip.text = element_text(size = general_size, face = "bold"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title = element_text(size = general_size + 2, face = "bold"))

fnr_heatmap


dev.off()
ggsave("~/Documents/plots_project/recall.svg", plot_recall, width = 10, height = 6)
dev.off()

dev.off()
ggsave("~/Documents/plots_project/fnr.svg", plot_fnr, width = 10, height = 6)
dev.off()

dev.off()
ggsave("~/Documents/plots_project/recall_heatmap.svg", recall_heatmap  , width = 14, height = 12)
dev.off()

dev.off()
ggsave("~/Documents/plots_project/fnr_heatmap.svg", fnr_heatmap, width = 14, height = 12)
dev.off()


JI %>% 
  mutate(x1 = factor(as.vector(alt_name_tools_rev[as.character(tool)]), levels = tools_levels)) %>% 
  mutate(x2 = factor(as.vector(alt_name_tools_rev[as.character(variable)]), levels = tools_levels)) %>%
  ggplot(aes(x = x1, y = x2, fill = p3)) +
  geom_tile() +
  #scale_fill_gradient2(low = "white",  high = "red") +
  scale_fill_viridis_c() + 
  theme_minimal() +
  labs(fill = "Jaccaard index") +
  xlab("") +
  theme(axis.title.y = element_blank(),
        legend.position = "bottom",
        panel.border = element_rect(fill = "transparent", color = "black", linewidth = 1),
        axis.text.x = element_text(angle = 90, size = general_size),
        axis.text.y = element_text(size = general_size),
        strip.text = element_text(size = general_size, face = "bold"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title = element_text(size = general_size + 2, face = "bold"))




idplot <- unigenes_per_tool %>% filter(tool %in% tool_2[1:4]) %>% 
  group_by(query) %>% 
  mutate(n = n_distinct(tool))  %>% 
  mutate(n = n > 1) %>% 
ggplot( aes(id, colour = tool, linetype = n)) +
  geom_freqpoly(binwidth = 1, linewidth = 2) +
  facet_wrap(. ~ tool, scales = "free_y") +
  scale_color_manual(values = rep(pal_10_d[9],10)) +
  theme_minimal() +
  labs(fill = "Jaccaard index") +
  xlab("Identity") +
  ylab("Unigenes") +
  theme(axis.title.y = element_blank(),
        legend.position = "none",
        panel.border = element_rect(fill = "transparent", color = "black", linewidth = 1),
        axis.text.x = element_text(angle = 90, size = general_size),
        axis.text.y = element_text(size = general_size),
        strip.text = element_text(size = general_size, face = "bold"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title = element_text(size = general_size + 2, face = "bold"))


dev.off()
ggsave("~/Documents/plots_project/ids.svg", idplot, width = 10, height = 8)


#######





table(lst$amrfinder.norm.prot$ARG.class, lst$amrfinder.norm.prot$Method)

unigenes_per_tool %>% filter(tool %in% tool_2 & id < 80) %>% 
  group_by(query) %>% 
  mutate(n = n_distinct(tool))  %>% 
  mutate(n = n > 1) %>% 
  ggplot( aes(id, colour = tool, linetype = n)) +
  geom_freqpoly(binwidth = 1) +
  facet_wrap(. ~ tool, scales = "free_y")





t70 <- unigenes_per_tool %>% 
  group_by(tool) %>% mutate(total_tool = n()) %>% ungroup() %>%
  group_by(tool, new_level) %>% mutate(total_class = n()) %>% ungroup() %>%
  filter(tool %in% c("RGI (DIAMOND - nt)", "DeepARG (nt)") & id < 70) %>% 
  group_by(query) %>% 
  mutate(n = n_distinct(tool))  %>% 
  mutate(n = n > 1) %>% 
  ungroup() %>% group_by(tool, new_level, n) %>% 
  summarise(total_tool = total_tool[1], total_class = total_class[1], N_tool_class_bool = n()) %>% 
  mutate(p_bool = N_tool_class_bool / sum(N_tool_class_bool)) %>%
  ungroup() %>% group_by(tool) %>% mutate(p_tool = N_tool_class_bool / sum(N_tool_class_bool)) %>%
  mutate(P = N_tool_class_bool / total_class) %>%
  arrange(tool, desc(N_tool_class_bool), desc(P), desc(p_tool), desc(p_bool)) 

t60 <- unigenes_per_tool %>% 
  group_by(tool) %>% mutate(total_tool = n()) %>% ungroup() %>%
  group_by(tool, new_level) %>% mutate(total_class = n()) %>% ungroup() %>%
  filter(tool %in% c("RGI (DIAMOND - nt)", "DeepARG (nt)") & id < 60) %>% 
  group_by(query) %>% 
  mutate(n = n_distinct(tool))  %>% 
  mutate(n = n > 1) %>% 
  ungroup() %>% group_by(tool, new_level, n) %>% 
  summarise(total_tool = total_tool[1], total_class = total_class[1], N_tool_class_bool = n()) %>% 
  mutate(p_bool = N_tool_class_bool / sum(N_tool_class_bool)) %>%
  ungroup() %>% group_by(tool) %>% mutate(p_tool = N_tool_class_bool / sum(N_tool_class_bool)) %>%
  mutate(P = N_tool_class_bool / total_class) %>%
  arrange(tool, desc(N_tool_class_bool), desc(P), desc(p_tool), desc(p_bool))

t50 <- unigenes_per_tool %>% 
  group_by(tool) %>% mutate(total_tool = n()) %>% ungroup() %>%
  group_by(tool, new_level) %>% mutate(total_class = n()) %>% ungroup() %>%
  filter(tool %in% c("RGI (DIAMOND - nt)", "DeepARG (nt)") & id < 50) %>% 
  group_by(query) %>% 
  mutate(n = n_distinct(tool))  %>% 
  mutate(n = n > 1) %>% 
  ungroup() %>% group_by(tool, new_level, n) %>% 
  summarise(total_tool = total_tool[1], total_class = total_class[1], N_tool_class_bool = n()) %>% 
  mutate(p_bool = N_tool_class_bool / sum(N_tool_class_bool)) %>%
  ungroup() %>% group_by(tool) %>% mutate(p_tool = N_tool_class_bool / sum(N_tool_class_bool)) %>%
  mutate(P = N_tool_class_bool / total_class) %>%
  arrange(tool, desc(N_tool_class_bool), desc(P), desc(p_tool), desc(p_bool))

data.frame(t60)
t70 %>% filter(!n) %>% ungroup() %>% group_by(tool) %>% summarise(n = n(), min_ptool = min(p_tool), min_pbool = min(p_bool), md_bool = median(p_bool), mn_bool = mean(p_bool))
t60 %>% filter(!n) %>% ungroup() %>% group_by(tool) %>% summarise(n = n(), min_ptool = min(p_tool), min_pbool = min(p_bool), md_bool = median(p_bool), mn_bool = mean(p_bool))

t60 %>% filter(!n, p_bool> 0.95) %>% ungroup() %>% group_by(tool) %>% summarise(N = sum(N_tool_class_bool), n = n(), min_ptool = min(p_tool), min_pbool = min(p_bool), md_bool = median(p_bool), mn_bool = mean(p_bool))
t60 %>% filter(n, p_bool> 0.95) %>% ungroup() %>% group_by(tool) %>% summarise(N = sum(N_tool_class_bool), n = n(), min_ptool = min(p_tool), min_pbool = min(p_bool), md_bool = median(p_bool), mn_bool = mean(p_bool))
t70 %>% filter(!n, p_bool> 0.95) %>% ungroup() %>% group_by(tool) %>% summarise(N = sum(N_tool_class_bool), n = n(), min_ptool = min(p_tool), min_pbool = min(p_bool), md_bool = median(p_bool), mn_bool = mean(p_bool))
t70 %>% filter(n, p_bool> 0.95) %>% ungroup() %>% group_by(tool) %>% summarise(N = sum(N_tool_class_bool), n = n(), min_ptool = min(p_tool), min_pbool = min(p_bool), md_bool = median(p_bool), mn_bool = mean(p_bool))
t70 %>% ungroup() %>% group_by(n, tool) %>% summarise(N = sum(N_tool_class_bool), classes = n(), min_ptool = min(p_tool), min_pbool = min(p_bool), md_bool = median(p_bool), mn_bool = mean(p_bool))
t70 %>% filter(p_bool> 0.9) %>% ungroup() %>% group_by(n, tool) %>% summarise(N = sum(N_tool_class_bool), classes = n(), min_ptool = min(p_tool), min_pbool = min(p_bool), md_bool = median(p_bool), mn_bool = mean(p_bool))

t50 
t50 %>% filter(p_bool> 0.9) %>% ungroup() %>% group_by(n, tool) %>% summarise(N = sum(N_tool_class_bool), classes = n(), min_ptool = min(p_tool), min_pbool = min(p_bool), md_bool = median(p_bool), mn_bool = mean(p_bool))

t70 %>% 
  ggplot(aes(x = new_level, y = p_bool)) +
  geom_point() + 
  xlab("") + 
  facet_grid(tool ~  . , scales = "free") +
  theme(axis.title.y = element_blank(),
        legend.position = "bottom",
        panel.border = element_rect(fill = "transparent", color = "black", linewidth = 1),
        axis.text.x = element_text(angle = 90, size = general_size),
        axis.text.y = element_text(size = general_size),
        strip.text = element_text(size = general_size, face = "bold"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title = element_text(size = general_size + 2, face = "bold"))


t70_plot <- t70 %>% 
  ggplot(aes(x = new_level, y = N_tool_class_bool)) +
  geom_col(aes(fill=n)) + 
  scale_fill_manual(values = pal15_rep[c(3,13)]) +  
  xlab("") + 
  labs(fill = "Present in other tool") +
  facet_grid(tool ~  . , scales = "free") +
  theme_minimal() +
  theme(legend.position = "bottom",
        axis.title.y = element_blank(),
        panel.border = element_rect(fill = "transparent", color = "black", linewidth = 1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = general_size + 2, face = "bold"),
        axis.title = element_text(size = general_size + 2, face = "bold"),
        axis.text.x = element_text(angle = 90, size = general_size),
        axis.text.y = element_text(size = general_size))

ggsave("~/Documents/plots_project/t70.svg", t70_plot, width = 10, height = 7)


data.frame(hm_recall_fnr_2 %>% filter(tool %in% "DeepARG (nt)") %>% 
  arrange(med_recall) %>% filter(new_level %in% (t70 %>% filter(p_bool> 0.9, !n, tool %in% "DeepARG (nt)") %>% select(new_level) %>% pull())))

data.frame(hm_recall_fnr_2 %>% filter(tool %in% "RGI (DIAMOND - nt)") %>% 
             arrange(med_recall) %>% filter(new_level %in% (t70 %>% filter(p_bool> 0.9, !n, tool %in% "RGI (DIAMOND - nt)") %>% select(new_level) %>% pull())))

data.frame(hm_recall_fnr_2 %>% filter(tool %in% "RGI (DIAMOND - nt)") %>% 
             arrange(med_recall) %>% filter(new_level %in% (t50 %>% filter(p_bool> 0.9, !n, tool %in% "RGI (DIAMOND - nt)") %>% select(new_level) %>% pull())))






####

t50_fa <- unigenes_per_tool %>% 
  group_by(tool) %>% mutate(total_tool = n()) %>% ungroup() %>%
  group_by(tool, new_level) %>% mutate(total_class = n()) %>% ungroup() %>%
  filter(tool %in% c("fARGene (nt)", "AMRFinderPlus (aa)") & id < 50) %>% 
  group_by(query) %>% 
  mutate(n = n_distinct(tool))  %>% 
  mutate(n = n > 1) %>% 
  ungroup() %>% group_by(tool, new_level, n) %>% 
  summarise(total_tool = total_tool[1], total_class = total_class[1], N_tool_class_bool = n()) %>% 
  mutate(p_bool = N_tool_class_bool / sum(N_tool_class_bool)) %>%
  ungroup() %>% group_by(tool) %>% mutate(p_tool = N_tool_class_bool / sum(N_tool_class_bool)) %>%
  mutate(P = N_tool_class_bool / total_class) %>%
  arrange(tool, desc(N_tool_class_bool), desc(P), desc(p_tool), desc(p_bool)) 

t70_fa <- unigenes_per_tool %>% 
  group_by(tool) %>% mutate(total_tool = n()) %>% ungroup() %>%
  group_by(tool, new_level) %>% mutate(total_class = n()) %>% ungroup() %>%
  filter(tool %in% c("fARGene (nt)", "AMRFinderPlus (aa)") & id < 70) %>% 
  group_by(query) %>% 
  mutate(n = n_distinct(tool))  %>% 
  mutate(n = n > 1) %>% 
  ungroup() %>% group_by(tool, new_level, n) %>% 
  summarise(total_tool = total_tool[1], total_class = total_class[1], N_tool_class_bool = n()) %>% 
  mutate(p_bool = N_tool_class_bool / sum(N_tool_class_bool)) %>%
  ungroup() %>% group_by(tool) %>% mutate(p_tool = N_tool_class_bool / sum(N_tool_class_bool)) %>%
  mutate(P = N_tool_class_bool / total_class) %>%
  arrange(tool, desc(N_tool_class_bool), desc(P), desc(p_tool), desc(p_bool)) 

t70_fa %>% filter(p_bool> 0.9) %>% ungroup() %>% group_by(n, tool) %>% summarise(N = sum(N_tool_class_bool), classes = n(), min_ptool = min(p_tool), min_pbool = min(p_bool), md_bool = median(p_bool), mn_bool = mean(p_bool))

t70_fa_plot <- t70_fa %>% 
  ggplot(aes(x = new_level, y = N_tool_class_bool)) +
  geom_col(aes(fill=n)) + 
  scale_fill_manual(values = pal15_rep[c(3,13)]) +  
  xlab("") + 
  labs(fill = "Present in other tool") +
  facet_grid(tool ~  . , scales = "free") +
  theme_minimal() +
  theme(legend.position = "bottom",
        axis.title.y = element_blank(),
        panel.border = element_rect(fill = "transparent", color = "black", linewidth = 1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = general_size + 2, face = "bold"),
        axis.title = element_text(size = general_size + 2, face = "bold"),
        axis.text.x = element_text(angle = 90, size = general_size),
        axis.text.y = element_text(size = general_size))
t70_fa_plot


ggsave("~/Documents/plots_project/t70_fa.svg", t70_fa_plot, width = 10, height = 7)


plot(lst$deeparg.norm$id, lst$deeparg.norm$probability)
plot(lst$deeparg.norm$id, lst$deeparg.norm$alignment.length)
plot(lst$deeparg.norm$id, lst$deeparg.norm$alignment.bitscore)
plot(lst$deeparg.norm$alignment.bitscore, lst$deeparg.norm$id)
plot(lst$deeparg.norm$alignment.length, lst$deeparg.norm$probability)
plot(lst$deeparg.norm$alignment.length, lst$deeparg.norm$alignment.bitscore)
plot(lst$deeparg.norm$alignment.bitscore, lst$deeparg.norm$probability)

sum(lst$deeparg.norm$id < 60 & lst$deeparg.norm$probability>.95)
hist(lst$deeparg.norm$alignment.bitscore[lst$deeparg.norm$id < 60 & lst$deeparg.norm$probability>.95])
hist(lst$deeparg.norm$alignment.length[lst$deeparg.norm$id < 60 & lst$deeparg.norm$probability>.95])





##################################################################################################################################################################
##################################################################################################################################################################
##################################################################################################################################################################
##################################################################################################################################################################
##################################################################################################################################################################
##################################################################################################################################################################


core <- readRDS(file = "code_R_analysis/output_abundance_diversity_resistome/core_resistome.rds")
pan <- readRDS(file = "code_R_analysis/output_abundance_diversity_resistome/pan_resistome.rds")


# abundance %>% filter(!habitat %in% c( "amplicon", "isolate",  "built-environment" )) %>% 
#  select(sample) %>% distinct() %>% mutate(n = n())

core <- core %>% mutate(habitat = factor(habitat, levels = EN))
pan <- pan %>% mutate(habitat = factor(habitat, levels = EN))
pan <- pan %>% mutate(tool = ifelse(tool == "RGI (DIAMOND nt)", "RGI (DIAMOND - nt)", tool))
core <- core %>% mutate(tool = ifelse(tool == "RGI (DIAMOND nt)", "RGI (DIAMOND - nt)", tool))

sumcore <- core %>% filter(cut %in% 0.5 & cnt > 900, !habitat %in% c( "amplicon", "isolate",  "built-environment" ),
                           tool %in% tool_2) %>% ungroup() %>% 
  group_by(new_level, tool, habitat) %>% summarise(unigenes = n_distinct(ARO))  %>% 
  mutate(tool = factor(tool, levels = tool_2))

sumpan <-pan %>% ungroup() %>% group_by(tool, habitat, aggregation, gene_class) %>% 
  summarise(md = median(unigenes), mn = mean(unigenes)) %>%
  filter(aggregation == "new_level", !habitat %in% c( "amplicon", "isolate",  "built-environment" ),
         tool %in% tool_2) %>% mutate(tool = factor(tool, levels = tool_2))

general_size <- 10

sumpan %>% filter(habitat %in% c("human gut")) %>% 
ggplot(aes(x = tool, y = mn)) +
  geom_col(aes(fill = tool)) +
  theme_minimal() +
  labs(fill = "Gene class") +
  #facet_grid(. ~ tool, scales = "free_x") +
  xlab("Tool") +
  xlab("Size pan-resistome") +
  scale_fill_manual(values = c(pal_10_d,pal_10_d,pal_10_d,pal_10_d,pal_10_d,pal_10_d)) +
  theme(axis.title.y = element_blank(),
        #legend.position = "bottom",
        panel.border = element_rect(fill = "transparent", color = "black", linewidth = 1),
        axis.text.x = element_text(angle = 90, size = general_size),
        axis.text.y = element_text(size = general_size),
        strip.text = element_text(size = general_size, face = "bold"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title = element_text(size = general_size + 2, face = "bold"))


sumpan %>% filter(habitat %in% c("wastewater")) %>% 
  ggplot(aes(x = tool, y = mn)) +
  geom_col(aes(fill = gene_class)) +
  theme_minimal() +
  labs(fill = "Gene class") +
  xlab("Tool") +
  xlab("Size pan-resistome") +
  theme(axis.title.y = element_blank(),
        #legend.position = "bottom",
        panel.border = element_rect(fill = "transparent", color = "black", linewidth = 1),
        axis.text.x = element_text(angle = 90, size = general_size),
        axis.text.y = element_text(size = general_size),
        strip.text = element_text(size = general_size, face = "bold"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title = element_text(size = general_size + 2, face = "bold"))


sumcore %>% filter(habitat %in% c("human gut")) %>% 
  ggplot(aes(x = tool, y = unigenes)) +
  geom_col(aes(fill = tool)) +
  theme_minimal() +
  labs(fill = "Gene class") +
  #facet_grid(. ~ tool, scales = "free_x") +
  xlab("Tool") +
  xlab("Core-resistome") +
  scale_fill_manual(values = c(pal_10_d,pal_10_d,pal_10_d,pal_10_d,pal_10_d,pal_10_d)) +
  theme(axis.title.y = element_blank(),
        #legend.position = "bottom",
        panel.border = element_rect(fill = "transparent", color = "black", linewidth = 1),
        axis.text.x = element_text(angle = 90, size = general_size),
        axis.text.y = element_text(size = general_size),
        strip.text = element_text(size = general_size, face = "bold"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title = element_text(size = general_size + 2, face = "bold"))


sumcore %>% filter(habitat %in% c("human gut")) %>% 
  ggplot(aes(x = tool, y = unigenes)) +
  geom_col(aes(fill = new_level)) +
  theme_minimal() +
  labs(fill = "Gene class") +
  #facet_grid(. ~ tool, scales = "free_x") +
  xlab("Tool") +
  xlab("Core-resistome") +
  scale_fill_manual(values = c(pal_10_d,pal_10_d,pal_10_d,pal_10_d,pal_10_d,pal_10_d)) +
  theme(axis.title.y = element_blank(),
        #legend.position = "bottom",
        panel.border = element_rect(fill = "transparent", color = "black", linewidth = 1),
        axis.text.x = element_text(angle = 90, size = general_size),
        axis.text.y = element_text(size = general_size),
        strip.text = element_text(size = general_size, face = "bold"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title = element_text(size = general_size + 2, face = "bold"))


sumcore %>% filter(habitat %in% c("pig gut")) %>% 
  ggplot(aes(x = tool, y = unigenes)) +
  geom_col(aes(fill = new_level)) +
  theme_minimal() +
  labs(fill = "Gene class") +
  xlab("Tool") +
  xlab("Size core-resistome") +
  theme(axis.title.y = element_blank(),
        #legend.position = "bottom",
        panel.border = element_rect(fill = "transparent", color = "black", linewidth = 1),
        axis.text.x = element_text(angle = 90, size = general_size),
        axis.text.y = element_text(size = general_size),
        strip.text = element_text(size = general_size, face = "bold"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title = element_text(size = general_size + 2, face = "bold"))

sumcore %>% filter(habitat %in% c("wastewater")) %>% 
  ggplot(aes(x = tool, y = unigenes)) +
  geom_col(aes(fill = new_level)) +
  theme_minimal() +
  labs(fill = "Gene class") +
  xlab("Tool") +
  xlab("Size core-resistome") +
  theme(axis.title.y = element_blank(),
        #legend.position = "bottom",
        panel.border = element_rect(fill = "transparent", color = "black", linewidth = 1),
        axis.text.x = element_text(angle = 90, size = general_size),
        axis.text.y = element_text(size = general_size),
        strip.text = element_text(size = general_size, face = "bold"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title = element_text(size = general_size + 2, face = "bold"))


##################################################################################################################################################################
##################################################################################################################################################################
##################################################################################################################################################################
##################################################################################################################################################################
##################################################################################################################################################################
##################################################################################################################################################################




fig1_div <- diversity_parent  %>% 
  group_by(habitat, habitat2, tool, sample) %>% summarise(total = sum(unigenes)) %>%
  filter(!habitat %in% not_env, tool %in% tool_selected, !sample %in% extreme_samples) %>% 
  ggplot(aes( x = habitat2)) +
  geom_boxplot(aes(y = total, fill = tool), outlier.shape = NA) +
  scale_fill_manual(values = pal) +
  theme_minimal() +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))  + 
  ylab("Diversity") +
  xlab("Source") +
  labs(fill = "") +
  theme(
    legend.position = "bottom",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(angle = 90))
#scale_y_continuous(labels = scales::log10_trans())
fig1_div


pdf("fig1_diversity.pdf", width = 14, height = 10)
fig1_div
dev.off()

png("fig1_abundance.png",  width= 1600,  height = 800, res = 150)
fig1
dev.off()

png("fig1_diversity.png",  width= 1600,  height = 800, res = 150)
fig1_div
dev.off()

# heatmap

factor_new_level_heatmap <- abundance_parent %>% ungroup() %>% 
  group_by(tool, new_level) %>% summarise(total = sum(normed10m)) %>%
  ungroup() %>% arrange(tool, desc(total)) %>% 
  group_by(tool, new_level) %>%  ungroup() %>% select(new_level) %>% distinct() %>% pull()


hm1 <- abundance_parent %>% filter(!habitat %in% not_env, tool %in% tool_selected, !sample %in% extreme_samples) %>% 
  ungroup() %>% mutate(N = n_distinct(sample)) %>% 
  ungroup() %>% group_by(tool, new_level) %>% summarise(mn = sum(normed10m)/N[1])  %>%
  #ungroup() %>% group_by(tool, new_level) %>% summarise(med = median(mn))  %>%
  ungroup() %>% group_by(tool) %>% mutate(p = log10(mn+1) / max(log10(mn+1))) %>%
  ungroup()  %>% mutate(tool = factor(as.character(tool), levels = tool_selected)) %>%
  ungroup()  %>% complete(tool, new_level, fill = list(p = 0)) %>%
  mutate(new_level = factor(as.character(new_level), levels = factor_new_level_heatmap)) %>% 
   ggplot(aes( x = tool, y = new_level, fill = p)) +
   geom_tile() +
   #scale_fill_gradient(low = "white", high = pal[length(pal)]) +
  scale_fill_viridis_c() + 
#   facet_grid(. ~ habitat) +
  xlab("Tool") + 
  ylab("Ontology") + 
   theme_minimal() +
  labs(fill = "") +
   theme(
     legend.position = "bottom",
     panel.grid.major = element_blank(), 
     panel.grid.minor = element_blank(),
     panel.background = element_blank(),
     axis.text.x = element_text(angle = 90)) 


hm2 <- abundance_parent %>% filter(!habitat %in% not_env, tool %in% tool_selected, !sample %in% extreme_samples) %>% 
  ungroup() %>% mutate(N = n_distinct(sample)) %>% 
  ungroup() %>% group_by(tool, new_level) %>% summarise(mn = sum(normed10m)/N[1])  %>%
  #ungroup() %>% group_by(tool, new_level) %>% summarise(med = median(mn))  %>%
  ungroup() %>% group_by(tool) %>% mutate(p = log10(mn+1) / max(log10(mn+1))) %>%
  ungroup()  %>% mutate(tool = factor(as.character(tool), levels = tool_selected)) %>%
  ungroup()  %>% complete(tool, new_level, fill = list(p = 0)) %>%
  mutate(new_level = factor(as.character(new_level), levels = factor_new_level_heatmap)) %>% 
  ggplot(aes( x = new_level, y = tool, fill = p)) +
  geom_tile() +
  #scale_fill_gradient(low = "white", high = pal[length(pal)]) +
  scale_fill_viridis_c() + 
  #   facet_grid(. ~ habitat) +
  ylab("Tool") + 
  xlab("Ontology") + 
  theme_minimal() +
  labs(fill = "") +
  theme(
    legend.position = "bottom",
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.text.x = element_text(angle = 90)) 

hm3 <- abundance_parent %>% filter(!habitat %in% not_env, tool %in% tool_selected, !sample %in% extreme_samples) %>% 
  ungroup() %>% mutate(N = n_distinct(sample)) %>% 
  ungroup() %>% group_by(tool, new_level) %>% summarise(mn = sum(normed10m)/N[1])  %>%
  #ungroup() %>% group_by(tool, new_level) %>% summarise(med = median(mn))  %>%
  ungroup() %>% group_by(new_level) %>% mutate(p = log10(mn+1) / max(log10(mn+1))) %>%
  ungroup()  %>% mutate(tool = factor(as.character(tool), levels = tool_selected)) %>%
  ungroup()  %>% complete(tool, new_level, fill = list(p = 0)) %>%
  mutate(new_level = factor(as.character(new_level), levels = factor_new_level_heatmap)) %>% 
  ggplot(aes( x = new_level, y = tool, fill = p)) +
  geom_tile() +
  #scale_fill_gradient(low = "white", high = pal[length(pal)]) +
  scale_fill_viridis_c() + 
  #   facet_grid(. ~ habitat) +
  ylab("Tool") + 
  xlab("Ontology") + 
  theme_minimal() +
  labs(fill = "") +
  theme(
    legend.position = "bottom",
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.text.x = element_text(angle = 90)) 


png("hm_by_tool.png",  width= 1600,  height = 800, res = 150)
hm2
dev.off()

png("hm_by_class.png" , width = 1600,  height = 800, res = 150)
hm3
dev.off()








##### abundance Radial per tool

rad0 <- abundance_parent %>% filter(!habitat %in% not_env, tool %in% tool_selected, !sample %in% extreme_samples) %>% 
  ungroup() %>% mutate(N = n_distinct(sample)) %>%  # total number of sampls
  group_by(tool, new_level) %>% summarise(mn = log10(sum(normed10m)/N[1] +1))  %>% # average 
  ungroup() %>% arrange(tool, desc(mn)) %>%
  pivot_longer(-c(tool, new_level), names_to = "variable", values_to = "value") %>%
  ungroup() %>%
  complete(tool, new_level, variable, fill = list(value = 0)) %>%
  group_by(tool, variable) %>%
  mutate(angle = 2 * pi * (row_number() - 1) / n(),  # angle for each variable
         x = sin(angle) * value,
         y = cos(angle) * value)


# adding the first point to close the polygone
df_poly_mn <- rad0 %>% filter(variable %in% "mn") %>%
  group_by(tool, variable) %>%
  do(rbind(., slice(., 1))) %>%
  mutate(group = rep(unique(df$group), each = nrow(rad0) + 1)) 


make_circle <- function(radius, n = 500, center_x = 0, center_y = 0) {
  tibble(
    x = cos(seq(0, 2 * pi, length.out = n)) * radius + center_x,
    y = sin(seq(0, 2 * pi, length.out = n)) * radius + center_y,
    r = radius
  )
}

# Combine all circles into one dataframe
circles_abundance <- bind_rows(lapply(c( 0.5, 1, 1.5, 2, 2.5, 3), make_circle))

df_radial_plot <- df_poly_mn %>% filter( tool %in% c("fargene","deeparg", "rgi.diamond", "resfinder"), value > 0) %>% 
  ungroup() %>% group_by(tool, new_level) %>% 
  mutate(label = if_else(value == max(value) & value > 1, new_level, "")) %>% 
  mutate(hjust_val = ifelse(x < 0, 1, 0)) %>%
  mutate(hjust_val = ifelse(angle %in% c(0, pi / 2, pi, pi*3/2), 0.5, hjust_val)) %>%
  mutate(angle_text = abs(atan2(x,y)* 180/pi) + 90) %>%
  mutate(angle_text = ifelse(x > 0 & y > 0, 180 - angle_text, 
                             ifelse(x > 0 & y < 0, 180 - angle_text, 
                                    ifelse( x < 0 & y > 0,  angle_text - 180, 
                                            ifelse(x < 0 & y < 0, angle_text - 180, angle_text )))))

df_radial_plot1 <- bind_rows(df_radial_plot, df_radial_plot %>% mutate(x = 0, y = 0, label = ""))



fig2_2 <- df_radial_plot1 %>%
  ggplot( aes(x = x, y = y)) +
  geom_polygon(aes(fill = tool), alpha = 0.4) +
  geom_point(aes(color = tool),  size = 2, alpha = 0.6) +
  geom_text(aes(label = label, hjust = hjust_val, angle = angle_text), vjust = -0.5, size = 2) +
  geom_path(data = circles_abundance, aes(x, y, group = r), color = "grey",  show.legend = F, alpha = 0.4) + 
  geom_line(aes(group = new_level), alpha = 0.4) +
  scale_fill_manual(values = pal4) +
  scale_color_manual(values = pal4) +
  coord_equal() +
  facet_grid(. ~ tool) +
  #geom_point(data = data.frame(x = 0, y = 0, tool="a"), aes(x = x, y = y), color = "black", size = 2, show.legend = F) +
  theme_minimal() +
  xlab("") +
  ylab("") +
  labs(fill = "", color="") +
  theme(
    strip.text = element_text(face = "bold"),
    strip.background = element_rect(fill = "lightgrey", color = NA),
    plot.title = element_text(face = "bold"),
    legend.position = "none",
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()) +
  labs(title = "")

fig2_2



tool_name = c("rgi.diamond", "rgi.diamond.prot", "rgi.blast", "deeparg", "deeparg.prot", 
  "fargene", "fargene.prot", "resfinder", "amrfinder", "amrfinder.prot", "abricate.argannot", 
  "abricate.card", "abricate.megares", "abricate.ncbi", "abricate.resfinder")

names(tool_name) <- c("rgi.diamond", "rgi.diamond.prot", "rgi.blast", "deeparg.norm", "deeparg.norm.prot",
                      "fargene", "fargene.prot", "resfinder.norm", "amrfinder.norm", "amrfinder.norm.prot", "abricate.argannot.norm", 
                      "abricate.card.norm", "abricate.megares.norm", "abricate.ncbi.norm", "abricate.resfinder.norm")







rad2 <- diversity_parent %>% filter(!habitat %in% not_env, tool %in% tool_selected, !sample %in% extreme_samples) %>% 
  ungroup() %>% mutate(N = n_distinct(sample)) %>%  # total number of sampls
  group_by(tool, new_level) %>% summarise(mn = log10(sum(unigenes)/N[1] +1))  %>% # average 
  ungroup() %>% arrange(tool, desc(mn)) %>%
  pivot_longer(-c(tool, new_level), names_to = "variable", values_to = "value") %>%
  ungroup() %>%
  complete(tool, new_level, variable, fill = list(value = 0)) %>%
  group_by(tool, variable) %>%
  mutate(angle = 2 * pi * (row_number() - 1) / n(),  # angle for each variable
         x = sin(angle) * value,
         y = cos(angle) * value)


# adding the first point to close the polygone
df_poly_mn_2 <- rad2 %>% filter(variable %in% "mn") %>%
  group_by(tool, variable) %>%
  do(rbind(., slice(., 1))) %>%
  mutate(group = rep(unique(df$group), each = nrow(rad2) + 1)) 



# Combine all circles into one dataframe
circles_diversity <- bind_rows(lapply(c( 0.5, 1, 1.5, 2, 2.5, 3), make_circle))

df_radial_plot_2 <- df_poly_mn_2 %>% filter( tool %in% c("fargene","deeparg", "rgi.diamond", "resfinder"), value > 0) %>% 
  ungroup() %>% group_by(tool, new_level) %>% 
  mutate(label = if_else(value == max(value) & value > 1, new_level, "")) %>% 
  select(-angle) %>%
  left_join(df_radial_plot %>% select(tool, new_level, angle, hjust_val, angle_text), by = c("tool", "new_level"))

df_radial_plot_21 <- bind_rows(df_radial_plot_2, df_radial_plot_2 %>% mutate(x = 0, y = 0, label = ""))



fig3_2 <- df_radial_plot_21 %>%
  ggplot( aes(x = x, y = y)) +
  geom_polygon(aes(fill = tool), alpha = 0.4) +
  geom_point(aes(color = tool),  size = 2, alpha = 0.6) +
  geom_text(aes(label = label, hjust = hjust_val, angle = angle_text), vjust = -0.5, size = 2) +
  geom_path(data = circles_abundance, aes(x, y, group = r), color = "grey",  show.legend = F, alpha = 0.4) + 
  geom_line(aes(group = new_level), alpha = 0.4) +
  scale_fill_manual(values = pal4) +
  scale_color_manual(values = pal4) +
  coord_equal() +
  facet_grid(. ~ tool) +
  #geom_point(data = data.frame(x = 0, y = 0, tool="a"), aes(x = x, y = y), color = "black", size = 2, show.legend = F) +
  theme_minimal() +
  xlab("") +
  ylab("") +
  labs(fill = "", color="") +
  theme(
    strip.text = element_text(face = "bold"),
    strip.background = element_rect(fill = "lightgrey", color = NA),
    plot.title = element_text(face = "bold"),
    legend.position = "none",
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()) +
  labs(title = "")

fig3_2




#abu_div <- abundance_parent %>% select(sample, tool, habitat, habitat2, parent_label, new_level, normed10m) %>% 
#  left_join(diversity_parent %>% select(sample, tool, new_level, habitat, habitat2, parent_label, unigenes)) %>%
#  ungroup() %>% mutate(N = n_distinct(sample)) %>%  group_by(sample, tool, new_level) %>% summarise(N = max(N), normed10m = sum(normed10m), unigenes = sum(unigenes)) %>% 
#  ungroup() %>% mutate(ratio = normed10m  / unigenes) %>% mutate(logratio = log10(ratio + 1))

#abu_div2 <- data.frame(abu_div %>% group_by(tool, new_level) %>% summarise(sr = sum(ratio), slr= sum(logratio), mn = sum(ratio)/max(N), mn_log = sum(logratio)/max(N)) %>% 
#  ungroup() %>% group_by(tool) %>% mutate(normratio = (mn - mean(mn)) /sd(mn), sd = sd(mn), normrlogratio = (mn_log - mean(mn_log)) /sd(mn_log), sdlog = sd(mn_log) ))






g_legend <- function(a.gplot){
  tmp <- ggplotGrob(a.gplot)
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

g_legend(fig2)

fig_ab_div <- grid.arrange(fig2_2 + ggtitle("Abundance") + theme(legend.position = "none"), 
             fig3_2 + ggtitle("Diversity") , 
             layout_matrix = rbind(rep(1, 5), rep(1, 5), rep(1, 5), rep(1, 5), 
                                               rep(2, 5), rep(2, 5), rep(2, 5), rep(2, 5)))



png("fig_ab_div.png", , width = 1600, height = 800, res = 150)
grid.arrange(fig2_2 + ggtitle("Abundance") + theme(legend.position = "none"), 
             fig3_2 + ggtitle("Diversity") , 
             layout_matrix = rbind(rep(1, 5), rep(1, 5), rep(1, 5), rep(1, 5), 
                                   rep(2, 5), rep(2, 5), rep(2, 5), rep(2, 5)))
dev.off()


# bar saggered plot 




########################################################################################
########################################################################################
########################################################################################
########################################################################################





## function

plot_hm_class_tool <- function(PM_new, cl, cl2, ngenes){
  cl3 <- diversity_parent %>% filter(tool %in% cl2) %>% 
    group_by(new_level) %>%
    summarise(m = n_distinct(unigenes)) %>% arrange(desc(m)) %>% 
    ungroup() %>% select(new_level) %>% mutate(as.character(new_level)) %>% pull() 
  cl3 <- as.character(c(cl3[1:ngenes],
                        cl3[(length(cl3)-(ngenes-1)):length(cl3)]))
  
  hmdeep <- ggplot(
    PM_new %>% 
      filter(new_level %in% cl3, 
             tool %in% cl | variable %in% cl) %>%
      mutate(relative = "Tool", new_level = factor(new_level, levels = cl3)) %>%
      mutate(relative = ifelse(variable %in% cl, cl, relative)) %>%
      mutate(x = ifelse(variable %in% cl, as.character(variable), as.character(tool)),
             y = ifelse(variable %in% cl, as.character(tool), as.character(variable))) %>%
      mutate(x = factor(x, levels = levels(tool)), y = factor(y, levels = levels(tool))) %>% 
      mutate(adjusted = ifelse(is.na(adjusted), 0, adjusted)) %>% 
      mutate(text = ifelse(relative == "Tool", col_sum - value, col_sum - value)) %>%
      mutate(text = ifelse(text == 0, "", text)),
    aes(x = new_level, y = y, fill = adjusted*100)) +
    geom_tile(color = "black") +
    scale_fill_gradient(low = "white", high = pal[length(pal)]) +
    geom_text(aes(label = text), color = "black", size = 3) +
    #coord_fixed() +
    ylab("") +
    xlab("") +
    labs(fill = "%") +
    facet_grid(relative ~ new_level, scales = "free_x") +
    theme_minimal() +
    theme(legend.position = "right",
          strip.text = element_text(face = "bold"),
          strip.background = element_rect(fill = "lightgrey", color = NA),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          axis.text.x = element_blank()) 
  return(hmdeep)
}


cl <- "deeparg.norm" #rgi.diamond, fargene
cl2 <- "deeparg"     #rgi.diamond, fargene

plot_hm_class_tool(PM_new, "deeparg.norm", "deeparg", ngenes=6)
plot_hm_class_tool(PM_new, "deeparg.norm", "fargene", ngenes=6)
plot_hm_class_tool(PM_new, "fargene", "fargene", ngenes=6)
plot_hm_class_tool(PM_new, "rgi.diamond", "rgi.diamond", ngenes=6)
plot_hm_class_tool(PM_new, "resfinder.norm", "resfinder", ngenes=6)


png("fargene_hm.png",  width= 2000,  height = 2000, res = 150)
plot_hm_class_tool(PM_new, "fargene", "fargene", ngenes=6)
dev.off()

png("rgi_hm.png",  width= 2000,  height = 1000, res = 150)
plot_hm_class_tool(PM_new, "rgi.diamond", "rgi.diamond", ngenes=6)
dev.off()

png("deeparg_hm.png",  width= 2000,  height = 1000, res = 150)
plot_hm_class_tool(PM_new, "deeparg.norm", "deeparg", ngenes=6)
dev.off()

png("deeparg_with_fargene_genes_hm.png",  width= 2000,  height = 1000, res = 150)
plot_hm_class_tool(PM_new, "deeparg.norm", "fargene", ngenes=6)
dev.off()



#########################################################################
#########################################################################
#########################################################################
#########################################################################


id.tools <- tools_per_unigene %>% select(tool, query, new_level, id) %>% distinct()
for(j in unique(id.tools$tool)){
  id.tools[,j] <- id.tools$query %in% id.tools$query[id.tools$tool==j]
}


ggplot(id.tools %>% filter(new_level %in% fg_classes, !fargene), aes(id, colour = tool, linetype = fargene)) +
  stat_ecdf() +
  facet_wrap(. ~ new_level, scales = "free_y")

ggplot(id.tools %>% filter(new_level %in% fg_classes, !fargene), aes(id, colour = tool, linetype = fargene)) +
  geom_freqpoly(binwidth = 1) +
  facet_wrap(. ~ new_level, scales = "free_y")


id.tools <- id.tools %>% mutate(tool = factor(tool_name[tool], levels = tool_selected))



plot_positive_id <- function(id.tools, cl, g, ngenes){
  cl3 <- diversity_parent %>% filter(tool %in% cl) %>% 
    group_by(new_level) %>%
    summarise(m = n_distinct(unigenes)) %>% arrange(desc(m)) %>% 
    ungroup() %>% select(new_level) %>% mutate(as.character(new_level)) %>% pull() 
  cl3 <- as.character(unique(cl3[1:ngenes]))
  
  ggplot(id.tools %>% 
           #mutate(fargene = ifelse(tool %in% "fargene", !fargene, fargene)) %>% 
           filter(new_level %in% cl3, !!sym(g)) %>% mutate(new_level = factor(new_level, levels = cl3)) %>%
           ungroup() %>% group_by(tool, new_level, !!sym(g)) %>% arrange(id), aes( id, colour = tool)) +
    geom_freqpoly(binwidth = 5, alpha = 0.7, linewidth = 1.5) +
    facet_wrap(. ~ new_level, scales = "free_y", nrow = 2) +  
    theme_minimal() +
    xlab("Identity level") +
    ylab("Unigenes") +
    labs(color="Tool") +
    theme(
      strip.text = element_text(face = "bold"),
      strip.background = element_rect(fill = "lightgrey", color = NA),
      plot.title = element_text(face = "bold"),
      legend.position = "bottom")
}



plot_negative_id <- function(id.tools, cl, cl2, g, ngenes){

  cl3 <- diversity_parent %>% filter(tool %in% cl) %>% 
    group_by(new_level) %>%
    summarise(m = n_distinct(unigenes)) %>% arrange(desc(m)) %>% 
    ungroup() %>% select(new_level) %>% mutate(as.character(new_level)) %>% pull() 
  cl3 <- as.character(unique(cl3[1:ngenes]))
  
  dftemp <- id.tools %>% 
    mutate(gr = ifelse(tool %in% cl2, !!sym(g), !!sym(g))) %>% 
    filter(new_level %in% cl3, !gr) %>% mutate(new_level = factor(new_level, levels = cl3)) %>%
    ungroup() %>% group_by(tool, new_level, gr) %>% arrange(id)
  
  ggplot(dftemp, aes( id, colour = tool)) +
    geom_freqpoly(binwidth = 5, alpha = 0.7,  linewidth = 1.5) +
    scale_color_discrete(drop = FALSE) +
    facet_wrap(. ~ new_level, scales = "free_y", nrow = 2) +  
    theme_minimal() +
    xlab("Identity level") +
    ylab("Unigenes") +
    labs(color="Tool") +
    theme(
      strip.text = element_text(face = "bold"),
      strip.background = element_rect(fill = "lightgrey", color = NA),
      plot.title = element_text(face = "bold"),
      legend.position = "bottom")
}

fgpos <- plot_positive_id(id.tools, "fargene", "fargene", 6)
fgneg <- plot_negative_id(id.tools, "fargene", "fargene", "fargene", 6)

rgipos <- plot_positive_id(id.tools, "rgi.diamond", "rgi.diamond", 6)
rgineg <- plot_negative_id(id.tools, "rgi.diamond", "rgi.diamond", "rgi.diamond", 6)

rgipos_fg <- plot_positive_id(id.tools, "fargene", "rgi.diamond", 6)
rgineg_fg <- plot_negative_id(id.tools, "fargene", "rgi.diamond", "rgi.diamond", 6)

deepargpos <- plot_positive_id(id.tools,  "fargene", "deeparg.norm", 6)
deepargneg <-plot_negative_id(id.tools, "fargene","deeparg", "deeparg.norm", 6)

deepargpos_fg <- plot_positive_id(id.tools,  "fargene", "deeparg.norm", 6)
deepargneg_fg <-plot_negative_id(id.tools, "fargene","deeparg", "deeparg.norm", 6)

deepargpos <- plot_positive_id(id.tools, "deeparg", "deeparg.norm", 6)
deepargneg <-plot_negative_id(id.tools, "deeparg", "deeparg", "deeparg.norm", 6)

rfpos <- plot_positive_id(id.tools, "resfinder", "resfinder.norm", 6)
rfneg <- plot_negative_id(id.tools, "resfinder", "resfinder", "resfinder.norm", 6)


he = 500
wi = 900

png("fargene_pos.png",  width= wi,  height = he, res = 150)
fgpos
dev.off()

png("fargene_neg.png",  width= wi,  height = he, res = 150)
fgneg
dev.off()

png("rgi_pos.png",  width= wi,  height = he, res = 150)
rgipos
dev.off()

png("rgi_neg.png",  width= wi,  height = he, res = 150)
rgineg
dev.off()

png("deep_pos.png",  width= wi,  height = he, res = 150)
deepargpos
dev.off()

png("deep_neg.png",  width= wi,  height = he, res = 150)
deepargneg
dev.off()

png("deep_pos_with_fg_gene.png",  width= wi,  height = he, res = 150)
deepargpos_fg
dev.off()

png("deep_neg_with_fg_gene.png.png",  width= wi,  height = he, res = 150)
deepargneg_fg
dev.off()





###

lst$rgi.diamond$new_level <- new_level_df$new[match(lst$rgi.diamond$parent_description, new_level_df$old)]
ggplot(lst$rgi.diamond %>% filter(new_level %in% c("GPA", "TET - RPG")),
       aes(Best_Hit_Bitscore, colour = new_level)) +
  geom_freqpoly(binwidth = 5, alpha = 0.7,  linewidth = 1.5) +
  scale_color_discrete(drop = FALSE) +
  facet_wrap(. ~ new_level, scales = "free_y", nrow = 2) +  
  scale_x_continuous(breaks = seq(0,1200, 100)) +
  theme_minimal() +
  xlab("Bitscore") +
  ylab("Unigenes") +
  labs(color="") +
  theme(
    strip.text = element_text(face = "bold"),
    strip.background = element_rect(fill = "lightgrey", color = NA),
    plot.title = element_text(face = "bold"),
    legend.position = "none")
  

ggplot(lst$rgi.diamond,
       aes(x = parent_description, fill = new_level)) +
  geom_bar() +
  scale_color_discrete(drop = FALSE) +
  facet_grid(. ~ new_level, scales = "free_x") +  
  theme_minimal() +
  xlab("") +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))  +
  ylab("Unigenes") +
  labs(color="") +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(angle = 90))



lst$deeparg.norm$new_level <- new_level_df$new[match(lst$deeparg.norm$parent_description, new_level_df$old)]
ggplot(lst$deeparg.norm,
       aes(x = parent_description, fill = new_level)) +
  geom_bar() +
  scale_color_discrete(drop = FALSE) +
  facet_grid(. ~ new_level, scales = "free_x") +  
  theme_minimal() +
  xlab("") +
  #scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
  #              labels = scales::trans_format("log10", scales::math_format(10^.x)))  +
  ylab("Unigenes") +
  labs(color="") +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(angle = 90))


lst$deeparg.norm %>% group_by(new_level, ARG.class) %>% summarise(n = n()) %>% arrange(desc(n))
##
##



saveRDS(df2, file = "~/df2.rds")
unique(df2$Parent_Label)




#########################################
#########################################


##### abundance Radial per tool

rad0_env <- abundance_parent %>% filter(habitat %in% c("human gut", "wastewater"), tool %in% c("deeparg","rgi.diamond","fargene"), !sample %in% extreme_samples) %>% 
  ungroup() %>% group_by(habitat) %>% mutate(N = n_distinct(sample)) %>%  # total number of sampls
  ungroup() %>% group_by(tool, habitat, new_level) %>% summarise(mn = log10(sum(normed10m)/N[1] +1))  %>% # average 
  ungroup() %>% arrange(tool, desc(mn)) %>%
  pivot_longer(-c(tool, habitat, new_level), names_to = "variable", values_to = "value") %>%
  ungroup() %>%
  complete(tool, habitat, new_level, variable, fill = list(value = 0)) %>%
  group_by(tool, habitat, variable) %>%
  mutate(angle = 2 * pi * (row_number() - 1) / n(),  # angle for each variable
         x = sin(angle) * value,
         y = cos(angle) * value)


# adding the first point to close the polygone
df_poly_mn_env <- rad0_env %>% filter(variable %in% "mn") %>%
  group_by(tool, habitat, variable) %>%
  do(rbind(., slice(., 1))) %>%
  mutate(group = rep(unique(df$group), each = nrow(rad0) + 1)) 

# Combine all circles into one dataframe
circles_abundance <- bind_rows(lapply(c( 0.5, 1, 1.5, 2, 2.5, 3), make_circle))

df_radial_plot_env <- df_poly_mn_env %>% filter( value > 0) %>% 
  ungroup() %>% group_by(tool, habitat, new_level) %>% 
  mutate(label = if_else(value == max(value) & value > 1, new_level, "")) %>% 
  mutate(hjust_val = ifelse(x < 0, 1, 0)) %>%
  mutate(hjust_val = ifelse(angle %in% c(0, pi / 2, pi, pi*3/2), 0.5, hjust_val)) %>%
  mutate(angle_text = abs(atan2(x,y)* 180/pi) + 90) %>%
  mutate(angle_text = ifelse(x > 0 & y > 0, 180 - angle_text, 
                             ifelse(x > 0 & y < 0, 180 - angle_text, 
                                    ifelse( x < 0 & y > 0,  angle_text - 180, 
                                            ifelse(x < 0 & y < 0, angle_text - 180, angle_text )))))

df_radial_plot1_env <- bind_rows(df_radial_plot_env, df_radial_plot_env %>% mutate(x = 0, y = 0, label = ""))



fig2_2_env <- df_radial_plot1_env %>%
  ggplot( aes(x = x, y = y)) +
  geom_polygon(aes(fill = tool), alpha = 0.4) +
  geom_point(aes(color = tool),  size = 2, alpha = 0.6) +
  geom_text(aes(label = label, hjust = hjust_val, angle = angle_text), vjust = -0.5, size = 2) +
  geom_path(data = circles_abundance, aes(x, y, group = r), color = "grey",  show.legend = F, alpha = 0.4) + 
  geom_line(aes(group = new_level), alpha = 0.4) +
  scale_fill_manual(values = pal4) +
  scale_color_manual(values = pal4) +
  coord_equal() +
  facet_grid(habitat ~ tool) +
  #geom_point(data = data.frame(x = 0, y = 0, tool="a"), aes(x = x, y = y), color = "black", size = 2, show.legend = F) +
  theme_minimal() +
  xlab("") +
  ylab("") +
  labs(fill = "", color="") +
  theme(
    strip.text = element_text(face = "bold"),
    strip.background = element_rect(fill = "lightgrey", color = NA),
    plot.title = element_text(face = "bold"),
    legend.position = "none",
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()) +
  labs(title = "")

fig2_2_env




rad2_env <- diversity_parent %>% filter(habitat %in% c("human gut", "wastewater"), tool %in% c("deeparg","rgi.diamond","fargene"), !sample %in% extreme_samples) %>% 
  ungroup() %>% group_by(habitat) %>% mutate(N = n_distinct(sample)) %>%  # total number of sampls
  group_by(tool, habitat, new_level) %>% summarise(mn = log10(sum(unigenes)/N[1] +1))  %>% # average 
  ungroup() %>% arrange(tool, desc(mn)) %>%
  pivot_longer(-c(tool, habitat, new_level), names_to = "variable", values_to = "value") %>%
  ungroup() %>%
  complete(tool, habitat, new_level, variable, fill = list(value = 0)) %>%
  group_by(tool, habitat, variable) %>%
  mutate(angle = 2 * pi * (row_number() - 1) / n(),  # angle for each variable
         x = sin(angle) * value,
         y = cos(angle) * value)


# adding the first point to close the polygone
df_poly_mn_2_env <- rad2_env %>% filter(variable %in% "mn") %>%
  group_by(tool, habitat, variable) %>%
  do(rbind(., slice(., 1))) %>%
  mutate(group = rep(unique(df$group), each = nrow(rad2) + 1)) 


df_radial_plot_21_env <- df_poly_mn_2_env %>% filter(  value > 0) %>% 
  ungroup() %>% group_by(tool, habitat, new_level) %>% 
  mutate(label = if_else(value == max(value) & value > 1, new_level, "")) %>% 
  select(-angle) %>%
  left_join(df_radial_plot_env %>% select(tool, habitat, new_level, angle, hjust_val, angle_text), by = c("tool", "habitat", "new_level"))

df_radial_plot_21_env <- bind_rows(df_radial_plot_21_env, df_radial_plot_21_env %>% mutate(x = 0, y = 0, label = ""))



fig3_2_env <- df_radial_plot_21_env %>%
  ggplot( aes(x = x, y = y)) +
  geom_polygon(aes(fill = tool), alpha = 0.4) +
  geom_point(aes(color = tool),  size = 2, alpha = 0.6) +
  geom_text(aes(label = label, hjust = hjust_val, angle = angle_text), vjust = -0.5, size = 2) +
  geom_path(data = circles_abundance, aes(x, y, group = r), color = "grey",  show.legend = F, alpha = 0.4) + 
  geom_line(aes(group = new_level), alpha = 0.4) +
  scale_fill_manual(values = pal4) +
  scale_color_manual(values = pal4) +
  coord_equal() +
  facet_grid(habitat ~ tool) +
  #geom_point(data = data.frame(x = 0, y = 0, tool="a"), aes(x = x, y = y), color = "black", size = 2, show.legend = F) +
  theme_minimal() +
  xlab("") +
  ylab("") +
  labs(fill = "", color="") +
  theme(
    strip.text = element_text(face = "bold"),
    strip.background = element_rect(fill = "lightgrey", color = NA),
    plot.title = element_text(face = "bold"),
    legend.position = "none",
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()) +
  labs(title = "")

fig3_2_env




#################
#################
# Venn diagrams


rgi.ven = ggVennDiagram(list(FNA = lst$rgi.diamond$query,
                             FNA.blast = lst$rgi.blast$query,
                             FAA = lst$rgi.diamond.prot$query),
                        color = 1, lwd = 0.7, label_size = 3, set_size = 3) + 
  scale_fill_gradient(low = "#F4FAFE", high = pal[length(pal)]) +
  theme(legend.position = "none") +
  ggtitle("RGI")

rgi.diamond.ven = ggVennDiagram(list(FNA = lst$rgi.diamond$query,
                                     FAA = lst$rgi.diamond.prot$query),
                                color = 1, lwd = 0.7, label_size = 3, set_size = 3) + 
  scale_fill_gradient(low = "#F4FAFE", high = pal[length(pal)]) +
  theme(legend.position = "none") +
  ggtitle("RGI")

deeparg.ven = ggVennDiagram(list(FNA = lst$deeparg.norm$query,
                                 FAA = lst$deeparg.norm.prot$query),
                            color = 1, lwd = 0.7, label_size = 3, set_size = 3) + 
  scale_fill_gradient(low = "#F4FAFE", high = pal[length(pal)]) +
  theme(legend.position = "none") +
  ggtitle("deepARG")

fargene.ven = ggVennDiagram(list(FNA = lst$fargene$query,
                                 FAA = lst$fargene.prot$query),
                            color = 1, lwd = 0.7, label_size = 3, set_size = 3) + 
  scale_fill_gradient(low = "#F4FAFE", high = pal[length(pal)]) +
  theme(legend.position = "none") +
  ggtitle("fARGene")
fargene.ven

amrfinder.ven = ggVennDiagram(list(FNA = lst$amrfinder.norm$query,
                                   FAA = lst$amrfinder.norm.prot$query),
                              color = 1, lwd = 0.7, label_size = 3, set_size = 3) + 
  scale_fill_gradient(low = "#F4FAFE", high = pal[length(pal)]) +
  theme(legend.position = "none") +
  ggtitle("amrfinder")

abricate.ven = ggVennDiagram(list(argannot = lst$abricate.argannot.norm$query,
                                  card = lst$abricate.card.norm$query,
                                  megares = lst$abricate.megares.norm$query,
                                  ncbi = lst$abricate.ncbi.norm$query,
                                  resfinder = lst$abricate.resfinder.norm$query),
                             color = 1, lwd = 0.7, label_size = 3, set_size = 3) + 
  scale_fill_gradient(low = "#F4FAFE", high = pal[length(pal)]) +
  theme(legend.position = "none") +
  ggtitle("abricate")



grid.arrange(rgi.ven, deeparg.ven, fargene.ven, amrfinder.ven, nrow = 2)

png("venn_fargene.png",  width= wi,  height = he, res = 150)
fargene.ven
dev.off()

png("venn_rgi.png",  width= wi,  height = he, res = 150)
rgi.ven
dev.off()

png("venn_deeparg.png",  width= wi,  height = he, res = 150)
deeparg.ven
dev.off()

png("venn_amrfinder.png",  width= wi,  height = he, res = 150)
amrfinder.ven
dev.off()

png("venn_abricate.png",  width= wi,  height = he, res = 150)
abricate.ven
dev.off()
