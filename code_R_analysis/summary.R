library(dplyr)
library(ggplot2)
library(ggVennDiagram)
library(gridExtra)
library(tidyverse)
library(RColorBrewer)

setwd("~/Documents/GitHub/arg_compare/")
df2 <- readRDS(file = "code_R_analysis/output_abundance_diversity_resistome/conversion_ARO_parent_new_level.rds")
lst <- readRDS("code_R_analysis/output_abundance_diversity_resistome/results_tools.rds")


abundance <- readRDS("code_R_analysis/output_abundance_diversity_resistome/abundance_diversity.rds")

#core <- readRDS("code_R_analysis/output_abundance_diversity_resistome/core_resistome.rds")
#pan <- readRDS("code_R_analysis/output_abundance_diversity_resistome/pan_resistome.rds")

#lapply(lst, function(x) sum(!x$query %in% genes_prot_dna$x))

#pal <- c("#543005", "#8c510a", "#bf812d", "#dfc27d", "#f6e8c3", "#f5f5f5", "#c7eae5", "#80cdc1", "#35978f", "#01665e", "#003c30")
pal_12 <- c("#a6cee3", "#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a","#ffff99","#b15928")
pal_8 <-  c("#a6cee3", "#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00")
pal_10 <-  c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a")


pal_10_q <- brewer.pal(8, "Dark2")
pal_10_d <- brewer.pal(10, "BrBG")
pal_9_d <- brewer.pal(9, "BrBG")


#pal_latent <- c("#00B5E2", "#013B67", "#CC0C01", "#985428", "#73006D", "#E08728", "#377E60", "#FD8CC0")


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



abundance <- abundance %>% mutate( tool = 
                     ifelse(tool == "RGI (BLAST nt)", "RGI (BLAST - nt)",
                     ifelse(tool == "RGI (BLAST aa)", "RGI (BLAST - aa)",
                     ifelse(tool == "RGI (DIAMOND aa)", "RGI (DIAMOND - aa)",
                     ifelse(tool == "RGI (DIAMOND nt)", "RGI (DIAMOND - nt)", tool)))))

lst$rgi.blast <- lst$rgi.blast %>% mutate(tool =  "RGI (BLAST - nt)")
lst$rgi.diamond <- lst$rgi.diamond %>% mutate(tool =  "RGI (DIAMOND - nt)")
lst$rgi.diamond.prot <- lst$rgi.diamond.prot %>% mutate(tool =  "RGI (DIAMOND - aa)")

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

#abundance_parent$new_level <- factor(abundance_parent$new_level, levels = factor_new_level2)
#diversity_parent$new_level <- factor(diversity_parent$new_level, levels = factor_new_level2)


# factor for habitat

abundance$habitat <- factor(abundance$habitat, levels = EN)
abundance$habitat <- as.character(abundance$habitat)

# environments that we are not interested in
not_env <- c("amplicon", "isolate")

# SUMMARIES
# SUMMARIES
# SUMMARIES


# total numbrer of unigenes captured with abundance 

args_abundances <- read.delim("data/abundances/args_abundances.tsv")

metadata <- read.delim("data/metadata_GMGC10.sample.meta.tsv")

args_abundances <- args_abundances %>% mutate(habitat = metadata$habitat[match(sample, metadata$sample_id)])

# name of unigenes that had abundance
unigenes_with_abundance <- args_abundances  %>% filter(!habitat %in% not_env) %>%
  select(X) %>% distinct()

length(unigenes_with_abundance$X)

# unigenes captured with the tools
unigenes <- do.call(rbind, lapply(lst, function(x) x[,c("query","tool")])) 

# total number of unique unigenes captured with the tools
length(unique(unigenes$query))


unigenes %>% group_by(tool) %>% summarise(n = n_distinct(query)) %>% ungroup() %>% arrange(n)

do.call(rbind, lapply(lst, function(x) x[,c("tool","ARO")])) %>% 
  group_by(tool) %>% mutate(N = n()) %>% filter(ARO == "") %>% mutate(n = n()) %>% 
  summarise( p = n[1] / N[1])

## UNIGENES missing aro
do.call(rbind, lapply(lst, function(x) x[,c("tool","query","ARO")])) %>% 
  group_by(tool) %>%
  mutate(N = n())  %>% 
  filter(ARO == "") %>% mutate(n = n()) %>% 
  summarise( p = n[1] / N[1])

## NUMBER OF PARENT CLASSES 
do.call(rbind, lapply(lst, function(x) x[,c("tool","parent_description")])) %>% 
  summarise(N = n_distinct(parent_description))

## NUMBER OF new levels
do.call(rbind, lapply(lst, function(x) x[,c("tool","new_level")])) %>% 
  summarise(N = n_distinct(new_level))



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

proportion_new_level_tool <- do.call(rbind, lapply(lst, function(x) x[,c("query","tool","new_level")])) %>% 
  mutate(N = n_distinct(query)) %>% 
  group_by(new_level) %>% mutate(M = n_distinct(query)) %>% ungroup() %>%
  group_by(tool) %>% mutate(ntool = n_distinct(query)) %>% ungroup() %>%
  mutate(P = M / N) %>%
  group_by(new_level, tool) %>% summarise(N = N[1], M = M[1], P = P[1], ntool = ntool[1], n = n_distinct(query)) %>% 
  arrange(desc(n)) %>% mutate(p = n/ntool)


# sum of propotions for GPA AND EFFLUX
proportion_new_level_tool %>% filter(new_level %in% c("GPA","Efflux p.")) %>% group_by(tool) %>% summarise(N = sum(N), M = sum(M), ntool = sum(ntool), n = sum(n), p = sum(p))

# sum of propotions for BETA-LACTAM
proportion_new_level_tool %>% filter(new_level %in% c("Class A", "Class B", "Class D", "Class C")) %>% group_by(tool) %>% summarise(N = sum(N), M = sum(M), P = sum(P), ntool = sum(ntool), n = sum(n), p = sum(p))


do.call(rbind, lapply(lst, function(x) x[,c("query","tool","new_level")])) %>% 
  mutate(N = n_distinct(query)) %>% 
  group_by(new_level) %>% summarise(N = N[1], n = n_distinct(query)) %>% 
  arrange(desc(n)) %>% mutate(p = n/N)

do.call(rbind, lapply(lst, function(x) x[is.na(x$new_level),c("query","tool","new_level")])) %>%
  group_by(tool) %>% summarise(n = n())

### conformity (overlap)


unigenes_per_tool <- do.call(rbind, lapply(seq_along(1:length(lst)), function(j) data.frame(tool = lst[[j]]$tool, query = lst[[j]]$query, 
                                                                                            ARO = lst[[j]]$ARO, parent = lst[[j]]$parent, 
                                                                                            parent_description = lst[[j]]$parent_description,
                                                                                            id = lst[[j]]$id, new_level = lst[[j]]$new_level)))
#unigenes_per_tool <- unigenes_per_tool #%>%
  #mutate(tool = ifelse(tool == "DeepARG (a.a.)", "DeepARG (aa)",
  #ifelse(tool == "RGI (BLAST a.a.)", "RGI (BLAST - aa)",
  #ifelse(tool == "RGI (DIAMOND a.a.)", "RGI (DIAMOND - aa)",
  #ifelse(tool == "RGI (DIAMOND  nt)", "RGI (DIAMOND - nt)",
  #ifelse(tool == "fARGene (a.a.)", "fARGene (aa)",
  #ifelse(tool == "AMRFinderPlus (a.a.)", "AMRFinderPlus (aa)", tool)))))))

#unigenes_per_tool <- unigenes_per_tool %>% filter(!is.na(parent)) %>% distinct()


tool_2 <- c("DeepARG (nt)", "RGI (DIAMOND - nt)", "fARGene (nt)", "AMRFinderPlus (aa)", "ResFinder (nt)",
            "ABRicate (ARG-ANNOT - nt)", "ABRicate (CARD - nt)", "ABRicate (MEGARES - nt)", "ABRicate (NCBI - nt)",
            "ABRicate (ResFinder - nt)")

tools_per_unigene <- unigenes_per_tool %>% ungroup() %>% filter(tool %in% tool_2) %>% 
  arrange(query) %>% 
  group_by(query) %>% 
  mutate(n_tools = n()) %>% 
  mutate(single = (n_tools ==1)) 

alt_name_tools <- c("DeepARG_nt", "RGI_DIAMOND_nt", "fARGene_nt", "AMRFinderPlus_aa", "ResFinder_nt",
                    "ABRicate_ARGANNOT_nt", "ABRicate_CARD_nt", "ABRicate_MEGARES_nt", "ABRicate_NCBI_nt",
                    "ABRicate_ResFinder_nt")

names(alt_name_tools) <- tool_2

tools_per_unigene <- tools_per_unigene %>% mutate(tool = alt_name_tools[tool])

gene.tool2 <- tools_per_unigene %>% select(tool, query, new_level) %>% distinct() %>% 
  mutate(v = T) %>% group_by(tool, new_level) %>% pivot_wider(names_from = tool, values_from = v) 

pairwise_match <- combn(colnames(gene.tool2[,-c(1,2)]), 2, function(cols) {
  sum(gene.tool2[[cols[1]]] & gene.tool2[[cols[2]]], na.rm = TRUE)
})

PM <- matrix(0, ncol(gene.tool2[,-c(1,2)]), ncol(gene.tool2[,-c(1,2)]),
             dimnames = list(colnames(gene.tool2[,-c(1,2)]), colnames(gene.tool2[,-c(1,2)])))
PM <- data.frame(PM)
PM <- PM[0,]
PM$tool <- NULL
PM$new_level <- NULL

for(j in unique(gene.tool2$new_level)){
  pairwise_match <- combn(colnames(gene.tool2[gene.tool2$new_level == j,-c(1,2)]), 2, function(cols) {
    sum(gene.tool2[gene.tool2$new_level ==j, cols[1]] & gene.tool2[gene.tool2$new_level ==j, cols[2]], na.rm = TRUE)
  })
  pairwise_matrix <- matrix(0, ncol(gene.tool2[,-c(1,2)]), ncol(gene.tool2[,-c(1,2)]),
                            dimnames = list(colnames(gene.tool2[,-c(1,2)]), colnames(gene.tool2[,-c(1,2)])))
  pairwise_matrix[lower.tri(pairwise_matrix)] <- pairwise_match
  pairwise_matrix[upper.tri(pairwise_matrix)] <- t(pairwise_matrix)[upper.tri(pairwise_matrix)]
  pairwise_matrix <- data.frame(pairwise_matrix)
  pairwise_matrix$tool <- rownames(pairwise_matrix)
  rownames(pairwise_matrix) <- NULL
  pairwise_matrix$new_level <- j
  PM <- rbind.data.frame(PM, pairwise_matrix)
}


PM_long <- PM %>% 
  pivot_longer(
    cols = -c(tool, new_level),
    names_to = "variable",
    values_to = "value"
  ) 


row_sums <- PM_long %>% ungroup() %>%
  group_by(new_level, tool) %>%
  summarise(row_sum = sum(value))


row_sums <- tools_per_unigene %>% group_by(tool, new_level) %>% 
  summarise(row_sum = n_distinct(query))

col_sums <- tools_per_unigene %>% group_by(tool, new_level) %>% 
  summarise(col_sum = n_distinct(query)) %>% 
  rename(variable = tool)

PM_long <- PM_long %>% mutate(tool = factor(tool, levels = unique(PM_long$tool)))
PM_long <- PM_long %>% mutate(variable = factor(variable, levels = unique(PM_long$tool)))
row_sums <- row_sums %>% ungroup() %>% mutate(tool = factor(tool, levels = unique(PM_long$tool)))
col_sums <- col_sums %>% ungroup() %>% mutate(variable = factor(variable, levels = unique(PM_long$tool)))

PM_new <- PM_long %>%
  left_join(row_sums, by = c("new_level", "tool")) %>%
  left_join(col_sums, by = c("new_level", "variable")) %>%
  mutate(
    i = as.numeric(tool),
    j = as.numeric(variable),
    adjusted = case_when(
      i < j ~ value / col_sum,     # upper triangle
      i > j ~ value / col_sum,     # lower triangle
      TRUE ~ NA                 # diagonal stays the same
    )
  )


PM_new %>% filter(new_level %in% "GPA", tool %in% c("RGI_DIAMOND_nt", "DeepARG_nt"), variable %in% c("RGI_DIAMOND_nt", "DeepARG_nt"))
PM_new %>% filter(new_level %in% "Efflux p.", tool %in% c("RGI_DIAMOND_nt", "DeepARG_nt"), variable %in% c("RGI_DIAMOND_nt", "DeepARG_nt"))

overlap_tool <- PM_new %>% group_by(tool, variable) %>% summarise(value = sum(value, na.rm = T), 
                                                  row_sum = sum(row_sum, na.rm = T),
                                                  col_sum = sum(col_sum, na.rm = T), i = i[1], j = j[1])
overlap_tool <- overlap_tool %>% 
  filter(tool != variable) %>%
  mutate(p1 = value/row_sum, p2 = value/col_sum, p3 = value / (row_sum + col_sum - value))

# jaccard index
JI <- overlap_tool %>% ungroup() %>% group_by(tool) %>% arrange(tool, desc(p3)) %>% select(tool, variable, p3) %>%
  mutate( pair = paste(pmin(as.character(tool), as.character(variable)), 
                       pmax(as.character(tool), as.character(variable)))) %>%
  group_by(pair) %>% slice_head(n = 1) %>% ungroup() %>% arrange(desc(p3) )

data.frame(JI)


# mean and median recall and FNR
PM_new %>% filter(tool != variable) %>% ungroup() %>% filter( !is.na(col_sum)) %>%
  mutate(recall = value / col_sum, fnr = (row_sum - value) / row_sum ) %>% 
  mutate(recall = ifelse(is.na(row_sum), NA, recall), fnr = ifelse(is.na(row_sum), NA, fnr)) %>%
  group_by(tool, new_level) %>% 
  summarise(mn_recall = mean(recall, na.rm = T), med_recall = median(recall, na.rm = T),
            mn_fnr = mean(fnr, na.rm = T), med_fnr = median(fnr, na.rm = T)) %>%
  ungroup() %>%
  group_by(tool) %>% 
  summarise(MN_recall = mean(mn_recall, na.rm = T), MED_recall = median(med_recall, na.rm = T),
            MN_fnr = mean(mn_fnr, na.rm = T), MED_fnr = median(med_fnr, na.rm = T))



# fargene mean and median recall and FNR
PM_new %>% filter(tool != variable) %>% ungroup() %>% filter( !is.na(col_sum)) %>%
  mutate(recall = value / col_sum, fnr = (row_sum - value) / row_sum ) %>% 
  mutate(recall = ifelse(is.na(row_sum), NA, recall), fnr = ifelse(is.na(row_sum), NA, fnr)) %>%
  group_by(tool, new_level) %>% 
  summarise(mn_recall = mean(recall, na.rm = T), med_recall = median(recall, na.rm = T),
            mn_fnr = mean(fnr, na.rm = T), med_fnr = median(fnr, na.rm = T)) %>%
  ungroup() %>% filter(tool == "fARGene_nt" & !is.na(mn_recall))

PM_new %>% filter(tool != variable) %>% ungroup() %>% filter( !is.na(col_sum)) %>%
  mutate(recall = value / col_sum, fnr = (row_sum - value) / row_sum ) %>% 
  mutate(recall = ifelse(is.na(row_sum), NA, recall), fnr = ifelse(is.na(row_sum), NA, fnr)) %>%
  group_by(tool, new_level) %>% 
  summarise(mn_recall = mean(recall, na.rm = T), med_recall = median(recall, na.rm = T),
            mn_fnr = mean(fnr, na.rm = T), med_fnr = median(fnr, na.rm = T)) %>%
  ungroup() %>% filter(tool == "AMRFinderPlus_aa" & !is.na(mn_recall))


as.data.frame(table(lst$amrfinder.norm.prot$new_level, lst$amrfinder.norm.prot$Method)) %>% 
  pivot_wider(names_from = Var2, values_from = Freq) %>% 
  arrange(desc(HMM))

data.frame(hm_recall_fnr_1 %>% filter(tool %in% "AMRFinderPlus_aa"  & value > 0) %>% 
             group_by(tool, new_level) %>% summarise(M_R = median(recall)))

data.frame(hm_recall_fnr_1 %>% filter(tool %in% "AMRFinderPlus_aa" & !variable %in% "fARGene_nt" & value > 0) %>% 
  group_by(tool, new_level) %>% summarise(M_R = median(recall)))

hm_recall_fnr_1 %>% filter(tool %in% c("AMRFinderPlus_aa", "ABRicate_MEGARES_nt") & !variable %in% "fARGene_nt" & value > 0) %>% 
             group_by(tool, new_level) %>% 
  summarise(mn_recall = mean(recall, na.rm = T), med_recall = median(recall, na.rm = T),
            mn_fnr = mean(fnr, na.rm = T), med_fnr = median(fnr, na.rm = T),
            iqr_recall = IQR(recall, na.rm = T), iqr_fnr = IQR(fnr, na.rm = T),
            sum_row = sum(row_sum - value, na.rm = T)) %>%
  ungroup() %>%
  group_by(tool) %>% 
  summarise(MN_recall = mean(mn_recall, na.rm = T), MED_recall = median(med_recall, na.rm = T),
            MN_fnr = mean(mn_fnr, na.rm = T), MED_fnr = median(med_fnr, na.rm = T),
            iqr_med_recall = IQR(med_recall, na.rm = T), iqr_med_fnr = IQR(med_fnr, na.rm = T),
            SUM_row = sum(sum_row, na.rm = T))

# summary 1 - calculation of recall and fnr per tool and gene class
hm_recall_fnr_1 <- PM_new %>% filter(tool != variable) %>% ungroup() %>% filter( !is.na(col_sum)) %>%
  mutate(recall = value / col_sum, fnr = (row_sum - value) / row_sum ) %>% 
  mutate(recall = ifelse(is.na(row_sum), NA, recall), fnr = ifelse(is.na(row_sum), NA, fnr))

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

 

hm_recall_fnr_2 <- hm_recall_fnr_2 %>% mutate(tool = as.vector(alt_name_tools_rev[as.character(tool)]))
hm_recall_fnr_2 <- hm_recall_fnr_2 %>% mutate(tool = factor(tool, levels = tools_levels))

general_size <- 10

plot_recall <- ggplot(hm_recall_fnr_2, aes(x = tool, y = med_recall, fill = tool)) +
  geom_boxplot() +
  scale_fill_manual(values = pal_10) +
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
  scale_fill_manual(values = pal_10) +
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


df$facet_var <- gsub(" ", "\n", df$facet_var)

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

pal15_rep <- c(rep("#a6cee3",2), rep("#1f78b4",3), rep("#b2df8a",2), 
  rep("#33a02c", 1), rep("#fb9a99",2), pal_10[6:10])

plot_count_genes_tool <- unigenes  %>% mutate(tool = factor(tool, levels = tools_levels)) %>%
  ggplot(aes( x = tool)) +
  geom_bar(aes(fill = tool), color = "black") +
  scale_fill_manual(values = pal15_rep) +
  theme_minimal() +
  #scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
  #              labels = scales::trans_format("log10", scales::math_format(10^.x)))  + 
  ylab("Number of detected ARGs") +
  xlab("Tool") +
  labs(fill = "") +
  scale_y_continuous(breaks = c(25000, 50000,75000,100000,125000), labels = c("25,000", "50,000", "75,000","100,000","125,000")) +
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




unigenes_per_tool %>% filter(tool %in% tool_2) %>% 
  group_by(query) %>% 
  mutate(n = n_distinct(tool))  %>% 
  mutate(n = n > 1) %>% 
ggplot( aes(id, colour = tool, linetype = n)) +
  geom_freqpoly(binwidth = 1) +
  facet_wrap(. ~ tool, scales = "free_y")


table(lst$amrfinder.norm.prot$ARG.class, lst$amrfinder.norm.prot$Method)

unigenes_per_tool %>% filter(tool %in% tool_2 & id < 60) %>% 
  group_by(query) %>% 
  mutate(n = n_distinct(tool))  %>% 
  mutate(n = n > 1) %>% 
  ggplot( aes(id, colour = tool, linetype = n)) +
  geom_freqpoly(binwidth = 1) +
  facet_wrap(. ~ tool, scales = "free_y")


unigenes_per_tool %>% filter(tool %in% tool_2 & id < 70) %>% 
  group_by(query) %>% 
  mutate(n = n_distinct(tool))  %>% 
  mutate(n = n > 1) %>% 
  ggplot( aes(new_level, fill = tool)) +
  geom_bar() +
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

t70 %>% filter(p_bool< 0.9)
data.frame(hm_recall_fnr_2 %>% filter(tool %in% "DeepARG (nt)") %>% 
  arrange(med_recall) %>% filter(new_level %in% (t70 %>% filter(p_bool> 0.9, !n, tool %in% "DeepARG (nt)") %>% select(new_level) %>% pull())))

data.frame(hm_recall_fnr_2 %>% filter(tool %in% "RGI (DIAMOND - nt)") %>% 
             arrange(med_recall) %>% filter(new_level %in% (t70 %>% filter(p_bool> 0.9, !n, tool %in% "RGI (DIAMOND - nt)") %>% select(new_level) %>% pull())))

data.frame(hm_recall_fnr_2 %>% filter(tool %in% "RGI (DIAMOND - nt)") %>% 
             arrange(med_recall) %>% filter(new_level %in% (t50 %>% filter(p_bool> 0.9, !n, tool %in% "RGI (DIAMOND - nt)") %>% select(new_level) %>% pull())))

sum(lst$rgi.diamond$id < 70)





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

sumcore <- core %>% filter(cut %in% 0.5 & cnt >= 400, !habitat %in% c( "amplicon", "isolate",  "built-environment" ),
                           tool %in% tool_2) %>% ungroup() %>% 
  group_by(new_level, tool, habitat) %>% summarise(unigenes = n_distinct(X))  %>% 
  mutate(tool = factor(tool, levels = tool_2))

sumpan <-pan %>% ungroup() %>% group_by(tool, habitat, aggregation, gene_class) %>% 
  summarise(md = median(unigenes), mn = mean(unigenes)) %>%
  filter(aggregation == "new_level", !habitat %in% c( "amplicon", "isolate",  "built-environment" ),
         tool %in% tool_2) %>% mutate(tool = factor(tool, levels = tool_2))

sumpan %>% filter(habitat %in% c("human gut")) %>% 
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


# tools selected to plot
tool_selected <- c("rgi.diamond",  "deeparg", "fargene", "amrfinder.prot", 
                   "abricate.argannot", "abricate.card", "abricate.megares", "abricate.ncbi",
                   "abricate.resfinder", "resfinder")


# number of samples per habitat
abundance_parent <- abundance_parent %>% group_by(habitat) %>% mutate(n1 = n_distinct(sample))
abundance_parent <- abundance_parent %>% group_by(habitat2) %>% mutate(n2 = n_distinct(sample))

diversity_parent <- diversity_parent %>% group_by(habitat) %>% mutate(n1 = n_distinct(sample))
diversity_parent <- diversity_parent %>% group_by(habitat2) %>% mutate(n2 = n_distinct(sample))



# Abundance boxplot  habitat 2

fig1 <- abundance_parent  %>% 
  group_by(habitat, habitat2, tool, sample) %>% summarise(total = sum(normed10m)+1) %>%
  filter(!habitat %in% not_env, tool %in% tool_selected, !sample %in% extreme_samples) %>% 
  ggplot(aes( x = habitat2)) +
  geom_boxplot(aes(y = total, fill = tool), outlier.shape = NA) +
  scale_fill_manual(values = pal) +
  theme_minimal() +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))  + 
  ylab("Normalized abundance") +
  xlab("Source") +
  labs(fill = "") +
  theme(
    legend.position = "bottom",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(angle = 90))
  #scale_y_continuous(labels = scales::log10_trans())
fig1


pdf("fig1_abundance.pdf", width = 14, height = 10)
fig1
dev.off()

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
