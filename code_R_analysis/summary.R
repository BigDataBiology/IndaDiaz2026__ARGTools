library(dplyr)
library(ggplot2)
library(ggVennDiagram)
library(gridExtra)
library(tidyverse)
library(RColorBrewer)
library(ggpattern)
library(grid)

# 130mm 1 col 185 2 col
# 180mm 1 col 210 2 col
# 220mm 1 col 225 2 col

# d is divergent / q is qualitative
pal_8_q <- brewer.pal(8, "Dark2")
pal_10_q <- brewer.pal(10, "Paired")
pal_11_q <- brewer.pal(11, "Paired")
pal_6_q <- brewer.pal(6, "Paired")
pal_10_d <- brewer.pal(10, "BrBG")
pal_4_d <- brewer.pal(4, "BrBG")
pal_11_d <- brewer.pal(11, "BrBG")
pal_9_d <- brewer.pal(9, "BrBG")

# pallet for 15 tool + format repeated by tool
pal15_rep <- c(rep(pal_10_d[1], 2), rep(pal_10_d[2], 3), rep(pal_10_d[3], 2),
               pal_10_d[4], rep(pal_10_d[5], 2), pal_10_d[6:10])

# size for font in plots 
general_size <- 8

# data

setwd("~/Documents/GitHub/arg_compare/")
df2 <- readRDS(file = "code_R_analysis/output_abundance_diversity_resistome/conversion_ARO_parent_new_level.rds")
lst <- readRDS("code_R_analysis/output_abundance_diversity_resistome/results_tools.rds")
metadata <- read.delim("data/metadata_GMGC10.sample.meta.tsv")
abundance <- readRDS("code_R_analysis/output_abundance_diversity_resistome/abundance_diversity.rds")

# args_abundances <- read.delim("data/abundances/args_abundances.tsv") 
# args_abundances <- args_abundances %>% mutate(habitat = metadata$habitat[match(sample, metadata$sample_id)])

# change the name of the rgi runs 

lst$rgi.blast <- lst$rgi.blast %>% mutate(tool =  "RGI (BLAST - nt)")
lst$rgi.diamond <- lst$rgi.diamond %>% mutate(tool =  "RGI (DIAMOND - nt)")
lst$rgi.diamond.prot <- lst$rgi.diamond.prot %>% mutate(tool =  "RGI (DIAMOND - aa)")


# HABITATS
EN <- c("human gut", "human oral",  "human skin", "human nose", "human vagina", 
        "dog gut", "cat gut", "mouse gut", "pig gut", "wastewater", "marine", 
        "freshwater", "soil" , "amplicon", "isolate",  "built-environment" )

# SOURCE FOR EACH HABITAT
SO <- c(rep("Humans", 5), rep("Mammals", 4),  
        "Wastewater", "Marine", "Freshwater", 
        "Soil", rep("Other", 2), "Built-environment")

names(SO) <- EN

tools_levels <- c("DeepARG (nt)", "fARGene (nt)",
                  "RGI (DIAMOND - nt)", "ABRicate (CARD - nt)",
                  "ResFinder (nt)", "ABRicate (ResFinder - nt)", 
                  "AMRFinderPlus (aa)", "ABRicate (NCBI - nt)",
                  "ABRicate (ARG-ANNOT - nt)", "ABRicate (MEGARES - nt)")

tool_2 <- c("DeepARG (nt)", "RGI (DIAMOND - nt)", "fARGene (nt)", "AMRFinderPlus (aa)", "ResFinder (nt)",
            "ABRicate (ARG-ANNOT - nt)", "ABRicate (CARD - nt)", "ABRicate (MEGARES - nt)", "ABRicate (NCBI - nt)",
            "ABRicate (ResFinder - nt)")

tool_label <- c(expression({DeepARG~""^"a"}), expression({fARGene~""^"a"}), 
                expression({RGI~""^"a,b"}),  expression({CARD~""^"a,d"}), 
                expression(ResFinder~{""^"a"}), expression({ResFinder~""^"a,d"}),
                expression({AMRFinderPlus~""^"c"}), expression({NCBI~""^"a,d"}),
                expression({ARG-ANNOT~""^"a,d"}), 
                expression({MEGARES~""^"a,d"}))


abundance <- abundance %>% mutate( tool = 
                     ifelse(tool == "RGI (BLAST nt)", "RGI (BLAST - nt)",
                     ifelse(tool == "RGI (BLAST aa)", "RGI (BLAST - aa)",
                     ifelse(tool == "RGI (DIAMOND aa)", "RGI (DIAMOND - aa)",
                     ifelse(tool == "RGI (DIAMOND nt)", "RGI (DIAMOND - nt)", tool)))))

# environments that we are not interested in
not_env <- c("amplicon", "isolate")


# changing habitats and tools to factor

abundance <- abundance %>% mutate(habitat = factor(habitat, levels = EN),
                                  habitat2 = factor(SO[habitat], levels = c("Humans","Mammals","Wastewater","Built-environment", "Soil","Marine","Freshwater","Other")),
                                  tool = factor(tool, levels = tools_levels))

abundance <- abundance %>% mutate(location = ifelse(habitat2 %in% c("Humans","Mammals","Wastewater","Built-environment"), "Human-related","External"))

abundance <- abundance %>% mutate(location = factor(location, levels = c("Human-related","External")))


# keep original abundance data frame in abundance0
abundance0 <- abundance

# filter for analysis of environments and tools wanted 

abundance <- abundance %>% filter(tool %in% tool_2 & !habitat %in% not_env)
abundance <- abundance %>% filter(aggregation %in% "new_level")



# SUMMARIES
# SUMMARIES
# SUMMARIES

# unigenes captured with the tools
unigenes <- do.call(rbind, lapply(lst, function(x) x[,c("query","tool", "ARO", "parent", "parent_description", "new_level", "id")])) 
unigenes <- unigenes %>% mutate(tool = factor(tool, levels = tools_levels))

# total number of unique unigenes captured with the tools
unigenes %>% select(query) %>% distinct() %>% summarise(n = n())

# unigenes per tool
unigenes  %>% group_by(tool) %>% summarise(n = n_distinct(query)) %>% ungroup() %>% arrange(n)

plot_count_genes_tool <- unigenes  %>% filter(tool %in% tool_2) %>% 
  ggplot(aes( x = tool)) +
  geom_bar(aes(fill = tool), color = "black") +
  scale_fill_manual(values = pal_10_q) +
  scale_x_discrete( labels = tool_label) +
  theme_minimal() +
  ylab("Unigenes") +
  xlab("") +
  labs(fill = "") +
  scale_y_continuous(breaks = c(0, 25000, 50000,75000,100000,125000), labels = c("0","25,000", "50,000", "75,000","100,000","125,000")) +
  theme(
    legend.position = "none",
    legend.text = element_text(size = general_size),
    panel.border = element_rect(fill = "transparent", color = "black", linewidth = 1),
    panel.grid.major.x = element_blank(),
    plot.margin = margin(0, 0, 0, 0, unit = "pt"),
    legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
    legend.margin = margin(0, 0, 0, 0, unit = "pt"),
    panel.spacing = unit(0, "pt"),
    panel.grid.minor.x = element_blank(),
    axis.title = element_text(size = general_size + 1, face = "bold"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = general_size),
    axis.text.y = element_text(size = general_size))

plot_count_genes_tool

dev.off()
ggsave("~/Documents/plots_project2/count_genes_per_tool_2.svg", plot_count_genes_tool, width = 130, height = 80, unit = "mm")
dev.off()  


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

# classes making up (together) at least 50%  of genes in each tool
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
  scale_fill_manual(values = pal_11_q) +
  scale_x_discrete(labels = tool_label) +
  theme_minimal() +
  ylab("Proportion") +
  xlab("") +
  labs(fill = "") +
  theme(
    legend.position = "right",
    legend.text = element_text(size = general_size),
    panel.border = element_rect(fill = "transparent", color = "black", linewidth = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    plot.margin = margin(0, 0, 0, 0, unit = "pt"),
    legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
    legend.margin = margin(0, 0, 0, 0, unit = "pt"),
    panel.spacing = unit(0, "pt"),
    axis.title = element_text(size = general_size + 1, face = "bold"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = general_size),
    axis.text.y = element_text(size = general_size))

top10class_plot

dev.off()
ggsave("~/Documents/plots_project2/top10classes.svg", top10class_plot, width = 130, height = 100, unit = "mm")
dev.off()  


##############################################################

## parragraph on abundance fold differences 

ab_tool_env <- abundance  %>% 
  group_by(habitat2, tool, sample, aggregation) %>% summarise(total = sum(normed10m)+1) %>%
  filter(!habitat2 %in% "Other") %>% ungroup() %>% group_by(tool, habitat2) %>% summarise(md = median(total), mn = mean(total)) 

ab_tool_env %>% filter(habitat2 %in% "Humans") %>% ungroup() %>% arrange(desc(md))
ab_tool_env %>% filter(habitat2 %in% "Mammals") %>% ungroup() %>% arrange(desc(md))
ab_tool_env %>% filter(habitat2 %in% "Wastewater") %>% ungroup() %>% arrange(desc(md))
ab_tool_env %>% filter(habitat2 %in% "Soil") %>% ungroup() %>% arrange(desc(md))
ab_tool_env %>% filter(habitat2 %in% "Marine") %>% ungroup() %>% arrange(desc(md))
ab_tool_env %>% filter(habitat2 %in% "Freshwater") %>% ungroup() %>% arrange(desc(md))
ab_tool_env %>% filter(habitat2 %in% "Built-environment") %>% ungroup() %>% arrange(desc(md))

##############################################################

abundance_plot <- abundance  %>% 
  group_by(habitat2, tool, sample, aggregation) %>% summarise(total = sum(normed10m)+1) %>%
  filter(!habitat2 %in% "Other") %>% 
  ggplot(aes( x = habitat2)) +
  geom_boxplot(aes(y = total, fill = tool), outlier.shape = NA) +
  scale_fill_manual(values = pal_10_q, labels = tool_label) +
  theme_minimal() +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)),
                limits = c(1, 20000)) + 
  ylab("Abundance") +
  xlab("") +
  labs(fill = "") +
  theme(
    legend.position = "right",
    legend.text = element_text(size = general_size ),
    panel.border = element_rect(fill = "transparent", color = "black", linewidth = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    plot.margin = margin(0, 0, 0, 0, unit = "pt"),
    legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
    legend.margin = margin(0, 0, 0, 0, unit = "pt"),
    panel.spacing = unit(0, "pt"),
    axis.title = element_text(size = general_size + 1, face = "bold"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = general_size),
    axis.text.y = element_text(size = general_size))

abundance_plot

dev.off()
ggsave("~/Documents/plots_project2/abundance.svg", abundance_plot, width = 130, height = 80, unit = "mm")
dev.off()  


abundance_plot_inv <- abundance  %>% 
  group_by(habitat2, tool, sample, aggregation) %>% summarise(total = sum(normed10m)+1, location = location[1]) %>%
  filter(!habitat2 %in% "Other") %>% 
  ggplot(aes( x = tool)) +
  scale_x_discrete(labels = tool_label) +
  geom_boxplot(aes(y = total, fill = habitat2), outlier.shape = NA) +
  facet_grid( . ~ location) +
  scale_fill_manual(values = pal_10_q) +
  theme_minimal() +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)),
                limits = c(1, 20000)) + 
  ylab("Abundance") +
  xlab("") +
  labs(fill = "") +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = general_size ),
    panel.border = element_rect(fill = "transparent", color = "black", linewidth = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    plot.margin = margin(0, 0, 0, 0, unit = "pt"),
    legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
    legend.margin = margin(0, 0, 0, 0, unit = "pt"),
    panel.spacing = unit(0, "pt"),
    strip.text = element_text(size = general_size, face = "bold"),
    axis.title = element_text(size = general_size + 1, face = "bold"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = general_size),
    axis.text.y = element_text(size = general_size))

abundance_plot_inv

dev.off()
ggsave("~/Documents/plots_project2/abundance_inverted.svg", abundance_plot_inv, width = 130, height = 80, unit = "mm")
dev.off()  



abundance_class <- abundance %>% 
  ungroup() %>% filter(tool %in% tool_2, !habitat %in% not_env, aggregation %in% "new_level") %>% 
  mutate(tool = as.character(tool), gene = as.character(gene), sample = as.character(sample)) %>%
  complete(sample, gene, tool, 
           fill = list(normed10m = 0, scaled = 0, raw = 0, raw_unique = 0, unigenes = 0)) %>%
  mutate(tool = factor(tool, levels = tools_levels))

abundance_class <- abundance_class %>% mutate(habitat = metadata$habitat[match(sample, metadata$sample_id)])
abundance_class <- abundance_class %>% mutate(habitat2 = metadata$habitat2[match(sample, metadata$sample_id)])
abundance_class <- abundance_class %>% mutate(location = metadata$location[match(sample, metadata$sample_id)])
abundance_class <- abundance_class %>% mutate(aggregation = abundance$aggregation[match(gene, abundance$gene)])


top20 <- abundance_class %>% 
  filter(habitat %in% c("human gut"), tool %in% tool_2) %>% 
  ungroup() %>% 
  group_by(tool, gene) %>%
  summarise(mn = median(normed10m), s = sum(normed10m)) %>% 
  arrange(desc(mn)) %>% 
  ungroup() %>%
  group_by(tool) %>% 
  arrange(desc(mn), desc(s)) %>% 
#  filter(mn > 0) %>% 
  slice_head(n = 15) %>% 
  ungroup() %>%
  select(gene) %>% 
  distinct() %>% 
  pull()

top20 <- abundance_class %>% 
  filter(habitat %in% c("human gut"), tool %in% tool_2) %>% 
  ungroup() %>% 
  group_by(tool, gene) %>%
  summarise(mn = median(normed10m), s = sum(normed10m)) %>% 
  arrange(desc(mn)) %>% 
  ungroup() %>%
  group_by(tool) %>% 
  arrange(desc(mn), desc(s)) %>% 
  #filter(mn > 0) %>% 
  slice_head(n = 5) %>% 
  ungroup() %>%
  select(gene) %>% 
  distinct() %>% 
  pull()

top20 <- unique(c(top20, "Cell wall charge", "ERM", "Class C"))
top20 <- factor(c(top20, "Other"), levels = c(top20, "Other"))


ab_class_tool <- abundance_class  %>% 
  filter(habitat %in% "human gut") %>% 
  group_by(tool, gene, sample) %>% 
  mutate(g2 = ifelse(gene %in% top20, gene, "Other")) %>%
  ungroup() %>% 
  group_by(tool, g2) %>% 
  summarise(total = sum(normed10m), mn = mean(normed10m) , md = median(normed10m) ) %>% 
  mutate(p_total = total / sum(total)) 
  
ab_class_tool %>% filter(g2 %in% "rpoB") %>% arrange(desc(md), desc(total))
ab_class_tool %>% filter(g2 %in% "TET - RPG") %>% arrange(desc(md), desc(total))
ab_class_tool %>% filter(g2 %in% "GPA") %>% arrange(desc(md), desc(total))
ab_class_tool %>% filter(g2 %in% "Efflux p.") %>% arrange(desc(md), desc(total))
ab_class_tool %>% filter(g2 %in% "Cell wall charge") %>% arrange(desc(md), desc(total))

ab_class_tool %>% filter(tool %in% "DeepARG (nt)") %>% arrange(desc(md), desc(total)) %>% print(n = 40)
ab_class_tool %>% filter(tool %in% "fARGene (nt)") %>% arrange(desc(md), desc(total)) %>% print(n = 40)
ab_class_tool %>% filter(tool %in% "ResFinder (nt)") %>% arrange(desc(md), desc(total))
ab_class_tool %>% filter(tool %in% "ResFinder (nt)") %>% arrange(desc(md), desc(total))
ab_class_tool %>% filter(tool %in% "ResFinder (nt)") %>% arrange(desc(md), desc(total))
ab_class_tool %>% filter(tool %in% "ResFinder (nt)") %>% arrange(desc(md), desc(total))
ab_class_tool %>% filter(tool %in% "ResFinder (nt)") %>% arrange(desc(md), desc(total))


human_abundance_class_plot <- abundance_class  %>% 
  filter(habitat %in% "human gut") %>% 
  group_by(habitat, tool, gene, sample) %>% summarise(total = normed10m + 1) %>%
  ungroup() %>%
  mutate(g2 = ifelse(gene %in% top20, gene, "Other")) %>% 
  mutate(g2 = factor(g2, levels = top20)) %>% 
  filter(g2 %in% top20) %>%
  ggplot(aes( x = g2)) +
  geom_boxplot(aes(y = total, fill = tool), position = position_dodge2(preserve = "single"), outlier.shape = NA) +
  scale_fill_manual(values = pal_10_q, labels = tool_label) +
  theme_minimal() +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)),
                limits = c(1, 20000)) + 
  ylab("Abundance") +
  xlab("") +
  labs(fill = "") +
  theme(
    legend.position = "right",
    legend.text = element_text(size = general_size ),
    panel.border = element_rect(fill = "transparent", color = "black", linewidth = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    plot.margin = margin(0, 0, 0, 0, unit = "pt"),
    legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
    legend.margin = margin(0, 0, 0, 0, unit = "pt"),
    panel.spacing = unit(0, "pt"),
    axis.title = element_text(size = general_size + 1, face = "bold"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = general_size),
    axis.text.y = element_text(size = general_size))


human_abundance_class_plot

dev.off()
ggsave("~/Documents/plots_project2/abundance_class_humangut.svg", human_abundance_class_plot, width = 130, height = 80, unit = "mm")
dev.off()  







##############################################################
##############################################################



div_tool_env <- abundance_class  %>% 
  filter(habitat %in% "human gut") %>% 
  group_by(tool, gene, sample) %>% 
  mutate(g2 = ifelse(gene %in% top20, gene, "Other")) %>%
  ungroup() %>% 
  group_by(tool, g2) %>% 
  summarise(total_div = sum(unigenes), mn_div = mean(unigenes) , md_div = median(unigenes) ) %>% 
  mutate(p_total_div = total_div / sum(total_div)) 

div_tool_env  %>% filter(tool %in% c("AMRFinderPlus (aa)", "ABRicate (NCBI - nt)")) %>% ungroup()
div_tool_env  %>% filter(tool %in% c("ResFinder (nt)", "ABRicate (ResFinder - nt)")) %>% ungroup() 
div_tool_env  %>% filter(tool %in% c("RGI (DIAMOND - nt)", "ABRicate (CARD - nt)")) %>% ungroup() 



div_tool_env %>% select(tool, g2, mn_div, md_div) %>% bind_cols(ab_class_tool)


human_diversity_class_plot <- abundance_class  %>% 
  filter(habitat %in% "human gut") %>% 
  group_by(habitat, tool, gene, sample) %>% summarise(total = unigenes + 1) %>%
  ungroup() %>%
  mutate(g2 = ifelse(gene %in% top20, gene, "Other")) %>% 
  mutate(g2 = factor(g2, levels = top20)) %>% 
  filter(g2 %in% top20) %>%
  ggplot(aes( x = g2)) +
  geom_boxplot(aes(y = total, fill = tool), position = position_dodge2(preserve = "single"), outlier.shape = NA) +
  scale_fill_manual(values = pal_10_q, labels = tool_label) +
  theme_minimal() +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)),
                limits = c(1, 2000)) + 
  ylab("Diversity") +
  xlab("") +
  labs(fill = "") +
  theme(
    legend.position = "right",
    legend.text = element_text(size = general_size ),
    panel.border = element_rect(fill = "transparent", color = "black", linewidth = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    plot.margin = margin(0, 0, 0, 0, unit = "pt"),
    legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
    legend.margin = margin(0, 0, 0, 0, unit = "pt"),
    panel.spacing = unit(0, "pt"),
    axis.title = element_text(size = general_size + 1, face = "bold"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = general_size),
    axis.text.y = element_text(size = general_size))

human_diversity_class_plot

diversity_plot <- abundance  %>% 
  group_by(habitat2, tool, sample, aggregation) %>% summarise(total = sum(unigenes) + 1) %>%
  filter(tool %in% tool_2, aggregation %in% "new_level") %>% 
  filter(!habitat2 %in% "Other") %>% 
  ggplot(aes( x = habitat2)) +
  geom_boxplot(aes(y = total, fill = tool), outlier.shape = NA) +
  scale_fill_manual(values = pal_10_q, labels = tool_label) +
  theme_minimal() +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)),
                limits = c(1, 6000)) + 
  ylab("Diversity") +
  xlab("") +
  labs(fill = "") +
  theme(
    legend.position = "right",
    legend.text = element_text(size = general_size),
    panel.border = element_rect(fill = "transparent", color = "black", linewidth = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title = element_text(size = general_size + 1, face = "bold"),
    plot.margin = margin(0, 0, 0, 0, unit = "pt"),
    legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
    legend.margin = margin(0, 0, 0, 0, unit = "pt"),
    panel.spacing = unit(0, "pt"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = general_size),
    axis.text.y = element_text(size = general_size)) #+


diversity_plot

dev.off()
ggsave("~/Documents/plots_project2/diversity.svg", diversity_plot, width = 130, height = 80, unit = "mm")
dev.off()  



diversity_plot_inv <- abundance  %>% 
  group_by(habitat2, tool, sample, aggregation) %>% summarise(total = sum(unigenes)+1, location = location[1]) %>%
  filter(!habitat2 %in% "Other") %>%
  ggplot(aes( x = tool)) +
  facet_grid( . ~ location) +
  scale_x_discrete(labels = tool_label) +
  geom_boxplot(aes(y = total, fill = habitat2), outlier.shape = NA) +
  scale_fill_manual(values = pal_10_q) +
  theme_minimal() +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)),
                limits = c(1, 5000)) + 
  ylab("Diversity") +
  xlab("") +
  labs(fill = "") +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = general_size ),
    panel.border = element_rect(fill = "transparent", color = "black", linewidth = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    plot.margin = margin(0, 0, 0, 0, unit = "pt"),
    legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
    legend.margin = margin(0, 0, 0, 0, unit = "pt"),
    panel.spacing = unit(0, "pt"),
    strip.text = element_text(size = general_size, face = "bold"),
    axis.title = element_text(size = general_size + 1, face = "bold"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = general_size),
    axis.text.y = element_text(size = general_size))


diversity_plot_inv

dev.off()  
ggsave("~/Documents/plots_project2/diversity_inverted.svg", diversity_plot_inv, width = 130, height = 80, unit = "mm")

abundance_plot_inv

### Correlations
abundance
cor_abu_div <- abundance_class %>% mutate(habitat2 = abundance$habitat2[match(habitat, abundance$habitat)]) %>%
  group_by(habitat, tool) %>% filter(normed10m != 0 & unigenes != 0) %>% summarise(corr = cor(normed10m, unigenes, method = "pearson"))
data.frame(cor_abu_div)


##############################################################
##############################################################

abu_tool_habitat <- abundance  %>% 
  group_by(habitat, tool, sample, aggregation, gene) %>% 
  summarise(abundance = sum(normed10m) + 1,
            diversity = sum(unigenes) + 1) %>%
  filter(!habitat %in% not_env, tool %in% tool_2) 

# tool_2[c(1:3,5)] = "DeepARG (nt)" "RGI (DIAMOND - nt)" "fARGene (nt)" "ResFinder (nt)" 
human.genes <- abu_tool_habitat  %>% 
  filter(habitat %in% c("human gut"), aggregation %in% "new_level", tool %in% tool_2[c(1:3,6)]) %>%
  ungroup() %>% group_by( tool, gene) %>% summarise(n = median(abundance)) %>%
  ungroup() %>% arrange( tool, desc(n)) %>% group_by(tool)  %>% slice_head(n = 5) %>% 
  ungroup() %>% 
  select(gene) %>% 
  distinct() %>% 
  pull()

abundance_plot_habitat <- 
  abu_tool_habitat  %>% 
  filter(habitat %in% c("human gut"), aggregation %in% "new_level", as.character(tool) %in% tool_2[c(1:3,5)], gene %in% human.genes) %>%
  mutate(gene = factor(gene, levels = factor_new_level)) %>% 
  ggplot(aes( x = gene)) +
  geom_boxplot(aes(y = abundance, fill = tool), outlier.shape = NA, position = position_dodge(preserve = "single")) +
  scale_fill_manual(values = pal_6_q) +
  theme_minimal() +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)),
                limits = c(1, 20000)) + 
  ylab("Abundance") +
  xlab("") +
  labs(fill = "") +
  facet_grid(habitat ~ gene, scales = "free_x") +
  theme(
    legend.position = "right",
    legend.text = element_text(size = general_size),
    panel.border = element_rect(fill = "transparent", color = "black", linewidth = 1),
    panel.grid.major.x = element_blank(),
    plot.margin = margin(0, 0, 0, 0, unit = "pt"),
    legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
    legend.margin = margin(0, 0, 0, 0, unit = "pt"),
    panel.spacing = unit(0, "pt"),
    panel.grid.minor.x = element_blank(),
    strip.text = element_blank(),
    axis.title = element_text(size = general_size + 1, face = "bold"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = general_size),
    axis.text.y = element_text(size = general_size))

abundance_plot_habitat

dev.off()
ggsave("~/Documents/plots_project/abundance_human_gut.svg", abundance_plot_habitat, width = 130, height = 80, unit = "mm")
dev.off()  


diversity_plot_habitat <- 
  abu_tool_habitat  %>% 
  filter(habitat %in% c("human gut"), aggregation %in% "new_level", tool %in% tool_2[c(1:3,5)], gene %in% human.genes) %>%
  mutate(gene = factor(gene, levels = factor_new_level)) %>% 
  ggplot(aes( x = gene)) +
  geom_boxplot(aes(y = diversity, fill = tool), outlier.shape = NA, position = position_dodge(preserve = "single")) +
  scale_fill_manual(values = pal_6_q) +
  theme_minimal() +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)),
                limits = c(1, 1300)) + 
  ylab("Diversity") +
  xlab("") +
  labs(fill = "") +
  facet_grid(habitat ~ gene, scales = "free_x") +
  theme(
    legend.position = "right",
    legend.text = element_text(size = general_size),
    panel.border = element_rect(fill = "transparent", color = "black", linewidth = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    plot.margin = margin(0, 0, 0, 0, unit = "pt"),
    legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
    legend.margin = margin(0, 0, 0, 0, unit = "pt"),
    panel.spacing = unit(0, "pt"),
    strip.text = element_blank(),
    axis.title = element_text(size = general_size + 1, face = "bold"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = general_size),
    axis.text.y = element_text(size = general_size))

diversity_plot_habitat

dev.off()
ggsave("~/Documents/plots_project/diversity_human_gut.svg", diversity_plot_habitat, width = 130, height = 80, unit = "mm")
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


# double new level
double_level <- do.call(rbind, lapply(lst, function(x) x[,c("query","tool","new_level")])) %>%
  ungroup() %>% 
  group_by(query) %>% 
  summarise(n = n_distinct(new_level)) %>% 
  filter(n>1) %>% select(query) %>% 
  pull()


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


###
# overlap all genes between tools and classes
# create lists

sets0 <- tools_per_unigene %>%
  group_by(tool) %>%
  summarise(query = list(query), .groups = "drop")

sets1 <- tools_per_unigene %>%
  group_by(new_level, tool) %>%
  summarise(query = list(query), .groups = "drop")

# Pairwise combinations 
pairwise <- sets1 %>%
  group_by(new_level) %>%
  summarise(pairs = list(expand_grid(tool_ref = tool, tool_comp = tool)), .groups = "drop") %>%
  unnest(pairs)

JI_class_other <- pairwise %>%
  left_join(sets1, by = c("new_level", "tool_ref" = "tool")) %>%
  rename(qc_ref = query) %>%
  left_join(sets1, by = c("new_level", "tool_comp" = "tool")) %>%
  rename(qc_comp = query) %>%
  left_join(sets0, by = c( "tool_ref" = "tool")) %>%
  rename(q_ref = query) %>% 
  left_join(sets0, by = c( "tool_comp" = "tool")) %>%
  rename(q_comp = query)


new_intersect <- function(qc_ref, q_ref, qc_comp, q_comp){
  A <-  unlist(qc_ref)
  B <- unlist(q_ref)
  C <- unlist(qc_comp)
  D <- unlist(q_comp)
  x <- intersect(A, setdiff(D, C))
  y <- setdiff(C, intersect(C, setdiff(B, A)))
  z <- union(x, y)
  d <- intersect(A, z)
  r <- ifelse(length(z) == 0, NA, length(d) / length(z))
  return(r)
}

new_difference <- function(qc_ref, q_ref, qc_comp, q_comp){
  A <-  unlist(qc_ref)
  B <- unlist(q_ref)
  C <- unlist(qc_comp)
  D <- unlist(q_comp)
  x <- intersect(A, setdiff(D, C))
  y <- setdiff(C, intersect(C, setdiff(B, A)))
  z <- union(x, y)
  d <- setdiff(A, z)
  r <- ifelse(length(A) == 0, NA, length(d) / length(A))
  return(r)
}


JI_class_other <- JI_class_other %>% rowwise() %>% mutate(recall = new_intersect(qc_ref, q_ref, qc_comp, q_comp))
JI_class_other <- JI_class_other %>% rowwise() %>% mutate(fnr = new_difference(qc_ref, q_ref, qc_comp, q_comp))
JI_class_other <- JI_class_other %>% 
  mutate(ref_n_class = length(qc_ref), comp_n_class = length(qc_comp), ref_n_all = length(q_ref), comp_n_all = length(q_comp))
JI_class_other_filter <- JI_class_other %>% 
  ungroup() %>% filter(tool_ref != tool_comp) %>% 
  select(-c(qc_ref, qc_comp, q_ref, q_comp))


JI_class_other_filter %>% filter(tool_ref %in% "RGI (DIAMOND - nt)")
JI_class_other_filter %>% filter(tool_ref %in% "fARGene (nt)")
JI_class_other_filter %>% filter(tool_ref %in% "fARGene (nt)", new_level %in% "TET - enzyme")

##

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
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = general_size),
        axis.text.y = element_text(size = general_size),
        strip.text = element_text(size = general_size, face = "bold"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title = element_text(size = general_size + 1, face = "bold"))


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
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = general_size),
        axis.text.y = element_text(size = general_size),
        strip.text = element_text(size = general_size, face = "bold"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title = element_text(size = general_size + 1, face = "bold"))


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
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = general_size),
        axis.text.y = element_text(size = general_size),
        strip.text = element_text(size = general_size, face = "bold"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title = element_text(size = general_size + 1, face = "bold"))


# summary 1 - calculation of recall and fnr per tool and gene class


plot_recall_3tools <- JI_class_other_filter %>% filter(tool_ref %in% c("RGI (DIAMOND - nt)", "DeepARG (nt)", "fARGene (nt)", "ResFinder (nt)"),
                                          tool_comp %in% tool_2,
                                          new_level %in% c(human.genes, "Class C", "Class D", "MPH", "APH", "QNR")) %>% 
  ggplot(aes(x = new_level, y = recall)) +
  geom_boxplot(aes(fill = tool_ref)) +
  facet_grid( tool_ref ~ new_level , scales = "free_x",  labeller = labeller(tool = tool_label)) +
  scale_fill_manual(values = pal_6_q) +
  theme_minimal() +
  ylab("Recall") +
  xlab("Class") +
  labs(fill = "") +
  theme(
    legend.position = "none",
    panel.border = element_rect(fill = "transparent", color = "black", linewidth = 1),
    #legend.text = element_text(size = general_size),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    plot.margin = margin(0, 0, 0, 0, unit = "pt"),
    legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
    legend.margin = margin(0, 0, 0, 0, unit = "pt"),
    panel.spacing = unit(2, "pt"),
    strip.text.x = element_blank(),
    strip.text.y = element_text(size = general_size , face = "bold"),
    axis.title = element_text(size = general_size + 1, face = "bold"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = general_size),
    axis.text.y = element_text(size = general_size))

plot_recall_3tools

dev.off()
ggsave("~/Documents/plots_project/recall_box_1.svg", plot_recall_3tools, width = 200, height = 120, unit = "mm")
dev.off()  




plot_fdr_3tools <- JI_class_other_filter %>% filter(tool_ref %in% c("RGI (DIAMOND - nt)", "DeepARG (nt)", "fARGene (nt)", "ResFinder (nt)"),
                                       tool_comp %in% tool_2,
                                       new_level %in% c(human.genes, "Class C", "Class D", "MPH", "APH", "QNR")) %>% 
  ggplot(aes(x = new_level, y = fnr)) +
  geom_boxplot(aes(fill = tool_ref)) +
  facet_grid( tool_ref ~ new_level , scales = "free_x") +
  scale_fill_manual(values = pal_6_q) +
  theme_minimal() +
  ylab("False discovery rate") +
  xlab("Class") +
  labs(fill = "") +
  theme(
    legend.position = "none",
    panel.border = element_rect(fill = "transparent", color = "black", linewidth = 1),
    #legend.text = element_text(size = general_size),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    plot.margin = margin(0, 0, 0, 0, unit = "pt"),
    legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
    legend.margin = margin(0, 0, 0, 0, unit = "pt"),
    panel.spacing = unit(2, "pt"),
    strip.text.x = element_blank(),
    strip.text.y = element_text(size = general_size , face = "bold"),
    axis.title = element_text(size = general_size + 1, face = "bold"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = general_size),
    axis.text.y = element_text(size = general_size))

plot_fdr_3tools

dev.off()
ggsave("~/Documents/plots_project/fnr_box_1.svg", plot_fdr_3tools, width = 200, height = 120, unit = "mm")
dev.off()  



data.frame(JI_class_other_filter %>% filter(tool_ref %in% "AMRFinderPlus (aa)", tool_comp %in% tool_2) %>% 
             group_by(tool_ref, new_level) %>% summarise(M_R = median(recall)))

data.frame(JI_class_other_filter %>% filter(tool_ref %in% "fARGene (nt)", tool_comp %in% tool_2) %>% 
             group_by(tool_ref, new_level) %>% summarise(M_R = median(recall)))

JI_class_other_filter %>% filter(tool_ref %in% "fARGene (nt)", tool_comp %in% tool_2, new_level %in% c("TET - enzyme"))
JI_class_other_filter %>% filter(tool_ref %in% "fARGene (nt)", tool_comp %in% tool_2, new_level %in% c("APH"))
JI_class_other_filter %>% filter(tool_ref %in% "fARGene (nt)", tool_comp %in% tool_2, new_level %in% c("AAC"))
JI_class_other_filter %>% filter(tool_ref %in% "RGI (DIAMOND - nt)", tool_comp %in% tool_2, new_level %in% c("TET - enzyme"))
JI_class_other_filter %>% filter(tool_ref %in% "ABRicate (MEGARES - nt)", tool_comp %in% tool_2, new_level %in% c("TET - enzyme"))

data.frame(JI_class_other_filter %>% filter(tool_ref %in% "AMRFinderPlus (aa)" & !tool_comp %in% "fARGene (nt)" ) %>% 
             group_by(tool_ref, new_level) %>% summarise(M_R = median(recall)))


# summary 2 -  medians of recall and fnr per tool and gene class, accross all other tools

hm_recall_fnr_2 <- JI_class_other_filter %>% filter(tool_ref %in% tool_2, tool_comp %in% tool_2) %>% 
  group_by(tool_ref, new_level) %>%
  summarise(mn_recall = mean(recall, na.rm = T), med_recall = median(recall, na.rm = T),
            mn_fnr = mean(fnr, na.rm = T), med_fnr = median(fnr, na.rm = T),
            iqr_recall = IQR(recall, na.rm = T), iqr_fnr = IQR(fnr, na.rm = T),
            ref_n_class = ref_n_class[1], ref_n_all = ref_n_all[1])

# summary 2 -  medians of the medians of recall and fnr per tool accross gene class after accross tool

 hm_recall_fnr_2 %>% filter(tool_ref %in% c("fARGene (nt)", "AMRFinderPlus (aa)"), !is.na(med_recall)) 
 hm_recall_fnr_2 %>% filter(tool_ref %in% c("fARGene (nt)"), !is.na(med_recall)) 
 hm_recall_fnr_2 %>% filter(tool_ref %in% c("AMRFinderPlus (aa)"), !is.na(med_recall)) 

hm_recall_fnr_3 <- hm_recall_fnr_2 %>%
  ungroup() %>%
  group_by(tool_ref) %>% 
  summarise(MN_recall = mean(mn_recall, na.rm = T), MED_recall = median(med_recall, na.rm = T),
            MN_fnr = mean(mn_fnr, na.rm = T), MED_fnr = median(med_fnr, na.rm = T),
            iqr_med_recall = IQR(med_recall, na.rm = T), iqr_med_fnr = IQR(med_fnr, na.rm = T),
            ref_n_class = ref_n_class[1], ref_n_all = ref_n_all[1])

# unique genes per tool
tools_per_unigene %>% group_by(query) %>% filter(tool == "AMRFinderPlus (aa)") %>% 
  ungroup() %>% group_by(n_tools) %>% summarise(n = n()) %>% mutate( p = n/sum(n)) %>% filter(n_tools == 1)

tools_per_unigene %>% group_by(query) %>% filter(tool == "ABRicate (ResFinder - nt)") %>% 
  ungroup() %>% group_by(n_tools) %>% summarise(n = n()) %>% mutate( p = n/sum(n)) %>% filter(n_tools == 1)

tools_per_unigene %>% group_by(query) %>% filter(tool == "ABRicate (CARD - nt)") %>% 
  ungroup() %>% group_by(n_tools) %>% summarise(n = n()) %>% mutate( p = n/sum(n)) %>% filter(n_tools == 1)

tools_per_unigene %>% group_by(query) %>% filter(tool == "ABRicate (ARG-ANNOT - nt)") %>% 
  ungroup() %>% group_by(n_tools) %>% summarise(n = n()) %>% mutate( p = n/sum(n)) %>% filter(n_tools == 1)

tools_per_unigene %>% group_by(query) %>% filter(tool == "ABRicate_NCBI_nt") %>% 
  ungroup() %>% group_by(n_tools) %>% summarise(n = n()) %>% mutate( p = n/sum(n)) %>% filter(n_tools == 1)

tools_per_unigene %>% group_by(query) %>% filter(tool == "ABRicate (MEGARES - nt)") %>% 
  ungroup() %>% group_by(n_tools) %>% summarise(n = n()) %>% mutate( p = n/sum(n)) %>% filter(n_tools == 1)

tools_per_unigene %>% group_by(query) %>% filter(tool == "RGI_DIAMOND_nt") %>% 
  ungroup() %>% group_by(n_tools) %>% summarise(n = n()) %>% mutate( p = n/sum(n)) %>% filter(n_tools == 1)

tools_per_unigene %>% group_by(query) %>% filter(tool == "DeepARG (nt)") %>% 
  ungroup() %>% group_by(n_tools, new_level) %>% summarise(n = n())  %>% filter(n_tools == 1) %>% arrange(desc(n))


plot_recall <- ggplot(hm_recall_fnr_2, aes(x = tool_ref, y = med_recall, fill = tool_ref)) +
  geom_boxplot() +
  scale_fill_manual(values = pal_10_q) +
  scale_x_discrete(labels = tool_label) +
  theme_minimal() +
  ylab("Recall") +
  xlab("") +
  labs(fill = "") +
  theme(
    legend.position = "none",
    legend.text = element_text(size = general_size),
    panel.border = element_rect(fill = "transparent", color = "black", linewidth = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    plot.margin = margin(0, 0, 0, 0, unit = "pt"),
    legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
    legend.margin = margin(0, 0, 0, 0, unit = "pt"),
    panel.spacing = unit(0, "pt"),
    strip.text = element_blank(),
    axis.title = element_text(size = general_size + 1, face = "bold"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = general_size),
    axis.text.y = element_text(size = general_size))

plot_recall

dev.off()
ggsave("~/Documents/plots_project/recall.svg", plot_recall, width = 130, height = 80, unit = "mm")
dev.off()


plot_fnr <- ggplot(hm_recall_fnr_2, aes(x = tool_ref, y = med_fnr, fill = tool_ref)) +
  geom_boxplot() +
  scale_fill_manual(values = pal_10_q) +
  scale_x_discrete(labels = tool_label) +
  theme_minimal() +
  ylab("False Negative Rate") +
  xlab("") +
  labs(fill = "") + 
  theme(
    legend.position = "none",
    legend.text = element_text(size = general_size),
    panel.border = element_rect(fill = "transparent", color = "black", linewidth = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    plot.margin = margin(0, 0, 0, 0, unit = "pt"),
    legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
    legend.margin = margin(0, 0, 0, 0, unit = "pt"),
    panel.spacing = unit(0, "pt"),
    strip.text = element_blank(),
    axis.title = element_text(size = general_size + 1, face = "bold"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = general_size),
    axis.text.y = element_text(size = general_size))

plot_fnr

dev.off()
ggsave("~/Documents/plots_project/fnr.svg", plot_fnr, width = 130, height = 80, unit = "mm")
dev.off()

#df$facet_var <- gsub(" ", "\n", df$facet_var)

hm_recall_fnr_2 <- hm_recall_fnr_2 %>% 
  filter(!is.na(new_level)) %>% 
  mutate(strip = gsub("\\(", "\n\\(", tool_ref)) %>%
  mutate(strip = factor(strip, levels = gsub("\\(", "\n\\(", levels(tool_ref))))

recall_heatmap <- ggplot(hm_recall_fnr_2, aes(x = med_recall, y = new_level, fill = 1-iqr_recall)) +
  geom_col() +
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
        axis.title = element_text(size = general_size + 1, face = "bold"))

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
        axis.title = element_text(size = general_size + 1, face = "bold"))

fnr_heatmap


dev.off()
ggsave("~/Documents/plots_project/recall_heatmap.svg", recall_heatmap  , width = 14, height = 12)
dev.off()

dev.off()
ggsave("~/Documents/plots_project/fnr_heatmap.svg", fnr_heatmap, width = 14, height = 12)
dev.off()



idplot <- unigenes %>% filter(tool %in% tool_2) %>% 
  group_by(query) %>% 
  mutate(n = n_distinct(tool))  %>% 
  mutate(n = n > 1) %>% 
ggplot( aes(id, colour = tool, linetype = n)) +
  geom_freqpoly(binwidth = 1, linewidth = 2) +
  facet_wrap(. ~ tool, scales = "free_y", labeller = labeller(tool = tool_label)) +
  scale_color_manual(values = rep(pal_10_d[9],10)) +
  theme_minimal() +
  labs(fill = "Jaccaard index") +
  xlab("Identity") +
  ylab("Unigenes") +
  theme(axis.title.y = element_blank(),
        legend.position = "none",
        panel.border = element_rect(fill = "transparent", color = "black", linewidth = 1),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = general_size),
        axis.text.y = element_text(size = general_size),
        strip.text = element_text(size = general_size, face = "bold"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title = element_text(size = general_size + 1, face = "bold"))


idplot

dev.off()
ggsave("~/Documents/plots_project/ids_all.svg", idplot, width = 130, height = 80, unit = "mm")


idplot2 <- unigenes %>% filter(tool %in% tool_2[1:4]) %>% 
  group_by(query) %>% 
  mutate(n = n_distinct(tool))  %>% 
  mutate(n = n > 1) %>% 
  ggplot( aes(id, colour = tool, linetype = n)) +
  geom_freqpoly(binwidth = 1, linewidth = 2) +
  facet_wrap(. ~ tool, scales = "free_y", labeller = labeller(tool = tool_label)) +
  scale_color_manual(values = rep(pal_10_d[9],10)) +
  theme_minimal() +
  labs(fill = "Jaccaard index") +
  xlab("Identity") +
  ylab("Unigenes") +
  theme(axis.title.y = element_blank(),
        legend.position = "none",
        panel.border = element_rect(fill = "transparent", color = "black", linewidth = 1),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = general_size),
        axis.text.y = element_text(size = general_size),
        strip.text = element_text(size = general_size, face = "bold"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title = element_text(size = general_size + 1, face = "bold"))


dev.off()
ggsave("~/Documents/plots_project/ids.svg", idplot2, width = 130, height = 80, unit = "mm")

#######





table(lst$amrfinder.norm.prot$ARG.class, lst$amrfinder.norm.prot$Method)

#unigenes %>% filter(tool %in% tool_2 & id < 80) %>% 
#  group_by(query) %>% 
#  mutate(n = n_distinct(tool))  %>% 
#  mutate(n = n > 1) %>% 
#  ggplot( aes(id, colour = tool, linetype = n)) +
#  geom_freqpoly(binwidth = 1) +
#  facet_wrap(. ~ tool, scales = "free_y")



t70 <- unigenes %>% 
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

t60 <- unigenes %>% 
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

t50 <- unigenes %>% 
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



t70_plot <- t70 %>% 
  ggplot(aes(x = new_level, y = N_tool_class_bool)) +
  geom_col(aes(fill=n)) + 
  scale_fill_manual(values = pal_6_q) +  
  xlab("") + 
  ylab("Unigenes") + 
  labs(fill = "In other tool") +
  facet_grid(tool ~  . , scales = "free") +
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.text = element_text(size = general_size),
        panel.border = element_rect(fill = "transparent", color = "black", linewidth = 1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        plot.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.margin = margin(0, 0, 0, 0, unit = "pt"),
        panel.spacing = unit(0, "pt"),
        strip.text = element_text(size = general_size, face = "bold"),
        axis.title = element_text(size = general_size + 1, face = "bold"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = general_size),
        axis.text.y = element_text(size = general_size))

ggsave("~/Documents/plots_project/t70.svg", t70_plot, width = 130, height = 120, unit = "mm")


t70_plot_proportion <- t70 %>% filter(n == FALSE) %>%  
  ggplot(aes(x = new_level, y = P)) +
  geom_col(aes(fill=n)) + 
  scale_fill_manual(values = pal_6_q) +  
  ylab("Proportion in gene class") + 
  xlab("") +
  labs(fill = "In other tool") +
  facet_grid(tool ~  . , scales = "free") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        panel.border = element_rect(fill = "transparent", color = "black", linewidth = 1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = general_size, face = "bold"),
        axis.title = element_text(size = general_size + 1, face = "bold"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = general_size),
        axis.text.y = element_text(size = general_size))

ggsave("~/Documents/plots_project/t70_proportion_class.svg", t70_plot_proportion, width = 130, height = 120, unit = "mm")


data.frame(hm_recall_fnr_2 %>% filter(tool_ref %in% "DeepARG (nt)") %>% 
  arrange(med_recall) %>% filter(new_level %in% (t70 %>% filter(p_bool > 0.9, !n, tool %in% "DeepARG (nt)") %>% select(new_level) %>% pull())))

data.frame(hm_recall_fnr_2 %>% filter(tool_ref %in% "RGI (DIAMOND - nt)") %>% 
             arrange(med_recall) %>% filter(new_level %in% (t70 %>% filter(p_bool> 0.9, !n, tool %in% "RGI (DIAMOND - nt)") %>% select(new_level) %>% pull())))

data.frame(hm_recall_fnr_2 %>% filter(tool_ref %in% "RGI (DIAMOND - nt)") %>% 
             arrange(med_recall) %>% filter(new_level %in% (t50 %>% filter(p_bool> 0.9, !n, tool %in% "RGI (DIAMOND - nt)") %>% select(new_level) %>% pull())))






####

t50_fa <- unigenes %>% 
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

t70_fa <- unigenes %>% 
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
  scale_fill_manual(values = pal_6_q) +  
  xlab("") + 
  ylab("Unigenes") + 
  labs(fill = "In other tool") +
  facet_grid(tool ~  . , scales = "free") +
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.text = element_text(size = general_size),
        panel.border = element_rect(fill = "transparent", color = "black", linewidth = 1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        plot.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.margin = margin(0, 0, 0, 0, unit = "pt"),
        panel.spacing = unit(0, "pt"),
        strip.text = element_text(size = general_size , face = "bold"),
        axis.title = element_text(size = general_size + 1, face = "bold"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = general_size),
        axis.text.y = element_text(size = general_size))
t70_fa_plot


ggsave("~/Documents/plots_project/t70_fa.svg", t70_fa_plot, width = 130, height = 120, unit = "mm")



t70_fa_plot_proportion <- t70_fa %>% filter(n == FALSE) %>% 
  ggplot(aes(x = new_level, y = P)) +
  geom_col(aes(fill=n)) + 
  scale_fill_manual(values = pal_6_q) +  
  xlab("") + 
  ylab("Proportion in gene class") + 
  labs(fill = "In other tool") +
  facet_grid(tool ~  . , scales = "free") +
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.text = element_text(size = general_size),
        panel.border = element_rect(fill = "transparent", color = "black", linewidth = 1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        plot.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.margin = margin(0, 0, 0, 0, unit = "pt"),
        panel.spacing = unit(0, "pt"),
        strip.text = element_text(size = general_size , face = "bold"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = general_size),
        axis.text.y = element_text(size = general_size))
t70_fa_plot_proportion


ggsave("~/Documents/plots_project/t70_fa_proportion_class.svg.svg", t70_fa_plot_proportion, width = 130, height = 120, unit = "mm")


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

core <- readRDS(file = "code_R_analysis/output_abundance_diversity_resistome/core_resistome_no_raw_unique_filter.rds")
pan <- readRDS(file = "code_R_analysis/output_abundance_diversity_resistome/pan_resistome_no_rawunique_filter.rds")


core <- core %>% mutate(habitat = factor(habitat, levels = EN))
pan <- pan %>% mutate(habitat = factor(habitat, levels = EN))
pan <- pan %>% mutate(tool = ifelse(tool == "RGI (DIAMOND nt)", "RGI (DIAMOND - nt)", tool))
core <- core %>% mutate(tool = ifelse(tool == "RGI (DIAMOND nt)", "RGI (DIAMOND - nt)", tool))

sumcore <- core %>% filter(cut %in% 0.5 & cnt > 900, !habitat %in% c( "amplicon", "isolate" ),
                           tool %in% tool_2) %>% ungroup() %>% 
  group_by(new_level, tool, habitat) %>% summarise(unigenes = n_distinct(X))  %>% 
  mutate(tool = factor(tool, levels = tools_levels))

sumcore <- sumcore %>% mutate(tool = factor(tool, levels = tools_levels))

sumpan <- pan %>% ungroup() %>% group_by(tool, habitat, aggregation, gene_class) %>% 
  summarise(md = median(unigenes), mn = mean(unigenes)) %>%
  filter(aggregation == "new_level", !habitat %in% c( "amplicon", "isolate"),
         tool %in% tool_2) %>% mutate(tool = factor(tool, levels = tool_2)) %>% 
  mutate(tool = factor(tool, levels = tools_levels))



#sumpan %>% filter(habitat %in% c("human gut")) %>% group_by(gene_class) %>% summarise(s_mn = median(mn)) %>% arrange(desc(s_mn)) %>% print(n=30)
#sumcore %>% filter(habitat %in% c("human gut")) %>% group_by(new_level) %>% summarise(s_mn = median(unigenes)) %>% arrange(desc(s_mn)) %>% print(n=30)

top20 <- sumpan %>% filter(habitat %in% c("human gut")) %>% ungroup() %>% 
  arrange(desc(mn)) %>% group_by(tool) %>% slice_head(n = 5) %>% ungroup() %>% select(gene_class) %>% distinct() %>% pull()

top20 <- sumcore %>% filter(habitat %in% c("human gut")) %>% ungroup() %>% 
  arrange(desc(unigenes)) %>% group_by(tool) %>% slice_head(n = 5) %>% ungroup() %>% select(new_level) %>% distinct() %>% pull()

top20 <- unique(c(top20, "Cell wall charge"))

top20 <- factor(c(top20, "Other"), levels = c(top20, "Other"))



human_pan <- sumpan %>% filter(habitat %in% c("human gut")) %>% ungroup() %>% 
  mutate(gene_class = ifelse(as.character(gene_class) %in% levels(top20), gene_class, "Other")) %>% 
  mutate(pattern = ifelse(as.character(gene_class) %in% levels(top20)[1:round(length(top20)/2)], "yes", "no")) %>% 
  mutate(gene_class = factor(gene_class, levels = top20)) %>% ungroup() %>% 
  mutate(d = "Pan-resistome") %>% 
  rename(new_level = gene_class, unigenes = md) %>% select(tool, habitat, new_level, unigenes, pattern, d)

human_core <- sumcore %>% filter(habitat %in% c("human gut")) %>% ungroup() %>% 
  mutate(unigenes = as.numeric(unigenes)) %>% 
  mutate(new_level = ifelse(as.character(new_level) %in% levels(top20), new_level, "Other")) %>% 
  mutate(pattern = ifelse(as.character(new_level) %in% levels(top20)[1:round(length(top20)/2)], "yes", "no")) %>% 
  mutate(d = "Core-resistome") %>%
  mutate(new_level = factor(new_level, levels = top20))

human_resistome <- bind_rows(human_pan, human_core) %>% 
  mutate(tool = as.character(tool)) %>% 
  mutate(tool = factor(tool, levels = tools_levels)) 

core_human_class <- human_resistome %>% filter(d %in% "Core-resistome") %>% 
  ungroup() %>% 
  group_by(tool) %>% mutate(N = sum(unigenes), p = unigenes / N) %>% arrange(tool, desc(p))


data.frame(core_human_class %>% slice_head(n = 3))

core_human_class %>% filter(tool %in% "DeepARG (nt)")
core_human_class %>% filter(tool %in% "RGI (DIAMOND - nt)")
core_human_class %>% filter(tool %in% "fARGene (nt)")
core_human_class %>% filter(tool %in% "AMRFinderPlus (aa)")
core_human_class %>% filter(tool %in% "ResFinder (nt)")

p1 <- human_resistome %>% filter(d %in% "Core-resistome") %>% 
  ggplot(aes(x = tool, y = unigenes, fill = tool)) +
  geom_col() +
  theme_minimal() +
  labs(fill = "") +
  scale_x_discrete( labels = tool_label) +
  xlab("") +
  ylab("Total") +
  scale_fill_manual(values = pal_10_q) +
  theme(legend.position = "none",
    legend.text = element_text(size = general_size),
    panel.border = element_rect(fill = "transparent", color = "black", linewidth = 1),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    plot.margin = margin(0, 0, 0, 0, unit = "pt"),
    legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
    legend.margin = margin(0, 0, 0, 0, unit = "pt"),
    panel.spacing = unit(0, "pt"),
    strip.text = element_blank(),
    axis.title = element_text(size = general_size + 1, face = "bold"),
    axis.text.x = element_text(size = general_size),
    axis.text.y = element_blank()) +
  coord_flip()

p2 <- human_resistome %>% filter(d %in% "Core-resistome") %>% 
  group_by(tool, habitat, new_level) %>% 
  summarise(unigenes = sum(unigenes), d = d[1], pattern = pattern[1]) %>%
  ggplot(aes(x = new_level, y = unigenes, fill = tool)) +
  geom_col(position = position_dodge2(preserve = "single")) +
  theme_minimal() +
  labs(fill = "") +
  xlab("") +
  ylab("Unigenes per class") +
  scale_fill_manual(values = c(pal_10_q), labels = tool_label) +
  theme(legend.position = "right",
    panel.border = element_rect(fill = "transparent", color = "black", linewidth = 1),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = general_size),
    axis.text.y = element_text(size = general_size),
    strip.text = element_text(size = general_size, face = "bold"),
    plot.margin = margin(0, 0, 0, 0, unit = "pt"),
    legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
    legend.margin = margin(0, 0, 0, 0, unit = "pt"),
    panel.spacing = unit(0, "pt"),
    legend.text = element_text(size = general_size),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title = element_text(size = general_size + 1, face = "bold"))

p3_core <- p2 +
  annotation_custom(
    grob = ggplotGrob(p1),
    xmin = 5, xmax = 13.5,   # positioning within the data coordinates of p_main
    ymin = 70, ymax = 179
  )

# p3_core <- p2 +
#   annotation_custom(
#     grob = ggplotGrob(p1),
#     xmin = 2, xmax = 3.5,   # positioning within the data coordinates of p_main
#     ymin = 3, ymax = 6
#   )

p3_core

dev.off()
ggsave("~/Documents/plots_project/core_human_gut.svg", p3_core, width = 130, height = 120, unit = "mm")  

ggsave("~/Documents/plots_project/core_human_gut_no_raw_unique_filter.svg", p3_core, width = 130, height = 120, unit = "mm")  



p1 <- human_resistome %>% filter(d %in% "Pan-resistome") %>% 
  ggplot(aes(x = tool, y = unigenes, fill = tool)) +
  geom_col() +
  theme_minimal() +
  labs(fill = "") +
  scale_x_discrete( labels = tool_label) +
  xlab("") +
  ylab("Total") +
  scale_fill_manual(values = pal_10_q) +
  theme(legend.position = "none",
        legend.text = element_text(size = general_size),
        panel.border = element_rect(fill = "transparent", color = "black", linewidth = 1),
        #panel.grid.major.x = element_blank(),
        #panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        plot.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.margin = margin(0, 0, 0, 0, unit = "pt"),
        panel.spacing = unit(0, "pt"),
        strip.text = element_blank(),
        axis.title = element_text(size = general_size + 1, face = "bold"),
        axis.text.x = element_text(size = general_size),
        axis.text.y = element_blank()) +
  #axis.text.y = element_text(size = general_size)) +
  coord_flip()

p2 <- human_resistome %>% filter(d %in% "Pan-resistome") %>% 
  group_by(tool, habitat, new_level) %>% 
  summarise(unigenes = sum(unigenes), d = d[1], pattern = pattern[1]) %>%
  ggplot(aes(x = new_level, y = unigenes, fill = tool)) +
  geom_col(position = position_dodge2(preserve = "single")) +
  theme_minimal() +
  labs(fill = "") +
  xlab("") +
  ylab("Unigenes per class") +
  scale_fill_manual(values = c(pal_10_q), labels = tool_label) +
  theme(legend.position = "right",
        panel.border = element_rect(fill = "transparent", color = "black", linewidth = 1),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = general_size),
        axis.text.y = element_text(size = general_size),
        strip.text = element_text(size = general_size, face = "bold"),
        plot.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.margin = margin(0, 0, 0, 0, unit = "pt"),
        panel.spacing = unit(0, "pt"),
        legend.text = element_text(size = general_size),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title = element_text(size = general_size + 1, face = "bold"))

p3_pan <- p2 +
  annotation_custom(
    grob = ggplotGrob(p1),
    xmin = 5, xmax = 13,   # positioning within the data coordinates of p_main
    ymin = 3500, ymax = 8000
  )

p3_pan <- p2 +
  annotation_custom(
    grob = ggplotGrob(p1),
    xmin = 5, xmax = 13,   # positioning within the data coordinates of p_main
    ymin = 2500, ymax = 6000
  )

p3_pan

dev.off()

ggsave("~/Documents/plots_project/pan_human_gut.svg", p3_pan, width = 130, height = 120, unit = "mm") 

ggsave("~/Documents/plots_project/pan_human_gut_no_raw_unique_filter.svg", p3_pan, width = 130, height = 120, unit = "mm") 



# human_resistome %>% filter(d %in% "Core-resistome") %>% 
#   ggplot(aes(x = tool, y = unigenes, fill = new_level, pattern = pattern)) +
#     geom_col_pattern(pattern_density = 0.1, pattern_spacing = 0.02, pattern_size = 0.05) +
#     theme_minimal() +
#     labs(fill = "") +
#     xlab("") +
#     scale_fill_manual(values = c(pal_10_d,pal_10_d)) +
#     scale_pattern_manual(values = c("stripe", "none")) + 
#     theme(
#         panel.border = element_rect(fill = "transparent", color = "black", linewidth = 1),
#         axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = general_size),
#         axis.text.y = element_text(size = general_size),
#         strip.text = element_text(size = general_size, face = "bold"),
#         panel.grid.major.x = element_blank(),
#         panel.grid.minor.x = element_blank(),
#         axis.title = element_text(size = general_size + 1, face = "bold"))



hbt <- c("human gut", "wastewater", "freshwater", "soil", "cat gut", "dog gut", "pig gut")

core_environment <- sumcore %>% filter(habitat %in% hbt) %>% ungroup() %>% group_by(habitat, tool) %>% summarise(unigenes = sum(unigenes)) %>% 
  ggplot(aes(x = tool, y = unigenes, fill = habitat)) +
  geom_col(position = position_dodge2(preserve = "single", width = 0.9), width = .6) +
  theme_minimal() +
  labs(fill = "") +
  scale_x_discrete( labels = tool_label) +
  xlab("") +
  ylab("Core-resistome") +
  scale_fill_manual(values = pal_8_q) +
  theme(legend.position = "right",
        legend.text = element_text(size = general_size),
        panel.border = element_rect(fill = "transparent", color = "black", linewidth = 1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        plot.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.margin = margin(0, 0, 0, 0, unit = "pt"),
        panel.spacing = unit(0, "pt"),
        strip.text = element_blank(),
        axis.title = element_text(size = general_size + 1, face = "bold"),
        axis.text.x = element_text(angle = 90, size = general_size),
        axis.text.y = element_text(size = general_size)) 


ggsave("~/Documents/plots_project/core_environment_no_raw_unique_filter.svg", core_environment, width = 130, height = 120, unit = "mm") 

pan_environment <- sumpan %>% filter(habitat %in% hbt) %>% ungroup() %>% group_by(habitat, tool) %>% summarise(unigenes = sum(md)) %>% ungroup() %>% 
  ggplot(aes(x = tool, y = unigenes, fill = habitat)) +
  geom_col(position = position_dodge2(preserve = "single", width = 0.9), width = .6) +
  theme_minimal() +
  labs(fill = "") +
  scale_x_discrete( labels = tool_label) +
  #facet_grid(. ~ tool,  scales="free_x", switch="x") + 
  xlab("") +
  ylab("Core-resistome") +
  scale_fill_manual(values = pal_8_q) +
  theme(legend.position = "right",
        legend.text = element_text(size = general_size),
        panel.border = element_rect(fill = "transparent", color = "black", linewidth = 1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        plot.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.margin = margin(0, 0, 0, 0, unit = "pt"),
        panel.spacing = unit(1, "pt"),
        strip.text = element_blank(),
        axis.title = element_text(size = general_size + 1, face = "bold"),
        axis.text.x = element_text(angle = 90, size = general_size),
        axis.text.y = element_text(size = general_size)) 



ggsave("~/Documents/plots_project/pan_environment_no_raw_unique_filter.svg", pan_environment, width = 130, height = 120, unit = "mm") 


##################################################################################################################################################################
##################################################################################################################################################################
##################################################################################################################################################################
##################################################################################################################################################################
##################################################################################################################################################################
##################################################################################################################################################################



sets0 <- core %>% filter(cut %in% 0.5 & cnt > 900, !habitat %in% c( "amplicon", "isolate" ), tool %in% tool_2) %>% 
  ungroup() %>% mutate(tool = factor(tool, levels = tool_2)) %>% 
  group_by(habitat) %>%
  summarise(query = list(X), .groups = "drop")

sets1 <- core %>% filter(cut %in% 0.5 & cnt > 900, !habitat %in% c( "amplicon", "isolate" ), tool %in% tool_2) %>% 
  ungroup() %>% mutate(tool = factor(tool, levels = tool_2)) %>% 
  group_by(habitat) %>% 
  group_by( tool, habitat) %>%
  summarise(query = list(X), .groups = "drop")

# Pairwise combinations 
pairwise <- sets1 %>%
  group_by(tool) %>%
  summarise(pairs = list(expand_grid(habitat_ref = habitat, habitat_comp = habitat)), .groups = "drop") %>%
  unnest(pairs)

core_class_tool <- pairwise %>%
  left_join(sets1, by = c("tool", "habitat_ref" = "habitat")) %>%
  rename(qc_ref = query) %>%
  left_join(sets1, by = c("tool", "habitat_comp" = "habitat")) %>%
  rename(qc_comp = query) %>%
  left_join(sets0, by = c( "habitat_ref" = "habitat")) 


intersect_core <- function(qc_ref, qc_comp){
  A <-  unlist(qc_ref)
  B <- unlist(qc_comp)
  d <- intersect(A, B)
  r <- ifelse(length(A) == 0 | length(B) == 0, NA, length(d))
  return(r)
}


core_class_tool <- core_class_tool %>% rowwise() %>% mutate(shared = intersect_core(qc_ref, qc_comp))
core_class_tool <- core_class_tool %>% 
  mutate(ref_n_class = length(qc_ref), comp_n_class = length(qc_comp))
core_class_tool <- core_class_tool %>% 
  ungroup() %>% mutate( shared = ifelse(habitat_ref == habitat_comp, 0, shared)) %>% 
  select(-c(qc_ref, qc_comp))

core_class_tool %>% ungroup() %>% arrange(desc(recall))


ggplot(core_class_tool %>% filter(tool %in% c("RGI (DIAMOND - nt)", "DeepARG (nt)", "fARGene (nt)", "ResFinder (nt)")) %>%
         filter(as.numeric(habitat_ref) < as.numeric(habitat_comp)), aes(x = habitat_ref, y = habitat_comp, fill = shared)) +
  geom_tile() +
  scale_fill_viridis_c(na.value = 0) + 
  facet_wrap(. ~ tool) +
  theme_minimal() +
  xlab("") +
  ylab("") +
  labs(fill = "Shared core genes") +
  theme(axis.title.y = element_blank(),
        legend.position = "bottom",
        panel.border = element_rect(fill = "transparent", color = "black", linewidth = 1),
        axis.text.x = element_text(angle = 90, size = general_size),
        axis.text.y = element_text(size = general_size),
        strip.text = element_text(size = general_size, face = "bold"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title = element_text(size = general_size + 1, face = "bold"))


ggplot(core_class_tool %>% filter(tool %in% "RGI (DIAMOND - nt)"), aes(x = habitat_ref, y = habitat_comp, fill = shared)) +
  geom_tile() +
  scale_fill_viridis_c(na.value = 0) + 
  facet_grid(. ~ tool, scales = "free_x") +
  theme_minimal() +
  xlab("") +
  ylab("") +
  labs(fill = "Shared core genes") +
  theme(axis.title.y = element_blank(),
        legend.position = "bottom",
        panel.border = element_rect(fill = "transparent", color = "black", linewidth = 1),
        axis.text.x = element_text(angle = 90, size = general_size),
        axis.text.y = element_text(size = general_size),
        strip.text = element_text(size = general_size, face = "bold"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title = element_text(size = general_size + 1, face = "bold"))


ggplot(core_class_tool %>% filter(tool %in% "fARGene (nt)"), aes(x = habitat_ref, y = habitat_comp, fill = shared)) +
  geom_tile() +
  scale_fill_viridis_c(na.value = 0) + 
  facet_grid(. ~ tool, scales = "free_x") +
  theme_minimal() +
  xlab("") +
  ylab("") +
  labs(fill = "Shared core genes") +
  theme(axis.title.y = element_blank(),
        legend.position = "bottom",
        panel.border = element_rect(fill = "transparent", color = "black", linewidth = 1),
        axis.text.x = element_text(angle = 90, size = general_size),
        axis.text.y = element_text(size = general_size),
        strip.text = element_text(size = general_size, face = "bold"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title = element_text(size = general_size + 1, face = "bold"))






top20 <- sumcore %>% filter(habitat %in% c("human gut")) %>% ungroup() %>% 
  arrange(desc(unigenes)) %>% group_by(tool) %>% slice_head(n = 5) %>% ungroup() %>% select(new_level) %>% distinct() %>% pull()

top20 <- unique(c(top20, "Cell wall charge"))

top20 <- factor(c(top20, "Other"), levels = c(top20, "Other"))



pig_pan <- sumpan %>% filter(habitat %in% c("pig gut")) %>% ungroup() %>% 
  mutate(gene_class = ifelse(as.character(gene_class) %in% levels(top20), gene_class, "Other")) %>% 
  mutate(pattern = ifelse(as.character(gene_class) %in% levels(top20)[1:round(length(top20)/2)], "yes", "no")) %>% 
  mutate(gene_class = factor(gene_class, levels = top20)) %>% ungroup() %>% 
  mutate(d = "Pan-resistome") %>% 
  rename(new_level = gene_class, unigenes = md) %>% select(tool, habitat, new_level, unigenes, pattern, d)

pig_core <- sumcore %>% filter(habitat %in% c("pig gut")) %>% ungroup() %>% 
  mutate(unigenes = as.numeric(unigenes)) %>% 
  mutate(new_level = ifelse(as.character(new_level) %in% levels(top20), new_level, "Other")) %>% 
  mutate(pattern = ifelse(as.character(new_level) %in% levels(top20)[1:round(length(top20)/2)], "yes", "no")) %>% 
  mutate(d = "Core-resistome") %>%
  mutate(new_level = factor(new_level, levels = top20))

pig_resistome <- bind_rows(pig_pan, pig_core) %>% 
  mutate(tool = as.character(tool)) %>% 
  mutate(tool = factor(tool, levels = tool_2)) 


p1 <- pig_resistome %>% filter(d %in% "Core-resistome") %>% 
  ggplot(aes(x = tool, y = unigenes, fill = tool)) +
  geom_col() +
  theme_minimal() +
  labs(fill = "") +
  scale_x_discrete( labels = tool_label) +
  xlab("") +
  ylab("Total") +
  scale_fill_manual(values = pal_10_q) +
  theme(legend.position = "none",
        legend.text = element_text(size = general_size),
        panel.border = element_rect(fill = "transparent", color = "black", linewidth = 1),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        plot.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.margin = margin(0, 0, 0, 0, unit = "pt"),
        panel.spacing = unit(0, "pt"),
        strip.text = element_blank(),
        axis.title = element_text(size = general_size + 1, face = "bold"),
        axis.text.x = element_text(size = general_size),
        axis.text.y = element_blank()) +
  coord_flip()

p2 <- pig_resistome %>% filter(d %in% "Core-resistome") %>% 
  group_by(tool, habitat, new_level) %>% 
  summarise(unigenes = sum(unigenes), d = d[1], pattern = pattern[1]) %>%
  ggplot(aes(x = new_level, y = unigenes, fill = tool)) +
  geom_col(position = position_dodge2(preserve = "single")) +
  theme_minimal() +
  labs(fill = "") +
  xlab("") +
  ylab("Unigenes per class") +
  scale_fill_manual(values = c(pal_10_q), labels = tool_label) +
  theme(legend.position = "right",
        panel.border = element_rect(fill = "transparent", color = "black", linewidth = 1),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = general_size),
        axis.text.y = element_text(size = general_size),
        strip.text = element_text(size = general_size, face = "bold"),
        plot.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.margin = margin(0, 0, 0, 0, unit = "pt"),
        panel.spacing = unit(0, "pt"),
        legend.text = element_text(size = general_size),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title = element_text(size = general_size + 1, face = "bold"))

p3_core <- p2 +
  annotation_custom(
    grob = ggplotGrob(p1),
    xmin = 5, xmax = 13.5,   # positioning within the data coordinates of p_main
    ymin = 350, ymax = 1100
  )

p3_core

ggsave("~/Documents/plots_project/core_pig_gut_no_raw_unique_filter.svg", p3_core, width = 130, height = 120, unit = "mm")  




p1 <- pig_resistome %>% filter(d %in% "Pan-resistome") %>% 
  ggplot(aes(x = tool, y = unigenes, fill = tool)) +
  geom_col() +
  theme_minimal() +
  labs(fill = "") +
  scale_x_discrete( labels = tool_label) +
  xlab("") +
  ylab("Total") +
  scale_fill_manual(values = pal_10_q) +
  theme(legend.position = "none",
        legend.text = element_text(size = general_size),
        panel.border = element_rect(fill = "transparent", color = "black", linewidth = 1),
        #panel.grid.major.x = element_blank(),
        #panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        plot.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.margin = margin(0, 0, 0, 0, unit = "pt"),
        panel.spacing = unit(0, "pt"),
        strip.text = element_blank(),
        axis.title = element_text(size = general_size + 1, face = "bold"),
        axis.text.x = element_text(size = general_size),
        axis.text.y = element_blank()) +
  #axis.text.y = element_text(size = general_size)) +
  coord_flip()

p2 <- pig_resistome %>% filter(d %in% "Pan-resistome") %>% 
  group_by(tool, habitat, new_level) %>% 
  summarise(unigenes = sum(unigenes), d = d[1], pattern = pattern[1]) %>%
  ggplot(aes(x = new_level, y = unigenes, fill = tool)) +
  geom_col(position = position_dodge2(preserve = "single")) +
  theme_minimal() +
  labs(fill = "") +
  xlab("") +
  ylab("Unigenes per class") +
  scale_fill_manual(values = c(pal_10_q), labels = tool_label) +
  theme(legend.position = "right",
        panel.border = element_rect(fill = "transparent", color = "black", linewidth = 1),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = general_size),
        axis.text.y = element_text(size = general_size),
        strip.text = element_text(size = general_size, face = "bold"),
        plot.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.margin = margin(0, 0, 0, 0, unit = "pt"),
        panel.spacing = unit(0, "pt"),
        legend.text = element_text(size = general_size),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title = element_text(size = general_size + 1, face = "bold"))

p3_pan <- p2 +
  annotation_custom(
    grob = ggplotGrob(p1),
    xmin = 5, xmax = 13,   # positioning within the data coordinates of p_main
    ymin = 3500, ymax = 8000
  )


dev.off()

ggsave("~/Documents/plots_project/pan_pig_gut.svg", p3_pan, width = 130, height = 120, unit = "mm") 

################


sumpan2 <- pan %>% ungroup() %>% group_by(tool, habitat, aggregation, epoch) %>% 
  summarise(s = sum(unigenes)) %>%
  ungroup() %>% group_by(tool, habitat, aggregation) %>% 
  summarise(md = median(s), mn = mean(s), sd = sd(s)) %>%
  filter(aggregation == "new_level", !habitat %in% c( "amplicon", "isolate"),
         tool %in% tool_2) %>% mutate(tool = factor(tool, levels = tool_2)) %>% 
  mutate(tool = factor(tool, levels = tools_levels))


#sumpan %>% filter(aggregation %in% "new_level", habitat %in% c("human gut", "wastewater", "pig gut", "marine", "soil", "built-environment")) %>% 
sumpan2 %>% filter(aggregation %in% "new_level") %>% #, habitat %in% c("human gut", "wastewater", "pig gut", "marine", "soil", "built-environment")) %>% 
  ggplot(aes(x = tool, y = mn, fill = habitat)) +
  geom_col(position = position_dodge2(preserve = "single")) +
  geom_errorbar(aes(ymin = mn-sd, ymax = mn+sd), width=.2,
                position = position_dodge(.9)) +
  theme_minimal() +
  labs(fill = "") +
  xlab("") +
  ylab("Pan-resistome per class") +
  scale_fill_manual(values = c(pal_10_q,pal_10_q)) +
  theme(legend.position = "right",
        panel.border = element_rect(fill = "transparent", color = "black", linewidth = 1),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = general_size),
        axis.text.y = element_text(size = general_size),
        strip.text = element_text(size = general_size, face = "bold"),
        plot.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.margin = margin(0, 0, 0, 0, unit = "pt"),
        panel.spacing = unit(0, "pt"),
        legend.text = element_text(size = general_size),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title = element_text(size = general_size + 1, face = "bold"))


sumpan2 %>% filter(tool %in% c("ResFinder (nt)", "ABRicate (ResFinder - nt)")) %>% arrange(habitat) %>% group_by(habitat) %>% summarise(d = mn[1] - mn[2], p = (mn[1] - mn[2])/mn[1] )
sumpan2 %>% filter(tool %in% c("AMRFinderPlus (aa)", "ABRicate (NCBI - nt)")) %>% arrange(habitat) %>% group_by(habitat) %>% summarise(d = mn[1] - mn[2], p = (mn[1] - mn[2])/mn[1] ) %>% mutate(m = mean(p))
sumpan2 %>% filter(tool %in% c("RGI (DIAMOND - nt)", "ABRicate (CARD - nt)")) %>% arrange(habitat) %>% group_by(habitat) %>% summarise(d = mn[1] - mn[2], p = (mn[1] - mn[2])/mn[1] )


sumcore2 <- core %>% filter(cut %in% 0.5 & cnt > 900, !habitat %in% c( "amplicon", "isolate" ),
                           tool %in% tool_2) %>% ungroup() %>% 
  group_by( tool, habitat) %>% summarise(unigenes = n_distinct(X))  %>% 
  mutate(tool = factor(tool, levels = tools_levels))


sumcore2 %>% filter(tool %in% c("ResFinder (nt)", "ABRicate (ResFinder - nt)")) %>% arrange(habitat) %>% group_by(habitat) %>% summarise(d = unigenes[1] - unigenes[2], p = (unigenes[1] - unigenes[2])/unigenes[1] )
sumcore2 %>% filter(tool %in% c("AMRFinderPlus (aa)", "ABRicate (NCBI - nt)")) %>% arrange(habitat) %>% group_by(habitat) %>% summarise(d = unigenes[1] - unigenes[2], p = (unigenes[1] - unigenes[2])/unigenes[1] )
sumcore2 %>% filter(tool %in% c("RGI (DIAMOND - nt)", "ABRicate (CARD - nt)")) %>% arrange(habitat) %>% group_by(habitat) %>% summarise(d = unigenes[1] - unigenes[2], p = (unigenes[1] - unigenes[2])/unigenes[1] )

sumcore2 %>% filter(habitat %in% c("human gut"))
sumcore2 %>% filter(habitat %in% c("pig gut"))
sumcore2 %>% filter(habitat %in% c("cat gut"))
sumcore2 %>% filter(habitat %in% c("dog gut"))
sumcore2 %>% filter(habitat %in% c("mouse gut"))
sumcore2 %>% filter(habitat %in% c("human skin"))
sumcore2 %>% filter(habitat %in% c("human oral"))
sumcore2 %>% filter(habitat %in% c("human nose"))
sumcore2 %>% filter(habitat %in% c("human vagina"))
sumcore2 %>% filter(habitat %in% c("wastewater"))
sumcore2 %>% filter(habitat %in% c("freshwater"))
sumcore2 %>% filter(habitat %in% c("marine"))
sumcore2 %>% filter(habitat %in% c("soil"))
sumcore2 %>% filter(habitat %in% c("built-environment"))


sumpan2 %>% ungroup() %>% arrange(habitat) %>%  
  group_by(habitat) %>%  
  mutate(deeparg = mn[tool == "DeepARG (nt)"][1],
         rgi = mn[tool == "RGI (DIAMOND - nt)"][1]) %>%
  mutate(dp = deeparg / mn, drgi = rgi / mn) %>% 
  filter(!tool %in% c("DeepARG (nt)", "RGI (DIAMOND - nt)")) %>% 
  summarise(m_da = mean(dp), m_rgi = mean(drgi))

sumcore2 %>% ungroup() %>% arrange(habitat) %>%  
  group_by(habitat) %>%  
  mutate(deeparg = unigenes[tool == "DeepARG (nt)"][1],
         rgi = unigenes[tool == "RGI (DIAMOND - nt)"][1]) %>%
  mutate(dp = deeparg / unigenes, drgi = rgi / unigenes) %>% 
  filter(!tool %in% c("DeepARG (nt)", "RGI (DIAMOND - nt)")) %>% 
  summarise(m_da = mean(dp), m_rgi = mean(drgi))



sumcore2 %>% #, habitat %in% c("human gut", "wastewater", "pig gut", "marine", "soil", "built-environment")) %>% 
  ggplot(aes(x = tool, y = unigenes, fill = habitat)) +
  geom_col(position = position_dodge2(preserve = "single")) +
  theme_minimal() +
  labs(fill = "") +
  xlab("") +
  ylab("Core-resistome") +
  scale_fill_manual(values = c(pal_10_q,pal_10_q)) +
  theme(legend.position = "right",
        panel.border = element_rect(fill = "transparent", color = "black", linewidth = 1),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = general_size),
        axis.text.y = element_text(size = general_size),
        strip.text = element_text(size = general_size, face = "bold"),
        plot.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.margin = margin(0, 0, 0, 0, unit = "pt"),
        panel.spacing = unit(0, "pt"),
        legend.text = element_text(size = general_size),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title = element_text(size = general_size + 1, face = "bold"))



top20 <- sumcore %>% filter(habitat %in% c("human gut")) %>% ungroup() %>% 
  arrange(desc(unigenes)) %>% group_by(tool) %>% slice_head(n = 5) %>% ungroup() %>% select(new_level) %>% distinct() %>% pull()
top20 <- unique(c(top20, "Cell wall charge"))

sumcore_human_class <- sumcore %>% filter(habitat %in% "human gut") %>% 
  ungroup() %>% group_by(tool) %>% mutate(total = sum(unigenes)) %>% 
  filter( new_level %in% top20) %>% mutate(p = unigenes / total )

data.frame(sumcore_human_class)


#


sets0 <- core %>% filter(cut %in% 0.5 & cnt > 900, habitat %in% c( "human gut" ), tool %in% tool_2) %>%
  group_by(tool) %>%
  summarise(query = list(X), .groups = "drop")

sets1 <- core %>% filter(cut %in% 0.5 & cnt > 900, habitat %in% c( "human gut" ), tool %in% tool_2) %>%
  group_by(new_level, tool) %>%
  summarise(query = list(X), .groups = "drop")

# Pairwise combinations 
pairwise <- sets1 %>%
  group_by(new_level) %>%
  summarise(pairs = list(expand_grid(tool_ref = tool, tool_comp = tool)), .groups = "drop") %>%
  unnest(pairs)

JI_core_class <- pairwise %>%
  left_join(sets1, by = c("new_level", "tool_ref" = "tool")) %>%
  rename(qc_ref = query) %>%
  left_join(sets1, by = c("new_level", "tool_comp" = "tool")) %>%
  rename(qc_comp = query) %>%
  left_join(sets0, by = c( "tool_ref" = "tool")) %>%
  rename(q_ref = query) %>% 
  left_join(sets0, by = c( "tool_comp" = "tool")) %>%
  rename(q_comp = query)


new_intersect <- function(qc_ref, q_ref, qc_comp, q_comp){
  A <-  unlist(qc_ref)
  B <- unlist(q_ref)
  C <- unlist(qc_comp)
  D <- unlist(q_comp)
  x <- intersect(A, setdiff(D, C))
  y <- setdiff(C, intersect(C, setdiff(B, A)))
  z <- union(x, y)
  d <- intersect(A, z)
  r <- ifelse(length(z) == 0, NA, length(d) / length(z))
  return(r)
}

new_difference <- function(qc_ref, q_ref, qc_comp, q_comp){
  A <-  unlist(qc_ref)
  B <- unlist(q_ref)
  C <- unlist(qc_comp)
  D <- unlist(q_comp)
  x <- intersect(A, setdiff(D, C))
  y <- setdiff(C, intersect(C, setdiff(B, A)))
  z <- union(x, y)
  d <- setdiff(A, z)
  r <- ifelse(length(A) == 0, NA, length(d) / length(A))
  return(r)
}

JI_core_class <- pairwise %>%
  left_join(sets1, by = c("new_level", "tool_ref" = "tool")) %>%
  rename(qc_ref = query) %>%
  left_join(sets1, by = c("new_level", "tool_comp" = "tool")) %>%
  rename(qc_comp = query) %>%
  left_join(sets0, by = c( "tool_ref" = "tool")) %>%
  rename(q_ref = query) %>% 
  left_join(sets0, by = c( "tool_comp" = "tool")) %>%
  rename(q_comp = query)

JI_core_class <- JI_core_class %>% rowwise() %>% mutate(recall = new_intersect(qc_ref, q_ref, qc_comp, q_comp))
JI_core_class <- JI_core_class %>% rowwise() %>% mutate(fnr = new_difference(qc_ref, q_ref, qc_comp, q_comp))
JI_core_class <- JI_core_class %>% 
  mutate(ref_n_class = length(qc_ref), comp_n_class = length(qc_comp), ref_n_all = length(q_ref), comp_n_all = length(q_comp))
JI_core_class <- JI_core_class %>% 
  ungroup() %>% filter(tool_ref != tool_comp) %>% 
  select(-c(qc_ref, qc_comp, q_ref, q_comp))


data.frame(JI_core_class %>% filter(new_level %in% top20) )





###

sets0 <- core %>% filter(cut %in% 0.5 & cnt > 900, habitat %in% c( "human gut" ), tool %in% tool_2) %>%
  group_by(tool) %>%
  summarise(query = list(X), .groups = "drop")

pairwise <- expand_grid(tool_ref = sets0$tool, tool_comp = sets0$tool) %>% filter(tool_ref != tool_comp)

new_intersect2 <- function( q_ref, q_comp){
  A <- unlist(q_ref)
  B <- unlist(q_comp)
  d <- intersect(A, B)
  r <- ifelse(length(d) == 0, 0, length(d))
  return(r)
}

JI_core_overlap <- pairwise %>%
  left_join(sets0, by = c("tool_ref" = "tool")) %>%
  rename(q_ref = query) %>%
  left_join(sets0, by = c("tool_comp" = "tool")) %>%
  rename(q_comp = query)

JI_core_overlap <- JI_core_overlap %>% rowwise() %>% mutate(inter = new_intersect2(q_ref, q_comp)) %>% 
  mutate(ref_n = length(q_ref), comp_n = length(q_comp)) %>% select(-c(q_ref, q_comp)) 



JI_core_overlap <- JI_core_overlap %>% mutate(recall = inter / comp_n, diff = (ref_n - inter)/ref_n)

data.frame(JI_core_overlap)

JI_core_overlap %>% filter(tool_comp %in% "fARGene (nt)") %>% select(-c(diff))
JI_core_overlap %>% filter(tool_ref %in% "fARGene (nt)") 



##### abundance Radial per tool


## factors for new_level, we take the highest abundance per ontology by tool

# factor_aro <- abundance %>% ungroup() %>% filter(aggregation %in% "ARO") %>%
#   group_by(aggregation, tool, gene ) %>% summarise(total = sum(normed10m)) %>%
#   ungroup() %>% arrange(aggregation, tool, desc(total)) %>% 
#   group_by(aggregation, tool, gene) %>%  ungroup() %>% select(gene) %>% distinct() %>% pull()
#   
# factor_parent <- abundance %>% ungroup() %>% filter(aggregation %in% "parent_description") %>%
#   group_by(aggregation, tool, gene ) %>% summarise(total = sum(normed10m)) %>%
#     ungroup() %>% arrange(aggregation, tool, desc(total)) %>% 
#     group_by(aggregation, tool, gene) %>%  ungroup() %>% select(gene) %>% distinct() %>% pull()
#   
# factor_new_level <- abundance %>% ungroup() %>% filter(aggregation %in% "new_level") %>%
#   group_by(aggregation, tool, gene ) %>% summarise(total = sum(normed10m)) %>%
#     ungroup() %>% arrange(aggregation, tool, desc(total)) %>% 
#     group_by(aggregation, tool, gene) %>%  ungroup() %>% select(gene) %>% distinct() %>% pull()

# this is mainly for the circular / radial plots
# factor_new_level2 <- c(factor_new_level[seq(1, length(factor_new_level), by = 3)],
#                        factor_new_level[seq(1, length(factor_new_level), by = 3) + 1],
#                        factor_new_level[seq(1, length(factor_new_level), by = 3) + 2])



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
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))



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
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))


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
