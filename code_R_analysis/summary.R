library(dplyr)
library(ggplot2)
library(ggVennDiagram)
library(gridExtra)
library(tidyverse)
library(RColorBrewer)
library(ggpattern)
library(grid)
library(eulerr)

# 130mm 1 col 185 2 col
# 180mm 1 col 210 2 col
# 220mm 1 col 225 2 col

# retrieve legends from plots 

g_legend <- function(a.gplot){
  tmp <- ggplotGrob(a.gplot)
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}


# d is divergent / q is qualitative
pal_8_q <- brewer.pal(8, "Dark2")
pal_10_q <- brewer.pal(10, "Paired")
pal_11_q <- brewer.pal(11, "Paired")
pal_6_q <- brewer.pal(6, "Paired")
pal_10_d <- brewer.pal(10, "BrBG")
pal_4_d <- brewer.pal(4, "BrBG")
pal_6_d1 <- brewer.pal(6, "PiYG")
pal_11_d <- brewer.pal(11, "BrBG")
pal_9_d <- brewer.pal(9, "BrBG")

# pallet for 15 tool + format repeated by tool
pal15_rep <- c(rep(pal_10_d[1], 2), rep(pal_10_d[2], 3), rep(pal_10_d[3], 2),
               pal_10_d[4], rep(pal_10_d[5], 2), pal_10_d[6:10])


pal_10_q <- pal_10_q[c(6,2,8,7,4,3,10,9,1,5)]

# size for font in plots 
general_size <- 5

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
                  "AMRFinderPlus (aa)", "ABRicate (NCBI - nt)",
                  "ResFinder (nt)", "ABRicate (ResFinder - nt)", 
                  "ABRicate (ARG-ANNOT - nt)", "ABRicate (MEGARES - nt)")

tool_2 <- c("DeepARG (nt)", "RGI (DIAMOND - nt)", "fARGene (nt)", "AMRFinderPlus (aa)", "ResFinder (nt)",
            "ABRicate (ARG-ANNOT - nt)", "ABRicate (CARD - nt)", "ABRicate (MEGARES - nt)", "ABRicate (NCBI - nt)",
            "ABRicate (ResFinder - nt)")

tool_label <- c(expression({DeepARG~""^"a"}), expression({fARGene~""^"a"}), 
                expression({RGI~""^"a,b"}),  expression({CARD~""^"a,d"}), 
                expression({AMRFinderPlus~""^"c"}), expression({NCBI~""^"a,d"}),
                expression(ResFinder~{""^"a"}), expression({ResFinder~""^"a,d"}),
                expression({ARG-ANNOT~""^"a,d"}), 
                expression({MEGARES~""^"a,d"}))


abundance <- abundance %>% mutate( tool = 
                     ifelse(tool == "RGI (BLAST nt)", "RGI (BLAST - nt)",
                     ifelse(tool == "RGI (BLAST aa)", "RGI (BLAST - aa)",
                     ifelse(tool == "RGI (DIAMOND aa)", "RGI (DIAMOND - aa)",
                     ifelse(tool == "RGI (DIAMOND nt)", "RGI (DIAMOND - nt)", tool)))))

# environments that we are not interested in
not_env <- c("amplicon", "isolate", "built-environment" )


# changing habitats and tools to factor

abundance <- abundance %>% mutate(habitat = factor(habitat, levels = EN),
                                  habitat2 = factor(SO[habitat], levels = c("Humans","Mammals","Wastewater", "Freshwater","Soil", "Marine", "Other")),
                                  tool = factor(tool, levels = tools_levels))

abundance <- abundance %>% mutate(location = ifelse(habitat2 %in% c("Humans","Mammals","Wastewater","Built-environment"), "Human-related","External"))

abundance <- abundance %>% mutate(location = factor(location, levels = c("Human-related","External")))


# keep original abundance data frame in abundance0
abundance0 <- abundance

# filter for analysis of environments and tools wanted 

abundance <- abundance %>% filter(tool %in% tool_2 & !habitat %in% not_env)
abundance <- abundance %>% filter(aggregation %in% "new_level")
abundance <- abundance %>% mutate(habitat2 = factor(as.character(habitat2), levels = c("Humans","Mammals","Wastewater", "Freshwater","Soil", "Marine")))


# SUMMARIES
# SUMMARIES
# SUMMARIES

# unigenes captured with the tools

unigenes0 <- do.call(rbind, lapply(lst, function(x) x[,c("query","tool", "ARO", "parent", "parent_description", "new_level", "id")])) 
unigenes0 %>% select(query) %>% distinct() %>% summarise(n = n())

unigenes <- do.call(rbind, lapply(lst, function(x) x[,c("query","tool", "ARO", "parent", "parent_description", "new_level", "id")])) 
unigenes <- unigenes %>% 
  filter(tool %in% tools_levels) %>% 
  mutate(tool = factor(tool, levels = tools_levels))

# unigenes per tool
unigenes  %>% group_by(tool) %>% summarise(n = n_distinct(query)) %>% ungroup() %>% arrange(n)

plot_count_genes_tool <- unigenes  %>% filter(tool %in% tool_2) %>% 
  ggplot(aes( x = tool)) +
  geom_bar(aes(fill = tool), color = "black", linewidth = 0.2) +
  scale_fill_manual(values = pal_10_q) +
  scale_x_discrete( labels = tool_label) +
  theme_minimal() +
  ylab("Genes") +
  xlab("") +
  ggtitle("A") +
  labs(fill = "") +
  #scale_y_continuous(expand = c(0, 0), limits = c(0, 126000)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 128000),  breaks = c(0, 25000, 50000,75000,100000,125000), labels = c("0","25,000", "50,000", "75,000","100,000","125,000")) +
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
    title = element_text(size = general_size + 2, face = "bold"),
    axis.title = element_text(size = general_size + 1, face = "bold"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = general_size),
    axis.text.y = element_text(size = general_size))

plot_count_genes_tool

dev.off()
ggsave("~/Documents/plots_project2/count_genes_per_tool_2.svg", plot_count_genes_tool  , width = 90, height = 70, unit = "mm")
dev.off()  

plot_count_genes_tool

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
  arrange(desc(n)) %>% mutate(p = n/N) %>% slice_head(n=11)  %>% 
  select(new_level)  %>% pull()

top10class_tool <- unigenes %>% 
  mutate(N = n_distinct(query)) %>% 
  group_by(tool, new_level) %>% summarise(n = n_distinct(query)) %>% mutate(N = sum(n)) %>%
  arrange(desc(n)) %>% mutate(p = n/N) %>% slice_head(n=11) 

top10class <- top10class_tool %>% filter(p >= 0.05) %>% ungroup() %>% select(new_level) %>% distinct()%>% pull()

top10class <- factor(c(top10class, "Other"), 
                     levels = c(top10class, "Other"))

top10class_plot <- bind_rows(proportion_new_level_tool %>% filter(new_level %in% top10class) %>% 
            select(new_level, tool, p) %>% ungroup() %>% group_by(tool) %>% mutate(P = sum(p)),
          proportion_new_level_tool %>% filter(!new_level %in% top10class, tool %in% tool_2) %>% 
            select(new_level, tool, p) %>% ungroup() %>% group_by(tool) %>% 
            summarise(p = sum(p)) %>% mutate(new_level = "Other")) %>% 
  mutate(tool = factor(tool, levels = tools_levels),
         new_level = factor(new_level, levels = top10class)) %>% 
  ggplot(aes( x = tool, y = p)) +
  geom_col(aes(fill = new_level), color = "black", linewidth = 0.2) +
  scale_fill_manual(values = c(pal_10_d, pal_6_d1)) +
  scale_x_discrete(labels = tool_label) +
  theme_minimal() +
  ylab("Proportion") +
  xlab("") +
  labs(fill = "") +
  ggtitle("B") + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.0001)) +
  guides(fill = guide_legend(ncol = 2)) +
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
    title = element_text(size = general_size + 2, face = "bold"),
    axis.title = element_text(size = general_size + 1, face = "bold"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = general_size),
    axis.text.y = element_text(size = general_size))

top10class_plot

dev.off()
ggsave("~/Documents/plots_project2/top10classes.svg", top10class_plot, width = 110, height = 70, unit = "mm")
dev.off()  


unigenes_class_plot <- proportion_new_level_tool %>% mutate(gene = ifelse(new_level %in% top20, new_level, "Other")) %>% 
  mutate(gene = factor(gene, levels = top20)) %>% 
  ungroup() %>%
  group_by(tool, gene) %>%
  summarise(tot = sum(n)) %>%
  ggplot(aes( x = gene, y = tot, fill = tool, color = tool)) +
  geom_col(position = position_dodge2(preserve = "single"), color = "black", linewidth = 0.1) +
  scale_y_continuous(breaks = c(10000, 20000,30000,40000,50000)) +
  facet_grid(. ~ gene, scales  = "free_x") +
  scale_fill_manual(values = pal_10_q, labels = tool_label) +
  scale_color_manual(values = pal_10_q, labels = tool_label) +
  xlab("") +
  ggtitle("A") + 
  ylab("Genes") +
  theme_minimal() +
  labs(fill = "", color = "") +
  #scale_y_continuous(limits = c(-20, 55000), expand = c(0, 0)) +
  theme(legend.position = "none",
        panel.border = element_rect(fill = "transparent", color = "black", linewidth = 1),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = general_size),
        axis.text.y = element_text(size = general_size),
        #strip.text = element_text(size = general_size, face = "bold"),
        strip.text = element_blank(),
        plot.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.margin = margin(0, 0, 0, 0, unit = "pt"),
        panel.spacing = unit(0, "pt"),
        legend.text = element_text(size = general_size),
        title = element_text(size = general_size + 2, face = "bold"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title = element_text(size = general_size + 1, face = "bold"))

unigenes_class_plot


ggsave("~/Documents/plots_project2/unigenes_class.svg", unigenes_class_plot, width = 180, height = 70, unit = "mm")

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
  geom_boxplot(aes(y = total, fill = tool), outlier.shape = NA, linewidth = 0.2) +
  scale_fill_manual(values = pal_10_q, labels = tool_label) +
  theme_minimal() +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)),
                limits = c(1, 20000)) + 
  ylab("Abundance") +
  xlab("") +
  labs(fill = "") +
  ggtitle("A") +
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
    title = element_text(size = general_size + 2, face = "bold"),
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
  geom_boxplot(aes(y = total, fill = habitat2), outlier.shape = NA, linewidth = 0.2) +
  facet_grid( . ~ location) +
  scale_fill_manual(values = pal_10_q) +
  theme_minimal() +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)),
                limits = c(1, 50000)) + 
  ylab("Abundance") +
  xlab("") +
  labs(fill = "") +
  ggtitle("A") +
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
    title = element_text(size = general_size + 2, face = "bold"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = general_size),
    axis.text.y = element_text(size = general_size))

abundance_plot_inv + theme( legend.position = "none")

dev.off()
ggsave("~/Documents/plots_project2/abundance_inverted.svg", abundance_plot_inv + theme( legend.position = "none"), width = 180, height = 70, unit = "mm")
dev.off()  


dev.off()
ggsave("~/Documents/plots_project2/abundance_legend.svg", g_legend(abundance_plot_inv), width = 180, height = 10, unit = "mm")
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

top20 <- unique(c(top20, "Cell wall charge", "ERM", "Class C", "Target-modifying enzyme", "ABC-F"))
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


abundance_medians_human <- abundance_class %>% ungroup() %>%
  group_by(tool, gene) %>%
  summarise(md = median(normed10m), q75 = quantile(normed10m, 0.75), q25 = quantile(normed10m, 0.25)) 

abundance_medians_human_non_zero <- abundance_medians_human %>% filter(q75 > 0) 

abundance_class %>% ungroup() %>% 
  group_by(tool, gene) %>% filter(gene %in% "ERM", habitat %in% "human gut") %>%
  summarise(md = median(normed10m), q75 = quantile(normed10m, 0.75), q25 = quantile(normed10m, 0.25),
            sn = sum(normed10m), md_genes = median(unigenes), su= sum(unigenes), n = n_distinct(sample)) 

abundance_class %>% ungroup() %>% 
  group_by(tool, gene) %>% filter(gene %in% "ERM") %>%
  summarise(md = median(normed10m), q75 = quantile(normed10m, 0.75), q25 = quantile(normed10m, 0.25),
            md_genes = median(unigenes)) 

abundance_human_medians_heatmap <- abundance_medians_human %>% 
  ungroup() %>% 
  group_by(gene) %>% 
  mutate(mx =  max(md)) %>% 
  filter(mx > 0) %>% 
  mutate(norm  = log(md + 1) / log(mx + 1)) %>%
ggplot(aes(x = tool, y = gene, fill = norm)) +
  geom_tile() +
  scale_x_discrete(labels = tool_label) +
  scale_fill_viridis_c() + 
  theme_minimal() +
  labs(fill = "") +
  xlab("") +
  ylab("Gene class") +
  ggtitle("A") + 
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
    axis.title = element_text(size = general_size + 1, face = "bold"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = general_size),
    axis.text.y = element_text(size = general_size))

abundance_human_medians_heatmap


human_abundance_class_plot <- abundance_class  %>% 
  filter(habitat %in% "human gut") %>% 
  ungroup() %>%
  mutate(g2 = ifelse(gene %in% top_abundance, gene, "Other")) %>% 
  mutate(g2 = factor(g2, levels = top_abundance)) %>% 
  group_by(habitat, tool, g2, sample) %>% summarise(total = sum(normed10m)) %>%
  filter(paste(g2, tool) %in% paste(abundance_medians_human_non_zero$gene, abundance_medians_human_non_zero$tool)) %>%
  ggplot(aes( x = g2, y = total +1, fill = tool)) +
  geom_boxplot(position = position_dodge2(preserve = "single"),  outlier.shape = NA, coef = 0, width = 1, linewidth = 0.2) +
  scale_fill_manual(values = pal_10_q, labels = tool_label) +
  facet_grid(. ~ g2, scales = "free_x") +
  theme_minimal() +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) + #,
                #limits = c(1, 10000)) + 
  coord_cartesian(ylim = c(1, 5000)) +
  ylab("Abundance") +
  xlab("") +
  labs(fill = "") +
  ggtitle("B") +
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
    strip.text.x = element_blank(),
    title = element_text(size = general_size + 2, face = "bold"),
    strip.text.y = element_blank(),
    axis.title = element_text(size = general_size + 1, face = "bold"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = general_size),
    axis.text.y = element_text(size = general_size)) 


human_abundance_class_plot 

dev.off()
ggsave("~/Documents/plots_project2/abundance_class_humangut.svg", human_abundance_class_plot + theme(legend.position = "none"), width = 180, height = 70, unit = "mm")
dev.off()  


ggsave("~/Documents/plots_project2/abundance_class_humangut_heatmap.svg", abundance_human_medians_heatmap, width = 50, height = 70, unit = "mm")



abundance_medians_human %>% 
  ungroup() %>% 
  group_by(gene) %>% 
  mutate(mx =  max(md)) %>% filter(gene %in% "ERM")
  
top_abundance <- abundance_medians_human %>% 
  ungroup() %>% 
  group_by(gene) %>% 
  mutate(mx =  max(md)) %>% 
  filter(md > 0) %>% select(gene) %>% distinct() %>% pull()

top_abundance <- unique(c(as.character(top20[1:(length(top20)-1)]), top_abundance))

abundance_medians_human_jitter <- abundance_medians_human %>% 
  ungroup() %>% 
  group_by(gene) %>% 
  mutate(mx =  max(md)) %>% 
  filter(q75 > 0) %>% 
  mutate(gene = factor(gene, levels = top_abundance)) %>%
  ggplot(aes(x = gene, y = md + 1, fill = tool, color = tool, group = tool)) +
  geom_pointrange(aes(ymin = q25 + 1, ymax = q75 + 1), 
                  position=position_jitter(width = 0.5)) +
  theme_minimal() +
  labs(fill = "", color = "") +
  facet_grid( . ~ gene , scales = "free_x" ) +
  xlab("") +
  ylab("") +
  ggtitle("A") + 
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))  +
  scale_fill_manual(values = pal_10_q, labels = tool_label) +
  scale_color_manual(values = pal_10_q, labels = tool_label) +
  theme(legend.position = "bottom",
        panel.border = element_rect(fill = "transparent", color = "black", linewidth = 1),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = general_size),
        axis.text.y = element_text(size = general_size),
        strip.text = element_blank(),
        plot.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.margin = margin(0, 0, 0, 0, unit = "pt"),
        panel.spacing = unit(0, "pt"),
        legend.text = element_text(size = general_size),
        panel.grid.major.x = element_blank(),
        title = element_text(size = general_size + 2, face = "bold"),
        panel.grid.minor.x = element_blank(),
        axis.title = element_text(size = general_size + 1, face = "bold"))

abundance_medians_human_jitter

ggsave("~/Documents/plots_project2/abundance_medians_human_jitter.svg", abundance_medians_human_jitter + theme(legend.position = "none"), width = 180, height = 73, unit = "mm")

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
  geom_boxplot(aes(y = total, fill = habitat2), outlier.shape = NA, linewidth = 0.2) +
  scale_fill_manual(values = pal_10_q) +
  theme_minimal() +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)),
                limits = c(1, 5000)) + 
  ylab("Diversity") +
  xlab("") +
  labs(fill = "") +
  ggtitle("B") +
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
    title = element_text(size = general_size + 2, face = "bold"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = general_size),
    axis.text.y = element_text(size = general_size))


diversity_plot_inv

dev.off()  
ggsave("~/Documents/plots_project2/diversity_inverted.svg", diversity_plot_inv + theme(legend.position = "none"), width = 180, height = 70, unit = "mm")

abundance_plot_inv

### Correlations

cor_abu_div <- abundance_class %>% mutate(habitat2 = abundance$habitat2[match(habitat, abundance$habitat)]) %>%
  group_by(habitat, tool) %>% filter(normed10m != 0 & unigenes != 0) %>% summarise(corr = cor(normed10m, unigenes, method = "pearson"))
data.frame(cor_abu_div)


##############################################################
##############################################################

abu_tool_habitat <- abundance  %>% 
  group_by(habitat, tool, sample, aggregation, gene) %>% 
  summarise(abundance = sum(normed10m),
            diversity = sum(unigenes) ) 


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

JI_class_other_filter %>% filter(tool_ref %in% "fARGene (nt)") %>% ungroup() %>% group_by(new_level) %>% summarise(m = median(recall)) %>%
  mutate(M = median(m))

JI_class_other_filter %>% filter(tool_ref %in% "RGI (DIAMOND - nt)") %>% ungroup() %>% group_by(new_level) %>% summarise(m = median(recall)) %>%
  mutate(M = median(m))

JI_class_other_filter %>% filter(tool_ref %in% "DeepARG (nt)") %>% ungroup() %>% group_by(new_level) %>% summarise(m = median(recall)) %>%
  mutate(M = median(m, na.rm = T))

JI_class_other_filter %>% filter(tool_ref %in% "ABRicate (MEGARES - nt)") %>% ungroup() %>% group_by(new_level) %>% summarise(m = median(recall)) %>%
  mutate(M = mean(m))

JI_class_other_filter %>% filter(tool_ref %in% "AMRFinderPlus (aa)") %>% ungroup() %>% group_by(new_level) %>% summarise(m = median(recall)) %>%
  mutate(M = median(m))

JI_class_other_filter %>% filter(tool_ref %in% "ABRicate (MEGARES - nt)") %>% filter(grepl("ABRicate", tool_comp)) %>% ungroup() %>% group_by(new_level) %>% summarise(m = median(recall)) %>%
  mutate(M = mean(m), mn = min(m)) %>% ggplot(aes(x = new_level, y = m)) + geom_point()


JI_class_other_filter %>% filter(tool_ref %in% "ABRicate (MEGARES - nt)") %>% filter(grepl("ABRicate", tool_comp)) %>% ungroup() %>% group_by(new_level) %>% summarise(m = median(recall)) %>%
  mutate(M = mean(m), mn = min(m)) %>% filter(m >= .97)

JI_class_other_filter %>% filter(tool_ref %in% "ABRicate (MEGARES - nt)") %>% filter(!grepl("ABRicate", tool_comp)) %>% ungroup() %>% group_by(new_level) %>% summarise(m = median(recall)) %>%
  mutate(M = median(m), mn = min(m)) 


JI_class_other_filter %>% 
  ungroup() %>% 
  group_by(tool_ref, new_level) %>% 
  summarise(md1 = median(recall, na.rm = T)) %>% 
  ungroup() %>% 
  group_by(tool_ref) %>% 
  summarise(md2 = median(md1, na.rm = T), iqr = IQR(md1, 0.25))


recall_plot <- JI_class_other_filter %>% filter(!is.na(recall)) %>%
  ggplot(aes(x = tool_ref, y = recall, fill = tool_ref)) + 
  geom_boxplot( linewidth = 0.2) +
  scale_fill_manual(values = pal_10_q, labels = tool_label) +
  scale_x_discrete( labels = tool_label) +
  theme_minimal() +
  coord_cartesian(ylim = c(0,1)) +
  ylab("Recall") +
  xlab("") +
  labs(fill = "") +
  ggtitle("A") +
  theme(
    legend.position = "none",
    legend.text = element_text(size = general_size ),
    panel.border = element_rect(fill = "transparent", color = "black", linewidth = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    plot.margin = margin(0, 0, 0, 0, unit = "pt"),
    legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
    legend.margin = margin(0, 0, 0, 0, unit = "pt"),
    panel.spacing = unit(0, "pt"),
    title = element_text(size = general_size + 2, face = "bold"),
    axis.title = element_text(size = general_size + 1, face = "bold"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = general_size),
    axis.text.y = element_text(size = general_size))

recall_plot

fnr_plot <- JI_class_other_filter %>% filter(!is.na(fnr)) %>%
  ggplot(aes(x = tool_ref, y = fnr, fill = tool_ref)) + 
  geom_boxplot(linewidth = 0.2) +
  scale_fill_manual(values = pal_10_q, labels = tool_label) +
  scale_x_discrete( labels = tool_label) +
  theme_minimal() +
  coord_cartesian(ylim = c(0,1)) +
  ylab("False negative rate") +
  xlab("") +
  labs(fill = "") +
  ggtitle("B") +
  theme(
    legend.position = "none",
    legend.text = element_text(size = general_size ),
    panel.border = element_rect(fill = "transparent", color = "black", linewidth = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    plot.margin = margin(0, 0, 0, 0, unit = "pt"),
    legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
    legend.margin = margin(0, 0, 0, 0, unit = "pt"),
    panel.spacing = unit(0, "pt"),
    title = element_text(size = general_size + 2, face = "bold"),
    axis.title = element_text(size = general_size + 1, face = "bold"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = general_size),
    axis.text.y = element_text(size = general_size))

fnr_plot


ggsave("~/Documents/plots_project2/recall.svg", recall_plot , width = 85, height = 70, unit = "mm")
ggsave("~/Documents/plots_project2/fdr.svg", fnr_plot   , width = 85, height = 70, unit = "mm")

##

JI_all %>% filter(as.numeric(tool_ref) >= as.numeric(tool_comp)) %>% 
  ggplot(aes(x = tool_comp, y = tool_ref, fill = jaccard)) +
  geom_tile() +
  #scale_fill_gradient2(low = "white",  high = "red") +
  scale_fill_viridis_c() + 
  theme_minimal() +
  ggtitle("B") + 
  labs(fill = "") +
  xlab("") +
  theme(axis.title.y = element_blank(),
        legend.position = "bottom",
        panel.border = element_rect(fill = "transparent", color = "black", linewidth = 1),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = general_size),
        axis.text.y = element_text(size = general_size),
        strip.text = element_text(size = general_size, face = "bold"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        title = element_text(size = general_size + 2, face = "bold"),
        axis.title = element_text(size = general_size + 1, face = "bold"))


JI_all %>% filter(as.numeric(tool_ref) >= as.numeric(tool_comp)) %>% filter(tool_ref != tool_comp) %>%
  summarise(mn = mean(jaccard))
  
JI_all %>% filter(tool_ref != tool_comp) %>% filter(tool_ref %in% "DeepARG (nt)", tool_comp %in% c("ABRicate (CARD - nt)",
                                                                                                   "ABRicate (NCBI - nt)",
                                                                                                   "ABRicate (ResFinder - nt)",
                                                                                                   "ABRicate (ARG-ANNOT - nt)",
                                                                                                   "ABRicate (MEGARES - nt)",
                                                                                                   "ResFinder (nt)")) %>% 
  summarise(mn = mean(recall))

JI_all %>% filter(tool_ref != tool_comp) %>% filter(tool_ref %in% "RGI (DIAMOND - nt)", tool_comp %in% c("ABRicate (CARD - nt)",
                                                                                                   "ABRicate (NCBI - nt)",
                                                                                                   "ABRicate (ResFinder - nt)",
                                                                                                   "ABRicate (ARG-ANNOT - nt)",
                                                                                                   "ABRicate (MEGARES - nt)",
                                                                                                   "ResFinder (nt)")) %>% 
  summarise(mn = mean(recall))


JI_all %>% filter(tool_ref != tool_comp) %>% filter(tool_ref %in% "ABRicate (MEGARES - nt)", tool_comp %in% c("ABRicate (CARD - nt)",
                                                                                                         "ABRicate (NCBI - nt)",
                                                                                                         "ABRicate (ResFinder - nt)",
                                                                                                         "ABRicate (ARG-ANNOT - nt)",
                                                                                                         "ABRicate (MEGARES - nt)",
                                                                                                         "ResFinder (nt)")) %>% 
  summarise(mn = mean(recall))




heatmap_proportion_unigenes_found_in_other <- JI_all %>% #filter(as.numeric(tool_ref) >= as.numeric(tool_comp)) %>% 
  ggplot(aes(x = tool_comp, y = tool_ref, fill = recall)) +
  geom_tile() +
  scale_fill_viridis_c() + 
  scale_x_discrete(labels = tool_label) +
  scale_y_discrete(labels = tool_label) +
  theme_minimal() +
  ggtitle("B") + 
  labs(fill = "") +
  xlab("") +
  theme(axis.title.y = element_blank(),
        legend.position = "right",
        panel.border = element_rect(fill = "transparent", color = "black", linewidth = 1),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = general_size),
        axis.text.y = element_text(size = general_size),
        strip.text = element_text(size = general_size, face = "bold"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        legend.text = element_text(size = general_size ),
        plot.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.margin = margin(0, 0, 0, 0, unit = "pt"),
        title = element_text(size = general_size + 2, face = "bold"),
        axis.title = element_text(size = general_size + 1, face = "bold"))

ggsave("~/Documents/plots_project2/proportion_found_in_other_tool.svg", heatmap_proportion_unigenes_found_in_other  , width = 90, height = 70, unit = "mm")

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


id3 <- unigenes %>% 
  filter(tool %in% c("DeepARG (nt)", "RGI (DIAMOND - nt)")) %>% 
  ungroup() %>% 
  group_by(tool, new_level) %>% 
  summarise(n = n()) %>% 
  arrange(desc(n)) %>% 
  slice_head(n = 6) %>% 
  ungroup() %>% 
  select(new_level) %>% distinct() %>% pull()

id3 <- unigenes %>% 
  filter(tool %in% c("DeepARG (nt)", "RGI (DIAMOND - nt)")) %>% 
  ungroup() %>% 
  group_by(tool, new_level) %>% 
  summarise(n = n()) %>% 
  arrange(desc(n)) %>% 
  slice_head(n = 20) %>% 
  ungroup() %>% 
  select(new_level) %>% distinct() %>% pull()

data_plot_rgi_deep <- bind_rows(
  unigenes %>% #filter(tool %in% c("DeepARG (nt)", "RGI (DIAMOND - nt)")) %>%
    mutate(id = id/100, dt = "unigenes") %>%
    ungroup() %>% group_by(tool, new_level) %>% 
    summarise(md = median(id, na.rm = T), q25 = quantile(id , 0.25, na.rm = T), q75 = quantile(id , 0.75, na.rm = T)) %>% 
    mutate(dt = "Identity"),
  JI_class_other_filter %>% #filter(tool_ref %in% c("DeepARG (nt)", "RGI (DIAMOND - nt)")) %>% 
  select(tool_ref, recall, new_level) %>% 
  ungroup() %>% group_by(tool_ref, new_level) %>% 
  summarise(md = median(recall, na.rm = T), q25 = quantile(recall , 0.25, na.rm = T), q75 = quantile(recall , 0.75, na.rm = T)) %>% 
  rename(tool = tool_ref) %>%
  mutate(dt  = "Recall"),
  JI_class_other_filter %>% #filter(tool_ref %in% c("DeepARG (nt)", "RGI (DIAMOND - nt)")) %>% 
    select(tool_ref, fnr, new_level) %>% 
    ungroup() %>% group_by(tool_ref, new_level) %>% 
    summarise(md = median(fnr, na.rm = T), q25 = quantile(fnr , 0.25, na.rm = T), q75 = quantile(fnr , 0.75, na.rm = T)) %>% 
    rename(tool = tool_ref) %>%
    mutate(dt  = "FNR"))

data_plot_rgi_deep %>% filter(new_level %in% id3) %>% 
  ggplot(aes(x = md, y = new_level, fill = tool, alpha = dt, color = tool)) +
  geom_point() + 
  facet_grid( . ~ tool)




rgi_id_recall_fnr <- data_plot_rgi_deep  %>% filter(tool %in% c("RGI (DIAMOND - nt)")) %>% 
  mutate(dt = factor(dt, levels = c("Identity", "Recall", "FNR"))) %>% 
  ungroup() %>% group_by(dt, tool) %>%
  arrange(dt, md) %>% 
  mutate(cl1 = factor(new_level, levels = unique(new_level))) %>%
  mutate(cl = as.numeric(cl1)) %>%
  mutate(cl = ifelse(dt %in% c("Identity"), cl - 0.15, cl)) %>%
  mutate(cl = ifelse(dt %in% c("FNR"), cl + 0.15, cl))


deep_id_recall_fnr <- data_plot_rgi_deep  %>% filter(tool %in% c("DeepARG (nt)")) %>% 
  mutate(dt = factor(dt, levels = c("Identity", "Recall", "FNR"))) %>% 
  ungroup() %>% group_by(dt, tool) %>%
  arrange(dt, md) %>% 
  mutate(cl1 = factor(new_level, levels = unique(new_level))) %>%
  mutate(cl = as.numeric(cl1)) %>%
  mutate(cl = ifelse(dt %in% c("Identity"), cl - 0.15, cl)) %>%
  mutate(cl = ifelse(dt %in% c("FNR"), cl + 0.15, cl))

fargene_id_recall_fnr <- data_plot_rgi_deep  %>% filter(tool %in% c("fARGene (nt)")) %>% 
  mutate(dt = factor(dt, levels = c("Identity", "Recall", "FNR"))) %>% 
  ungroup() %>% group_by(dt, tool) %>%
  arrange(dt, md) %>% 
  mutate(cl1 = factor(new_level, levels = unique(new_level))) %>%
  mutate(cl = as.numeric(cl1)) %>%
  mutate(cl = ifelse(dt %in% c("Identity"), cl - 0.15, cl)) %>%
  mutate(cl = ifelse(dt %in% c("FNR"), cl + 0.15, cl))


fargene_id_recall_fnr_plot <- fargene_id_recall_fnr %>% 
  ggplot(aes( y = cl, x = md, fill = dt, color = dt)) +
  geom_rect(aes(
    ymin = as.numeric(cl) - 0.06,
    ymax = as.numeric(cl) + 0.06,
    xmin = q25,
    xmax = q75), color = "black", linewidth = 0.2) + 
  geom_point(shape = 15, aes(color = dt), alpha = 0.75, size = 3) +
  scale_color_manual(values = pal_10_q) +
  scale_fill_manual(values = pal_10_q) +
  theme_minimal() +
  scale_y_continuous(limits = c(0, length(levels(fargene_id_recall_fnr$cl1)) + 0.3), 
                     expand = c(0, 0), 
                     breaks = 1:length(levels(fargene_id_recall_fnr$cl1)), 
                     labels = levels(fargene_id_recall_fnr$cl1)) + 
  facet_grid(. ~ tool, scales = "free_y") +
  ylab("") +
  xlab("") +
  labs(fill = "", color = "") +
  ggtitle("B") +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = general_size ),
    panel.border = element_rect(fill = "transparent", color = "black", linewidth = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    plot.margin = margin(0, 0, 0, 0, unit = "pt"),
    strip.text = element_text(size = general_size , face = "bold"),
    legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
    legend.margin = margin(0, 0, 0, 0, unit = "pt"),
    panel.spacing = unit(0, "pt"),
    title = element_text(size = general_size + 2, face = "bold"),
    axis.title = element_text(size = general_size + 1, face = "bold"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = general_size),
    axis.text.y = element_text(size = general_size))

fargene_id_recall_fnr_plot


deep_id_recall_fnr_plot <- deep_id_recall_fnr %>% 
  ggplot(aes( y = cl, x = md, fill = dt, color = dt)) +
  geom_rect(aes(
    ymin = as.numeric(cl) - 0.09,
    ymax = as.numeric(cl) + 0.09,
    xmin = q25,
    xmax = q75), color = "black", linewidth = 0.1) + 
  geom_point(shape = 15, aes(color = dt), alpha = 0.75) +
  scale_color_manual(values = pal_10_q) +
  scale_fill_manual(values = pal_10_q) +
  theme_minimal() +
  scale_y_continuous(limits = c(0, length(levels(deep_id_recall_fnr$cl1)) + 0.4), 
                     expand = c(0, 0), 
                     breaks = 1:length(levels(deep_id_recall_fnr$cl1)), 
                     labels = levels(deep_id_recall_fnr$cl1)) + 
  facet_grid(. ~ tool, scales = "free_y") +
  ylab("") +
  xlab("") +
  labs(fill = "", color = "") +
  ggtitle("B") +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = general_size ),
    panel.border = element_rect(fill = "transparent", color = "black", linewidth = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    plot.margin = margin(0, 0, 0, 0, unit = "pt"),
    strip.text = element_text(size = general_size , face = "bold"),
    legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
    legend.margin = margin(0, 0, 0, 0, unit = "pt"),
    panel.spacing = unit(0, "pt"),
    title = element_text(size = general_size + 2, face = "bold"),
    axis.title = element_text(size = general_size + 1, face = "bold"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = general_size),
    axis.text.y = element_text(size = general_size))


rgi_id_recall_fnr_plot <- rgi_id_recall_fnr %>% 
  ggplot(aes( y = cl, x = md, fill = dt, color = dt)) +
  geom_rect(aes(
    ymin = as.numeric(cl) - 0.09,
    ymax = as.numeric(cl) + 0.09,
    xmin = q25,
    xmax = q75), color = "black", linewidth = 0.1) + 
  geom_point(shape = 15, aes(color = dt), alpha = 0.75) +
  scale_color_manual(values = pal_10_q) +
  scale_fill_manual(values = pal_10_q) +
  theme_minimal() +
  scale_y_continuous(limits = c(0, length(levels(rgi_id_recall_fnr$cl1)) + 0.4), 
                     expand = c(0, 0), 
                     breaks = 1:length(levels(rgi_id_recall_fnr$cl1)), 
                     labels = levels(rgi_id_recall_fnr$cl1)) + 
  facet_grid(. ~ tool, scales = "free_y") +
  ylab("") +
  xlab("") +
  labs(fill = "", color = "") +
  ggtitle("C") +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = general_size ),
    panel.border = element_rect(fill = "transparent", color = "black", linewidth = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    plot.margin = margin(0, 0, 0, 0, unit = "pt"),
    strip.text = element_text(size = general_size , face = "bold"),
    legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
    legend.margin = margin(0, 0, 0, 0, unit = "pt"),
    panel.spacing = unit(0, "pt"),
    title = element_text(size = general_size + 2, face = "bold"),
    axis.title = element_text(size = general_size + 1, face = "bold"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = general_size),
    axis.text.y = element_text(size = general_size))

rgi_id_recall_fnr_plot
deep_id_recall_fnr_plot
fargene_id_recall_fnr_plot

ggsave("~/Documents/plots_project2/deep_id_recall.svg", deep_id_recall_fnr_plot + theme( legend.position = "none"), width = 87, height = 135, unit = "mm")
ggsave("~/Documents/plots_project2/rgi_id_recall.svg", rgi_id_recall_fnr_plot + theme( legend.position = "none"), width = 87, height = 135, unit = "mm")
ggsave("~/Documents/plots_project2/legend_id_recall.svg", g_legend(rgi_id_recall_fnr_plot), width = 180, height = 10, unit = "mm")

JI_class_other_filter %>% filter(tool_ref %in% c("DeepARG (nt)", "RGI (DIAMOND - nt)")) %>% 
  select(tool_ref, recall, new_level) %>% 
  mutate(tool = tool_ref, dt  = "recall")
  
JI_class_other_filter %>% filter(tool_ref %in% c("DeepARG (nt)", "RGI (DIAMOND - nt)")) %>% 
  ggplot(aes(x = recall, y = new_level, fill = tool_ref)) +
  geom_boxplot(outlier.shape = NA, outlier.size = NA)  + 
  facet_grid( . ~ tool_ref)


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


id_plot <- unigenes %>% ungroup() %>% filter(tool %in% tool_2) %>% 
  mutate(t = ifelse(tool %in% c("fARGene (nt)", "AMRFinderPlus (aa)", "RGI (DIAMOND - nt)", "DeepARG (nt)", "ResFinder (nt)"), "n", "y")) %>%
  ggplot(aes(x = id, color = tool, y=..count../sum(count))) +
  geom_freqpoly(bins = 35, linewidth = .7) +
  facet_grid(. ~ t ) +
  theme_minimal() + 
  scale_fill_manual(values = pal_10_q, labels = tool_label) +
  scale_color_manual(values = pal_10_q, labels = tool_label) +
  ggtitle("A") + 
  ylab("Proportion of genes") + 
  xlab("Identity") + 
  theme(legend.position = "none",
        panel.border = element_rect(fill = "transparent", color = "black", linewidth = 1),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = general_size),
        axis.text.y = element_text(size = general_size),
        strip.text = element_blank(),
        plot.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.margin = margin(0, 0, 0, 0, unit = "pt"),
        panel.spacing = unit(0, "pt"),
        legend.text = element_text(size = general_size),
        title = element_text(size = general_size + 2, face = "bold"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title = element_text(size = general_size + 1, face = "bold"))

ggsave("~/Documents/plots_project2/id_tool.svg", id_plot + theme( legend.position = "none"), width = 180, height = 60, unit = "mm")

  
unigenes %>% ungroup() %>% filter(tool %in% tool_2) %>% 
  mutate(t = ifelse(tool %in% c("fARGene (nt)", "AMRFinderPlus (aa)", "RGI (DIAMOND - nt)", "DeepARG (nt)"), "n", "y")) %>%
  ggplot(aes(x = id, color = tool)) +
  stat_ecdf() +
  facet_grid(. ~ t )

##################################################################################################################################################################
##################################################################################################################################################################
##################################################################################################################################################################
##################################################################################################################################################################
##################################################################################################################################################################
##################################################################################################################################################################

#core <- readRDS(file = "code_R_analysis/output_abundance_diversity_resistome/core_resistome.rds")
#pan <- readRDS(file = "code_R_analysis/output_abundance_diversity_resistome/pan_resistome.rds")

core <- readRDS(file = "code_R_analysis/output_abundance_diversity_resistome/core_resistome_no_raw_unique_filter.rds")
pan <- readRDS(file = "code_R_analysis/output_abundance_diversity_resistome/pan_resistome_no_rawunique_filter.rds")

core <- core %>% 
  mutate(tool = ifelse(tool == "RGI (DIAMOND nt)", "RGI (DIAMOND - nt)", tool)) %>%
  filter(tool %in% tool_2, !habitat %in% not_env) %>% 
  mutate(habitat = factor(habitat, levels = EN), tool = factor(tool, levels =  tools_levels))

pan <- pan %>% 
  mutate(tool = ifelse(tool == "RGI (DIAMOND nt)", "RGI (DIAMOND - nt)", tool)) %>% 
  filter(tool %in% tool_2, !habitat %in% not_env, aggregation %in% "new_level") %>% 
  mutate(habitat = factor(habitat, levels = EN), tool = factor(tool, levels =  tools_levels))

sumcore <- core %>% 
  filter(cut %in% 0.5 & cnt > 900) %>% 
  ungroup() %>% 
  group_by(new_level, tool, habitat) %>% 
  summarise(unigenes = n_distinct(X))  

sumpan <- pan %>% ungroup() %>% group_by(tool, habitat, aggregation, gene_class) %>% 
  summarise(md = median(unigenes), mn = mean(unigenes)) 

top20_c <- sumcore %>% filter(habitat %in% c("human gut")) %>% ungroup() %>% 
  arrange(desc(unigenes)) %>% group_by(tool) %>% slice_head(n = 6) %>% ungroup() %>% select(new_level) %>% distinct() %>% pull()

top20_c %in% top20
top20 %in% top20_c

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
  ylab("Core-resistome") +
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
  ylab("ARGs per class") +
  scale_fill_manual(values = c(pal_10_q), labels = tool_label) +
  theme(legend.position = "bottom",
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


p2.2 <- human_resistome %>% filter(d %in% "Core-resistome") %>% 
  group_by(tool, habitat, new_level) %>% 
  summarise(unigenes = sum(unigenes), d = d[1], pattern = pattern[1]) %>%
  ggplot(aes(x = unigenes + 1, y = new_level, color = tool)) +
  geom_jitter(height = 0.6, size = 2) +
  theme_minimal() +
  labs(color = "") +
  facet_grid(new_level ~ ., scales = "free_y") + 
  xlab("") +
  ylab("") +
  ggtitle("B") +
  scale_color_manual(values = c(pal_10_q), labels = tool_label) +
  scale_x_continuous(limits = c(0, 200), expand = c(0, 0)) +
  #scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
  #              labels = scales::trans_format("log10", scales::math_format(10^.x)),
  #              limits = c(1, 200)) + 
  theme(legend.position = "bottom",
        panel.border = element_rect(fill = "transparent", color = "black", linewidth = 1),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = general_size),
        axis.text.y = element_text(size = general_size),
        #strip.text = element_text(size = general_size, face = "bold"),
        strip.text = element_blank(),
        plot.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.margin = margin(0, 0, 0, 0, unit = "pt"),
        panel.spacing = unit(0, "pt"),
        legend.text = element_text(size = general_size),
        title = element_text(size = general_size + 2, face = "bold"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title = element_text(size = general_size + 1, face = "bold"))

p2.2


p2.21 <- human_resistome %>% filter(d %in% "Core-resistome") %>% 
  group_by(tool, habitat, new_level) %>% 
  summarise(unigenes = sum(unigenes), d = d[1], pattern = pattern[1]) %>%
  ggplot(aes(y = unigenes, x = new_level, fill = tool)) +
  #geom_jitter(height = 0.6, size = 2) +
  geom_col(position = position_dodge2(preserve = "single"), color = "black", linewidth = 0.1) +
  theme_minimal() +
  labs(color = "") +
  facet_grid(. ~ new_level, scales = "free_x") + 
  xlab("") +
  ylab("") +
  ggtitle("C") +
  scale_fill_manual(values = c(pal_10_q), labels = tool_label) +
  scale_y_continuous(limits = c(-2, 175), expand = c(0, 0)) +
  theme(legend.position = "none",
        panel.border = element_rect(fill = "transparent", color = "black", linewidth = 1),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = general_size),
        axis.text.y = element_text(size = general_size),
        #strip.text = element_text(size = general_size, face = "bold"),
        strip.text = element_blank(),
        plot.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.margin = margin(0, 0, 0, 0, unit = "pt"),
        panel.spacing = unit(0, "pt"),
        legend.text = element_text(size = general_size),
        title = element_text(size = general_size + 2, face = "bold"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title = element_text(size = general_size + 1, face = "bold"))

p2.21

p3_core <- p2 +
  annotation_custom(
    grob = ggplotGrob(p1),
    xmin = 5, xmax = 16.5,   # positioning within the data coordinates of p_main
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
ggsave("~/Documents/plots_project2/core_human_gut.svg", p3_core, width = 130, height = 120, unit = "mm")  
ggsave("~/Documents/plots_project2/core_class_human_gut.svg", p2.21, width = 180, height = 70, unit = "mm")  

ggsave("~/Documents/plots_project2/core_human_gut_no_raw_unique_filter.svg", p3_core, width = 130, height = 120, unit = "mm")  



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

p2.3 <- human_resistome %>% filter(d %in% "Pan-resistome") %>% 
  group_by(tool, habitat, new_level) %>% 
  summarise(unigenes = sum(unigenes), d = d[1], pattern = pattern[1]) %>%
  ggplot(aes(x = unigenes + 1, y = new_level, col = tool)) +
  geom_jitter(height = 0.6, size = 2) +
  theme_minimal() +
  labs(color = "") +
  facet_grid(new_level ~ ., scales = "free_y") + 
  xlab("") +
  ylab("") +
  ggtitle("A") +
  scale_color_manual(values = c(pal_10_q), labels = tool_label) +
  #scale_x_continuous(limits = c(0, 200), expand = c(0, 0)) +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)),
                limits = c(1, 9000)) + 
  theme(legend.position = "bottom",
        panel.border = element_rect(fill = "transparent", color = "black", linewidth = 1),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = general_size),
        axis.text.y = element_text(size = general_size),
        #strip.text = element_text(size = general_size, face = "bold"),
        strip.text = element_blank(),
        plot.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.margin = margin(0, 0, 0, 0, unit = "pt"),
        panel.spacing = unit(0, "pt"),
        legend.text = element_text(size = general_size),
        title = element_text(size = general_size + 2, face = "bold"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title = element_text(size = general_size + 1, face = "bold"))

p2.3

p3_pan <- p2 +
  annotation_custom(
    grob = ggplotGrob(p1),
    xmin = 5, xmax = 18,   # positioning within the data coordinates of p_main
    ymin = 3000, ymax = 8000
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

ggsave("~/Documents/plots_project2/pan_human_gut_no_raw_unique_filter.svg", p3_pan, width = 130, height = 120, unit = "mm") 


plot_env <- c("human gut", "human skin", "pig gut", "wastewater", "marine", "freshwater", "soil")


core_environment <- sumcore  %>% ungroup() %>% 
  group_by(habitat, tool) %>% 
  summarise(unigenes = sum(unigenes)) %>% 
  filter(habitat %in% plot_env) %>% 
  ggplot(aes(x = unigenes, y = habitat, fill = tool)) +
  geom_col(position = position_dodge2(preserve = "single", width = 1), width = 1) +
  theme_minimal() +
  labs(fill = "") +
  facet_grid(habitat ~ ., scales = "free_y") +
  scale_x_continuous(limits = c(-10, 1350), expand = c(0, 0)) +
  ylab("") +
  xlab("") +
  ggtitle("B") +
  scale_fill_manual(values = pal_10_q, labels = tool_label) +
  theme(legend.position = "right",
        legend.text = element_text(size = general_size),
        panel.border = element_rect(fill = "transparent", color = "black", linewidth = 1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        plot.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.margin = margin(0, 0, 0, 0, unit = "pt"),
        panel.spacing = unit(0, "pt"),
        title = element_text(size = general_size + 2, face = "bold"),
        strip.text = element_blank(),
        axis.title = element_text(size = general_size + 1, face = "bold"),
        axis.text.x = element_text(angle = 90, size = general_size),
        axis.text.y = element_text(size = general_size)) 

core_environment

ggsave("~/Documents/plots_project2/core_environment_no_raw_unique_filter.svg", core_environment + theme(legend.position = "none"), width = 90, height = 140, unit = "mm") 


core_environment_supplement <- sumcore  %>% ungroup() %>% 
  group_by(habitat, tool) %>% 
  summarise(unigenes = sum(unigenes)) %>% 
  filter(!habitat %in% plot_env) %>% 
  ggplot(aes(x = unigenes, y = habitat, fill = tool)) +
  geom_col(position = position_dodge2(preserve = "single", width = 1), width = 1) +
  theme_minimal() +
  labs(fill = "") +
  facet_grid(habitat ~ ., scales = "free_y") +
  scale_x_continuous(limits = c(-10, 1350), expand = c(0, 0)) +
  ylab("") +
  xlab("") +
  ggtitle("B") +
  scale_fill_manual(values = pal_10_q, labels = tool_label) +
  theme(legend.position = "right",
        legend.text = element_text(size = general_size),
        panel.border = element_rect(fill = "transparent", color = "black", linewidth = 1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        plot.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.margin = margin(0, 0, 0, 0, unit = "pt"),
        panel.spacing = unit(0, "pt"),
        title = element_text(size = general_size + 2, face = "bold"),
        strip.text = element_blank(),
        axis.title = element_text(size = general_size + 1, face = "bold"),
        axis.text.x = element_text(angle = 90, size = general_size),
        axis.text.y = element_text(size = general_size)) 

core_environment_supplement

ggsave("~/Documents/plots_project2/core_environment_no_raw_unique_filter.svg", core_environment + theme(legend.position = "none"), width = 90, height = 140, unit = "mm") 
ggsave("~/Documents/plots_project2/sup_core_environment_no_raw_unique_filter.svg", core_environment + theme(legend.position = "none"), width = 90, height = 140, unit = "mm") 

plot_env
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



################


sumpan2 <- pan %>% ungroup() %>% group_by(tool, habitat, aggregation, epoch) %>% 
  summarise(s = sum(unigenes)) %>%
  ungroup() %>% group_by(tool, habitat, aggregation) %>% 
  summarise(md = median(s), mn = mean(s), sd = sd(s)) 


#sumpan %>% filter(aggregation %in% "new_level", habitat %in% c("human gut", "wastewater", "pig gut", "marine", "soil", "built-environment")) %>% 

plot_env <- c("human gut", "human skin", "pig gut", "wastewater", "marine", "freshwater", "soil")

pan_resistome_plot_environment <- sumpan2 %>% 
  filter(habitat %in% plot_env) %>% 
  filter(aggregation %in% "new_level") %>% #, habitat %in% c("human gut", "wastewater", "pig gut", "marine", "soil", "built-environment")) %>% 
  ggplot(aes(x = mn + 1, y = habitat, fill = tool)) +
  geom_col(position = position_dodge2(preserve = "single", width = 1), width = 1) +
  geom_errorbar(aes(xmin = mn-sd, xmax = mn+sd), width=.2,
                position = position_dodge(.9)) +
  theme_minimal() +
  labs(fill = "") +
  facet_grid(habitat ~ ., scales = "free_y") + 
  ylab("") +
  xlab("") +
  ggtitle("A") + 
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) + 
  coord_cartesian(xlim = c(100,24000)) + 
  scale_fill_manual(values = pal_10_q, labels = tool_label) +
  theme(legend.position = "bottom",
        panel.border = element_rect(fill = "transparent", color = "black", linewidth = 1),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = general_size),
        axis.text.y = element_text(size = general_size),
        strip.text = element_blank(),
        plot.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.margin = margin(0, 0, 0, 0, unit = "pt"),
        title = element_text(size = general_size + 2, face = "bold"),
        panel.spacing = unit(0, "pt"),
        legend.text = element_text(size = general_size),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title = element_text(size = general_size + 1, face = "bold"))

pan_resistome_plot_environment 

ggsave("~/Documents/plots_project2/pan_environment_no_raw_unique_filter.svg", pan_resistome_plot_environment + theme(legend.position = "none"), width = 90, height = 140, unit = "mm") 
ggsave("~/Documents/plots_project2/pan_legend.svg", g_legend(pan_resistome_plot_environment), width = 180, height = 10, unit = "mm")


pan_resistome_plot_environment_supplement <- sumpan2 %>% 
  filter(!habitat %in% plot_env) %>% 
  filter(aggregation %in% "new_level") %>% #, habitat %in% c("human gut", "wastewater", "pig gut", "marine", "soil", "built-environment")) %>% 
  ggplot(aes(x = mn + 1, y = habitat, fill = tool)) +
  geom_col(position = position_dodge2(preserve = "single", width = 1), width = 1) +
  geom_errorbar(aes(xmin = mn-sd, xmax = mn+sd), width=.2,
                position = position_dodge(.9)) +
  theme_minimal() +
  labs(fill = "") +
  facet_grid(habitat ~ ., scales = "free_y") + 
  ylab("") +
  xlab("") +
  ggtitle("A") + 
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) + 
  coord_cartesian(xlim = c(100,5000)) + 
  scale_fill_manual(values = pal_10_q, labels = tool_label) +
  theme(legend.position = "bottom",
        panel.border = element_rect(fill = "transparent", color = "black", linewidth = 1),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = general_size),
        axis.text.y = element_text(size = general_size),
        strip.text = element_blank(),
        plot.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.margin = margin(0, 0, 0, 0, unit = "pt"),
        title = element_text(size = general_size + 2, face = "bold"),
        panel.spacing = unit(0, "pt"),
        legend.text = element_text(size = general_size),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title = element_text(size = general_size + 1, face = "bold"))

pan_resistome_plot_environment_supplement


ggsave("~/Documents/plots_project2/sup_pan_environment_no_raw_unique_filter.svg", pan_resistome_plot_environment_supplement + theme(legend.position = "none"), width = 90, height = 140, unit = "mm") 

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





# core_resistome_plot_environment_jitter <- sumcore2 %>% #, habitat %in% c("human gut", "wastewater", "pig gut", "marine", "soil", "built-environment")) %>% 
#   ggplot(aes(x = unigenes, y = habitat, color = tool)) +
#   geom_jitter(width = 0.2, size = 2) +
#   labs(fill = "", color = "") +
#   facet_grid( habitat ~ . , scales = "free_y" ) +
#   theme_minimal() +
#   xlab("") +
#   ylab("") +
#   ggtitle("B") + 
#   scale_color_manual(values = pal_10_q, label = tool_label) +
#   scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
#                 labels = scales::trans_format("log10", scales::math_format(10^.x)))  +
#   theme(legend.position = "bottom",
#         panel.border = element_rect(fill = "transparent", color = "black", linewidth = 1),
#         axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = general_size),
#         axis.text.y = element_text(size = general_size),
#         strip.text = element_blank(),
#         plot.margin = margin(0, 0, 0, 0, unit = "pt"),
#         legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
#         legend.margin = margin(0, 0, 0, 0, unit = "pt"),
#         panel.spacing = unit(0, "pt"),
#         legend.text = element_text(size = general_size),
#         panel.grid.major.x = element_blank(),
#         panel.grid.minor.x = element_blank(),
#         axis.title = element_text(size = general_size + 1, face = "bold"))
# 
# core_resistome_plot_environment_jitter
# 
# ggsave("~/Documents/plots_project2/core_environment.svg", core_resistome_plot_environment, width = 300, height = 100, unit = "mm")
# ggsave("~/Documents/plots_project2/pan_environment.svg", pan_resistome_plot_environment, width = 300, height = 100, unit = "mm")
# 
# ggsave("~/Documents/plots_project2/core_environment_jitter.svg", core_resistome_plot_environment_jitter, width = 80, height = 120, unit = "mm")
# ggsave("~/Documents/plots_project2/pan_environment_jitter.svg", pan_resistome_plot_environment_jitter, width = 80, height = 120, unit = "mm")



sumcore2 %>% ungroup() %>% arrange(habitat) %>%  
  group_by(habitat) %>%  
  mutate(deeparg = unigenes[tool == "DeepARG (nt)"][1],
         rgi = unigenes[tool == "RGI (DIAMOND - nt)"][1]) %>%
  mutate(dp = deeparg / unigenes, drgi = rgi / unigenes) %>% 
  filter(!tool %in% c("DeepARG (nt)", "RGI (DIAMOND - nt)")) %>% 
  summarise(m_da = mean(dp), m_rgi = mean(drgi))


sumcore_human_class <- sumcore %>% filter(habitat %in% "human gut") %>% 
  ungroup() %>% group_by(tool) %>% mutate(total = sum(unigenes)) %>% 
  mutate(p = unigenes / total )

data.frame(sumcore_human_class %>% 
  arrange(tool, desc(unigenes)) %>%  ungroup() %>% group_by(tool) %>% 
  slice_head(n = 5))

sumcore2 %>% filter(habitat %in% "human gut")


sumcore_human_class %>% filter(new_level %in% "rpoB") %>% ungroup() %>%
  mutate(deeparg = unigenes[tool == "DeepARG (nt)"][1],
         rgi = unigenes[tool == "RGI (DIAMOND - nt)"][1]) %>%
  mutate(dp = deeparg / unigenes, drgi = rgi / unigenes)


sumcore_human_class %>% filter(new_level %in% "GPA") %>% ungroup() %>%
  mutate(deeparg = unigenes[tool == "DeepARG (nt)"][1],
         rgi = unigenes[tool == "RGI (DIAMOND - nt)"][1]) %>%
  mutate(dp = deeparg / unigenes, drgi = rgi / unigenes)

sumcore_human_class %>% filter(new_level %in% "Efflux p.") %>% ungroup() %>%
  mutate(deeparg = unigenes[tool == "DeepARG (nt)"][1],
         rgi = unigenes[tool == "RGI (DIAMOND - nt)"][1]) %>%
  mutate(dp = deeparg / unigenes, drgi = rgi / unigenes)


sumcore_human_class %>% filter(new_level %in% "TET - RPG") %>% ungroup() %>%
  mutate(deeparg = unigenes[tool == "DeepARG (nt)"][1],
         rgi = unigenes[tool == "RGI (DIAMOND - nt)"][1]) %>%
  mutate(dp = deeparg / unigenes, drgi = rgi / unigenes) %>%
  summarise(m = mean(drgi))

#

sets0 <- core %>% filter(cut %in% 0.5 & cnt > 900, habitat %in% c( "human gut" )) %>%
  group_by(tool) %>%
  summarise(query = list(X), .groups = "drop")

sets1 <- core %>% filter(cut %in% 0.5 & cnt > 900, habitat %in% c( "human gut" )) %>%
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

new_union <- function(qc_ref, q_ref, qc_comp, q_comp){
  A <-  unlist(qc_ref)
  B <- unlist(q_ref)
  C <- unlist(qc_comp)
  D <- unlist(q_comp)
  x <- intersect(A, setdiff(D, C))
  y <- setdiff(C, intersect(C, setdiff(B, A)))
  z <- union(x, y)
  d <- union(A, z)
  r <- ifelse(length(d) == 0, 0, length(d))
  return(r)
}

new_intersect2 <- function(qc_ref, q_ref, qc_comp, q_comp){
  A <-  unlist(qc_ref)
  B <- unlist(q_ref)
  C <- unlist(qc_comp)
  D <- unlist(q_comp)
  x <- intersect(A, setdiff(D, C))
  y <- setdiff(C, intersect(C, setdiff(B, A)))
  z <- union(x, y)
  d <- intersect(A, z)
  r <- ifelse(length(z) == 0, NA, length(d))
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
JI_core_class <- JI_core_class %>% rowwise() %>% mutate(union = new_union(qc_ref, q_ref, qc_comp, q_comp))
JI_core_class <- JI_core_class %>% rowwise() %>% mutate(intersect = new_intersect2(qc_ref, q_ref, qc_comp, q_comp))
JI_core_class <- JI_core_class %>% 
  mutate(ref_n_class = length(qc_ref), comp_n_class = length(qc_comp), ref_n_all = length(q_ref), comp_n_all = length(q_comp))
JI_core_class <- JI_core_class %>% 
  mutate(jaccard = intersect / union)
JI_core_class <- JI_core_class %>% 
  ungroup() %>% filter(tool_ref != tool_comp) %>% 
  select(-c(qc_ref, qc_comp, q_ref, q_comp))


data.frame(JI_core_class %>% filter(new_level %in% top20) )


JI_core_class %>% ungroup() %>% group_by(new_level, tool_ref) %>%
  rowwise() %>%
  mutate(t1 = min(as.character(tool_ref), as.character(tool_comp)), 
         t2 = max(as.character(tool_ref), as.character(tool_comp))) %>%
  ungroup() %>% group_by(new_level, t1, t2) %>% slice_head(n = 1) %>% 
  filter(new_level %in% c("TET - RPG", "ERM", "Class A")) %>%
  ungroup() %>% group_by(new_level) %>% summarise(mn = mean(jaccard), md = median(jaccard))

data.frame(JI_core_class %>% filter(new_level %in% c("TET - RPG", "ERM", "Class A")))

JI_core_class %>% ungroup() %>% group_by(new_level, tool_ref) %>%
  rowwise() %>%
  mutate(t1 = min(as.character(tool_ref), as.character(tool_comp)), 
         t2 = max(as.character(tool_ref), as.character(tool_comp))) %>%
  ungroup() %>% group_by(new_level, t1, t2) %>% slice_head(n = 1) %>% 
  ungroup() %>% group_by(new_level) %>% summarise(mn = mean(jaccard), md = median(jaccard))


core_jaccard_human <- JI_core_class %>% ungroup() %>% #complete(tool, new_level, fill = list(unigenes = 0, total = 0)) %>%
  ungroup() %>% 
  filter(as.numeric(tool_ref) > as.numeric(tool_comp)) %>% 
  ggplot(aes(x = tool_ref, y = tool_comp, fill = jaccard)) +
  geom_tile() +
  scale_x_discrete(labels = tool_label) +
  scale_y_discrete(labels = tool_label) +
  facet_wrap(new_level ~ .) +
  scale_fill_viridis_c() + 
  theme_minimal() +
  labs(fill = "") +
  xlab("") +
  ylab("") +
  ggtitle("") + 
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = general_size ),
    panel.border = element_rect(fill = "transparent", color = "black", linewidth = 1),
    #panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    plot.margin = margin(0, 0, 0, 0, unit = "pt"),
    legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
    legend.margin = margin(0, 0, 0, 0, unit = "pt"),
    panel.spacing = unit(0, "pt"),
    strip.text = element_text(size = general_size + 1, face = "bold"),
    axis.title = element_text(size = general_size + 1, face = "bold"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = general_size),
    axis.text.y = element_text(size = general_size))

core_jaccard_human

ggsave("~/Documents/plots_project2/jaccard_class_core.svg", core_jaccard_human, width = 300, height = 300, unit = "mm") 


###

sets0 <- core %>% filter(cut %in% 0.5 & cnt > 900, habitat %in% c( "human gut" )) %>%
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

new_union2 <- function( q_ref, q_comp){
  A <- unlist(q_ref)
  B <- unlist(q_comp)
  d <- union(A, B)
  r <- ifelse(length(d) == 0, 0, length(d))
  return(r)
}

JI_core_overlap <- pairwise %>%
  left_join(sets0, by = c("tool_ref" = "tool")) %>%
  rename(q_ref = query) %>%
  left_join(sets0, by = c("tool_comp" = "tool")) %>%
  rename(q_comp = query)

JI_core_overlap <- JI_core_overlap %>% rowwise() %>% 
  mutate(inter = new_intersect2(q_ref, q_comp), union = new_union2(q_ref, q_comp)) %>% 
  mutate(ref_n = length(q_ref), comp_n = length(q_comp)) %>% select(-c(q_ref, q_comp)) 

JI_core_overlap <- JI_core_overlap %>% 
  mutate(recall = inter / comp_n, diff = (ref_n - inter)/ref_n, Jac = inter / union)

data.frame(JI_core_overlap)

JI_core_overlap %>% filter(tool_comp %in% "fARGene (nt)") %>% select(-c(diff))
JI_core_overlap %>% filter(tool_ref %in% "fARGene (nt)") 

JI_core_overlap %>% ungroup() %>% group_by(tool_ref) %>% summarise(md = median(Jac))



sumcore_human_class

sumpan_human_class <- sumpan %>% filter(habitat %in% "human gut") %>% 
  ungroup() %>% group_by(tool) %>% mutate(total = sum(mn)) %>% 
  mutate(p = mn / total )

complete(sample, gene, tool, 
         fill = list(normed10m = 0, scaled = 0, raw = 0, raw_unique = 0, unigenes = 0))



core_human_medians_heatmap <- sumcore_human_class %>% ungroup() %>% complete(tool, new_level, fill = list(unigenes = 0, total = 0)) %>%
  ungroup() %>% 
  group_by(new_level) %>% 
  mutate(mx =  max(unigenes, na.rm = T)) %>% 
  filter(mx > 0) %>% 
  mutate(norm  = log(unigenes + 1) / log(mx + 1)) %>%
  ggplot(aes(x = tool, y = new_level, fill = norm)) +
  geom_tile() +
  scale_x_discrete(labels = tool_label) +
  scale_fill_viridis_c() + 
  theme_minimal() +
  labs(fill = "") +
  xlab("") +
  ylab("Gene class") +
  ggtitle("B") + 
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
    axis.title = element_text(size = general_size + 1, face = "bold"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = general_size),
    axis.text.y = element_text(size = general_size))



pan_human_medians_heatmap <- sumpan_human_class %>% ungroup() %>% complete(tool, gene_class, fill = list(md = 0, total = 0)) %>%
  ungroup() %>% 
  group_by(gene_class) %>% 
  mutate(mx =  max(md, na.rm = T)) %>% 
  filter(mx > 0) %>% 
  mutate(norm  = log(md + 1) / log(mx + 1)) %>%
  ggplot(aes(x = tool, y = gene_class, fill = norm)) +
  geom_tile() +
  scale_x_discrete(labels = tool_label) +
  scale_fill_viridis_c() + 
  theme_minimal() +
  labs(fill = "") +
  xlab("") +
  ylab("Gene class") +
  ggtitle("C") + 
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
    axis.title = element_text(size = general_size + 1, face = "bold"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = general_size),
    axis.text.y = element_text(size = general_size))

ggsave("~/Documents/plots_project2/core_human_heatmap.svg", core_human_medians_heatmap, width = 60, height = 100, unit = "mm")



sumcore_pig_class <- sumcore %>% filter(habitat %in% "pig gut") %>% 
  ungroup() %>% group_by(tool) %>% mutate(total = sum(unigenes)) %>% 
  mutate(p = unigenes / total ) %>% arrange(desc(unigenes))

sumcore_pig_class %>% filter(new_level %in% "GPA")

sumcore_human_class %>% filter(new_level %in% "GPA")

sumcore_human_class %>% arrange(desc(unigenes))

sumcore %>% filter(habitat %in% c("human gut", "pig gut")) %>% ungroup() %>% 
  group_by(tool, habitat) %>% summarise(t = sum(unigenes)) %>% mutate(pig = t[2]) %>% 
  mutate(pig / t)

####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################


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





#################
#################
# Venn diagrams


card.ven <- ggVennDiagram(list(a1 = lst$rgi.diamond$query,
                             a2 = lst$abricate.card.norm$query),
                        color = 1, lwd = 0.7, label_size = 2, set_size = 2) + 
  scale_fill_gradientn(colors = c(pal_10_q[3], pal_10_q[4]), "white") +
  theme(legend.position = "none") +
  ggtitle("") +
  theme(
    legend.text = element_text(size = general_size),
    legend.title = element_text(size = general_size + 1),
    text = element_text(size = general_size),
    axis.text  = element_blank(),
    title = element_text(size = general_size + 2, , face = "bold"))

card.ven + coord_flip()

card.ven@region$fill_color <- c(
  pal_10_q[3],     # color for a1 only
  pal_10_q[4],     # color for a2 only
  "white"          # color for a1 ∩ a2
)

rf.ven = ggVennDiagram(list(a1 = lst$resfinder.norm$query,
                              a2 = lst$abricate.resfinder.norm$query),
                         color = 1, lwd = 0.7, label_size = 2, set_size = 2) + 
  scale_fill_gradientn(colors = c(pal_10_q[7], pal_10_q[8]), "white") +
  theme(legend.position = "none") +
  ggtitle("") +
  theme(
    legend.text = element_text(size = general_size),
    legend.title = element_text(size = general_size + 1),
    text = element_text(size = general_size),
    axis.text  = element_blank(),
    title = element_text(size = general_size + 2, , face = "bold"))

rf.ven + coord_flip()
ggsave("~/Documents/plots_project2/rf.ven.svg", rf.ven + coord_flip(), width = 70, height = 20, unit = "mm")

fit <- euler(c('A' = sum(!lst$abricate.resfinder.norm$query %in% lst$resfinder.norm$query), 
               'B' = sum(!lst$resfinder.norm$query %in% lst$abricate.resfinder.norm$query), 
               'A&B' = sum(lst$resfinder.norm$query %in% lst$abricate.resfinder.norm$query)))

ggsave(file="~/Documents/plots_project2/rf.ven.svg", 
       plot=plot(fit, fill=c(pal_10_q[8], pal_10_q[7], "white"), 
                 labels = list(cex = 0.5), quantities = list(cex = 0.5)), width = 150, height = 20, unit = "mm")

fit <- euler(c('A' = sum(!lst$abricate.card.norm$query %in% lst$rgi.diamond$query), 
               'B' = sum(!lst$rgi.diamond$query %in% lst$abricate.card.norm$query), 
               'A&B' = sum(lst$rgi.diamond$query %in% lst$abricate.card.norm$query)))

ggsave(file="~/Documents/plots_project2/card.ven.svg", 
       plot=plot(fit, fill=c(pal_10_q[4], pal_10_q[3], "white"), 
                 labels = list(cex = 0.5), quantities = list(cex = 0.5)), width = 150, height = 20, unit = "mm")



fit <- euler(c('A' = sum(!lst$abricate.ncbi.norm$query %in% lst$amrfinder.norm.prot$query), 
               'B' = sum(!lst$amrfinder.norm.prot$query %in% lst$abricate.ncbi.norm$query), 
               'A&B' = sum(lst$amrfinder.norm.prot$query %in% lst$abricate.ncbi.norm$query)))

ggsave(file="~/Documents/plots_project2/ncbi.ven.svg", 
       plot=plot(fit, fill=c(pal_10_q[6], pal_10_q[5], "white"), 
                 labels = list(cex = 0.5), quantities = list(cex = 0.5)), width = 150, height = 20, unit = "mm")



rgi.ven = ggVennDiagram(list(a1 = lst$rgi.diamond$query,
                             a2 = lst$rgi.blast$query,
                             b1 = lst$rgi.diamond.prot$query),
                        color = 1, lwd = 0.7, label_size = 2, set_size = 2) + 
  scale_fill_gradient(low = "#F4FAFE", high = pal_10_q[3]) +
  theme(legend.position = "none") +
  ggtitle("A") +
  theme(
    legend.text = element_text(size = general_size),
    legend.title = element_text(size = general_size + 1),
    text = element_text(size = general_size),
    axis.text  = element_blank(),
    title = element_text(size = general_size + 2, , face = "bold"))

rgi.ven

rgi.diamond.ven = ggVennDiagram(list(a = lst$rgi.diamond$query,
                                     b = lst$rgi.diamond.prot$query),
                                color = 1, lwd = 0.7, label_size = 2, set_size = 2) + 
  scale_fill_gradient(low = "#F4FAFE", high = pal_10_q[3]) +
  theme(legend.position = "none") +
  ggtitle("RGI") +
  theme(
    legend.text = element_text(size = general_size),
    legend.title = element_text(size = general_size + 1),
    text = element_text(size = general_size),
    axis.text  = element_blank(),
    title = element_text(size = general_size + 2, , face = "bold"))

rgi.diamond.ven

deeparg.ven = ggVennDiagram(list(a = lst$deeparg.norm$query,
                                 b = lst$deeparg.norm.prot$query),
                                 color = 1, lwd = 0.1, label_size = 2, set_size = 2) + 
  scale_fill_gradient(low = "#F4FAFE", high = pal_10_q[1]) +
  theme(legend.position = "none") +
  ggtitle("B") +
  theme(
    legend.text = element_text(size = general_size),
    legend.title = element_text(size = general_size + 1),
    text = element_text(size = general_size),
    axis.text  = element_blank(),
    title = element_text(size = general_size + 2, , face = "bold"))

deeparg.ven

fargene.ven = ggVennDiagram(list(a = lst$fargene$query,
                                 b = lst$fargene.prot$query),
                            color = 1, lwd = 0.1, label_size = 2, set_size = 2) + 
  scale_fill_gradient(low = "#F4FAFE", high = pal_10_q[2]) +  theme(legend.position = "none") +
  ggtitle("C") +
  theme(
    legend.text = element_text(size = general_size),
    legend.title = element_text(size = general_size + 1),
    text = element_text(size = general_size),
    axis.text  = element_blank(),
    title = element_text(size = general_size + 2, , face = "bold"))

fargene.ven

amrfinder.ven = ggVennDiagram(list(a = lst$amrfinder.norm$query,
                                   b = lst$amrfinder.norm.prot$query),
                              color = 1, lwd = 0.1, label_size = 2, set_size = 2) + 
  scale_fill_gradient(low = "#F4FAFE", high = pal_10_q[4]) +  theme(legend.position = "none") +
  ggtitle("D") +
  theme(
    legend.text = element_text(size = general_size),
    legend.title = element_text(size = general_size + 1),
    text = element_text(size = general_size),
    axis.text  = element_blank(),
    title = element_text(size = general_size + 2, , face = "bold"))

amrfinder.ven

abricate.ven = ggVennDiagram(list(argannot = lst$abricate.argannot.norm$query,
                                  card = lst$abricate.card.norm$query,
                                  megares = lst$abricate.megares.norm$query,
                                  ncbi = lst$abricate.ncbi.norm$query,
                                  resfinder = lst$abricate.resfinder.norm$query),
                                  color = 1, lwd = 0.1, label_size = 2, set_size = 2) + 
  scale_fill_gradient(low = "#F4FAFE", high = pal_10_q[2]) +  theme(legend.position = "none") +
  ggtitle("abricate") +
  theme(
    legend.text = element_text(size = general_size),
    legend.title = element_text(size = general_size + 1),
    text = element_text(size = general_size),
    axis.text  = element_blank(),
    title = element_text(size = general_size + 2, , face = "bold"))


ggsave("~/Documents/plots_project2/venn_diag.svg", grid.arrange(rgi.ven, deeparg.ven, fargene.ven, amrfinder.ven, nrow = 2), width = 70, height = 70, unit = "mm")




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







########


abundance  %>% 
  group_by(tool, sample, aggregation) %>% 
  summarise(total = sum(normed10m)) %>%
  ungroup() %>% group_by(tool, aggregation) %>%
  summarise(md = median(total)+1, q25 = quantile(total, 0.25) + 1, q75 = quantile(total, 0.75) + 1) 


abundance  %>% 
  group_by(tool, sample, aggregation) %>% 
  summarise(total = sum(unigenes)) %>%
  ungroup() %>% group_by(tool, aggregation) %>%
  summarise(md = median(total)+1, q25 = quantile(total, 0.25) + 1, q75 = quantile(total, 0.75) + 1) 


abundance_medians_h2 <- abundance  %>% 
  group_by(habitat2, location, tool, sample, aggregation) %>% 
  summarise(total = sum(normed10m)) %>%
  ungroup() %>% group_by(habitat2, location, tool, aggregation) %>%
  summarise(md = median(total)+1, q25 = quantile(total, 0.25) + 1, q75 = quantile(total, 0.75) + 1) 

diversity_medians_h2 <- abundance  %>% 
  group_by(habitat2, location, tool, sample, aggregation) %>% 
  summarise(total = sum(unigenes)) %>%
  ungroup() %>% group_by(habitat2, location, tool, aggregation) %>%
  summarise(md = median(total)+1, q25 = quantile(total, 0.25) + 1, q75 = quantile(total, 0.75) + 1)


abundance_medians_h2 <- abundance_medians_h2 %>% mutate(habitat2_factor = as.numeric(habitat2)) %>%
  mutate(habitat2_factor = ifelse(tool %in% c("RGI (DIAMOND - nt)"), habitat2_factor + 0.3, habitat2_factor)) %>%
  mutate(habitat2_factor = ifelse(tool %in% c("ABRicate (CARD - nt)"), habitat2_factor + 0.35, habitat2_factor)) %>%
  mutate(habitat2_factor = ifelse(tool %in% c("AMRFinderPlus (aa)"), habitat2_factor + 0.5, habitat2_factor)) %>%
  mutate(habitat2_factor = ifelse(tool %in% c("ABRicate (NCBI - nt)"), habitat2_factor + 0.55, habitat2_factor)) %>%
  mutate(habitat2_factor = ifelse(tool %in% c("ResFinder (nt)"), habitat2_factor + 0.65, habitat2_factor)) %>% 
  mutate(habitat2_factor = ifelse(tool %in% c("ABRicate (ResFinder - nt)"), habitat2_factor + 0.7, habitat2_factor)) %>% 
  mutate(habitat2_factor = ifelse(tool %in% c("fARGene (nt)"), habitat2_factor + 0.05, habitat2_factor)) %>%
  mutate(habitat2_factor = ifelse(tool %in% c("ABRicate (ARG-ANNOT - nt)"), habitat2_factor + 0.1, habitat2_factor))  %>%
  mutate(habitat2_factor = ifelse(tool %in% c("ABRicate (MEGARES - nt)"), habitat2_factor + 0.15, habitat2_factor)) 

diversity_medians_h2 <- diversity_medians_h2 %>% mutate(habitat2_factor = as.numeric(habitat2)) %>%
  mutate(habitat2_factor = ifelse(tool %in% c("RGI (DIAMOND - nt)"), habitat2_factor + 0.3, habitat2_factor)) %>%
  mutate(habitat2_factor = ifelse(tool %in% c("ABRicate (CARD - nt)"), habitat2_factor + 0.35, habitat2_factor)) %>%
  mutate(habitat2_factor = ifelse(tool %in% c("AMRFinderPlus (aa)"), habitat2_factor + 0.5, habitat2_factor)) %>%
  mutate(habitat2_factor = ifelse(tool %in% c("ABRicate (NCBI - nt)"), habitat2_factor + 0.55, habitat2_factor)) %>%
  mutate(habitat2_factor = ifelse(tool %in% c("ResFinder (nt)"), habitat2_factor + 0.65, habitat2_factor)) %>% 
  mutate(habitat2_factor = ifelse(tool %in% c("ABRicate (ResFinder - nt)"), habitat2_factor + 0.7, habitat2_factor)) %>% 
  mutate(habitat2_factor = ifelse(tool %in% c("fARGene (nt)"), habitat2_factor + 0.05, habitat2_factor)) %>%
  mutate(habitat2_factor = ifelse(tool %in% c("ABRicate (ARG-ANNOT - nt)"), habitat2_factor + 0.1, habitat2_factor))  %>%
  mutate(habitat2_factor = ifelse(tool %in% c("ABRicate (MEGARES - nt)"), habitat2_factor + 0.15, habitat2_factor)) 

abundance_plot_2 <- abundance_medians_h2  %>% 
  ggplot(aes( x = habitat2_factor, y = md, fill = tool, color = tool, group = interaction(tool, location))) +
  geom_rect(data = rects, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, alpha = alpha), fill = "grey", inherit.aes = FALSE) +
  geom_line(linetype = 3, linewidth = 0.7, show.legend = FALSE) +
  geom_rect(aes(
    xmin = as.numeric(habitat2_factor) - 0.03,
    xmax = as.numeric(habitat2_factor) + 0.03,
    ymin = q25,
    ymax = q75), color = "black", linewidth = 0.2) + 
  geom_point(shape = 15, color = "black") +
  scale_color_manual(values = pal_10_q, labels = tool_label) +
  scale_fill_manual(values = pal_10_q, labels = tool_label) +
  theme_minimal() +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)),
                limits = c(1, 27300), expand = c(0, 0)) + 
  scale_x_continuous(limits = c(0.9, 7), expand = c(0, 0), breaks = c(1.5, 2.5, 3.5, 4.5, 5.5, 6.5), labels = levels(abundance_medians_h2$habitat2)) + 
  scale_alpha(range = c(0, 0.3), guide = "none") +
  ylab("Abundance") +
  xlab("") +
  labs(fill = "") +
  ggtitle("B") +
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
    title = element_text(size = general_size + 2, face = "bold"),
    axis.title = element_text(size = general_size + 1, face = "bold"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = general_size),
    axis.text.y = element_text(size = general_size))

abundance_plot_2



diversity_plot_2 <- diversity_medians_h2  %>% 
  ggplot(aes( x = habitat2_factor, y = md, fill = tool, color = tool, group = interaction(tool, location))) +
  geom_rect(data = rects2, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, alpha = alpha), fill = "grey", inherit.aes = FALSE) +
  geom_line(linetype = 3, linewidth = 0.7, show.legend = FALSE) +
  geom_rect(aes(
    xmin = as.numeric(habitat2_factor) - 0.03,
    xmax = as.numeric(habitat2_factor) + 0.03,
    ymin = q25,
    ymax = q75), color = "black", linewidth = 0.2) + 
  geom_point(shape = 15, color = "black") +
  scale_color_manual(values = pal_10_q, labels = tool_label) +
  scale_fill_manual(values = pal_10_q, labels = tool_label) +
  theme_minimal() +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)),
                limits = c(1, 1400), expand = c(0, 0)) + 
  scale_x_continuous(limits = c(0.9, 7), expand = c(0, 0), breaks = c(1.5, 2.5, 3.5, 4.5, 5.5, 6.5), labels = levels(abundance_medians_h2$habitat2)) + 
  scale_alpha(range = c(0, 0.3), guide = "none") +
  ylab("Diversity") +
  xlab("") +
  labs(fill = "") +
  ggtitle("C") +
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
    title = element_text(size = general_size + 2, face = "bold"),
    axis.title = element_text(size = general_size + 1, face = "bold"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = general_size),
    axis.text.y = element_text(size = general_size))

diversity_plot_2


df <- data.frame(
  x = 1:7,
  y = c(1:7)
)

max_y <- max(27300)

# Create a data frame for the rectangles
rects <- data.frame(
  xmin = df$x - .1,
  xmax = df$x + .9,
  ymin = 1,
  ymax = max_y,
  alpha = rep( c(0, .1), length.out = nrow(df))
)

rects2 <- data.frame(
  xmin = df$x - .1,
  xmax = df$x + .9,
  ymin = 1,
  ymax = 1400,
  alpha = rep( c(0, .1), length.out = nrow(df))
)

unigenes



abundance_3 <- abundance  %>% 
  group_by(tool, aggregation) %>% 
  summarise(total = sum(normed10m), 
            md = median(normed10m), 
            mn = mean(normed10m), 
            q25 = quantile(normed10m, 0.25) , 
            q75 = quantile(total, 0.75))

abundance_3 <- abundance_3 %>% left_join(unigenes  %>% group_by(tool) %>% summarise(n = n_distinct(query)) %>% ungroup() %>% arrange(n))


ggplot(abundance_3, aes(x = log(n), y = log(total), color = tool, fill = tool)) +
  geom_point() + 
  geom_smooth(aes(x = log(n), y = log(total)), method = "lm", se = FALSE, color = "black", linetype = "dashed", inherit.aes = FALSE) 

plot(log(abundance_3$n), log(abundance_3$total))
abline(a = 0, b = 1, col = "red", lty = 2)




ggsave("~/Documents/plots_project2/abundance_2.svg", abundance_plot_2 + theme( legend.position = "none"), width = 180, height = 70, unit = "mm")
ggsave("~/Documents/plots_project2/diversity_2.svg", diversity_plot_2 + theme( legend.position = "none"), width = 180, height = 70, unit = "mm")
ggsave("~/Documents/plots_project2/abundance_legend_2.svg", g_legend(abundance_plot_2), width = 180, height = 15, unit = "mm")




hab_plot <- unique(abundance_class$habitat)
hab_plot <- hab_plot[hab_plot != "human gut"] 
names(hab_plot) <- stringr::str_replace_all(hab_plot, " ", "")


for(j in 1:length(hab_plot)){
  abundance_medians_env <- abundance_class %>% ungroup() %>%
    group_by(tool, gene) %>% filter(habitat %in% hab_plot[j]) %>%
    summarise(md = median(normed10m), q75 = quantile(normed10m, 0.75), q25 = quantile(normed10m, 0.25))
  
  top_abundance_env <- abundance_medians_env %>% 
    ungroup() %>% 
    group_by(gene) %>% 
    mutate(mx =  max(md)) %>% 
    filter(md > 0) %>% select(gene) %>% distinct() %>% pull()
  
  top20_env <- abundance_class %>% 
    filter(habitat %in% hab_plot[j], tool %in% tool_2) %>% 
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
  
  top20_env <- unique(c(top20_env, "Cell wall charge", "ERM", "Class C", "Target-modifying enzyme", "ABC-F"))
  top20_env <- factor(c(top20_env, "Other"), levels = c(top20_env, "Other"))
  
  top_abundance_env <- unique(c(top_abundance_env, as.character(top20_env)))
  
  abundance_medians_env_non_zero <- abundance_medians_env %>% filter(q75 > 0) 
  
  abundance_class_plot_env <- abundance_class  %>% 
    filter(habitat %in% hab_plot[j]) %>% 
    ungroup() %>%
    mutate(g2 = ifelse(gene %in% top_abundance_env, gene, "Other")) %>% 
    mutate(g2 = factor(g2, levels = top_abundance_env)) %>% 
    group_by(habitat, tool, g2, sample) %>% summarise(total = sum(normed10m)) %>%
    filter(paste(g2, tool) %in% paste(abundance_medians_env_non_zero$gene, abundance_medians_env_non_zero$tool)) %>%
    ggplot(aes( x = g2, y = total +1, fill = tool)) +
    geom_boxplot(position = position_dodge2(preserve = "single"),  outlier.shape = NA, coef = 0, width = 1, linewidth = 0.2) +
    scale_fill_manual(values = pal_10_q, labels = tool_label) +
    facet_grid(. ~ g2, scales = "free_x") +
    theme_minimal() +
    scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x))) + #,
    #limits = c(1, 10000)) + 
    coord_cartesian(ylim = c(1, 5000)) +
    ylab("Abundance") +
    xlab("") +
    labs(fill = "") +
    ggtitle("") +
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
      strip.text.x = element_blank(),
      title = element_text(size = general_size + 2, face = "bold"),
      strip.text.y = element_blank(),
      axis.title = element_text(size = general_size + 1, face = "bold"),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = general_size),
      axis.text.y = element_text(size = general_size)) 
  
  abundance_class_plot_env
  
  ggsave(paste0("~/Documents/plots_project2/fig_s_abundance_",names(hab_plot)[j],".svg"), abundance_class_plot_env, width = 180, height = 70, unit = "mm")
}
 


for(j in 1:length(hab_plot)){
  env_pan <- sumpan %>% filter(habitat %in% hab_plot[j]) %>% ungroup() %>% 
    mutate(gene_class = ifelse(as.character(gene_class) %in% levels(top20), gene_class, "Other")) %>% 
    mutate(pattern = ifelse(as.character(gene_class) %in% levels(top20)[1:round(length(top20)/2)], "yes", "no")) %>% 
    mutate(gene_class = factor(gene_class, levels = top20)) %>% ungroup() %>% 
    mutate(d = "Pan-resistome") %>% 
    rename(new_level = gene_class, unigenes = md) %>% select(tool, habitat, new_level, unigenes, pattern, d)
  
  env_core <- sumcore %>% filter(habitat %in% hab_plot[j]) %>% ungroup() %>% 
    mutate(unigenes = as.numeric(unigenes)) %>% 
    mutate(new_level = ifelse(as.character(new_level) %in% levels(top20), new_level, "Other")) %>% 
    mutate(pattern = ifelse(as.character(new_level) %in% levels(top20)[1:round(length(top20)/2)], "yes", "no")) %>% 
    mutate(d = "Core-resistome") %>%
    mutate(new_level = factor(new_level, levels = top20))
  
  top20_env_2 <- sumcore %>% filter(habitat %in% hab_plot[j]) %>% ungroup() %>% 
    arrange(desc(unigenes)) %>% group_by(tool) %>% slice_head(n = 6) %>% ungroup() %>% select(new_level) %>% distinct() %>% pull()
  
  env_resistome <- bind_rows(env_pan, env_core) %>% 
    mutate(tool = as.character(tool)) %>% 
    mutate(tool = factor(tool, levels = tools_levels)) 
  
  lim_min_plot <- ifelse(max(env_resistome$unigenes[env_resistome$d %in% "Core-resistome"]) < 50, 0, 
                         ifelse(max(env_resistome$unigenes[env_resistome$d %in% "Core-resistome"]) > 200, -4, 0))
  
  lim_max_plot <- ifelse(max(env_resistome$unigenes[env_resistome$d %in% "Core-resistome"]) < 20, max(env_resistome$unigenes[env_resistome$d %in% "Core-resistome"]), 
                         max(env_resistome$unigenes[env_resistome$d %in% "Core-resistome"]) + 2)
  
  p_core_env <- env_resistome %>% filter(d %in% "Core-resistome") %>% 
    group_by(tool, habitat, new_level) %>% 
    summarise(unigenes = sum(unigenes), d = d[1], pattern = pattern[1]) %>%
    ggplot(aes(y = unigenes, x = new_level, fill = tool)) +
    #geom_jitter(height = 0.6, size = 2) +
    geom_col(position = position_dodge2(preserve = "single"), color = "black", linewidth = 0.1) +
    theme_minimal() +
    labs(color = "") +
    facet_grid(. ~ new_level, scales = "free_x") + 
    xlab("") +
    ylab("") +
    ggtitle("") +
    scale_fill_manual(values = c(pal_10_q), labels = tool_label) +
    scale_y_continuous(limits = c(lim_min_plot, lim_max_plot), expand = c(0, 0)) +
    theme(legend.position = "bottom",
          panel.border = element_rect(fill = "transparent", color = "black", linewidth = 1),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = general_size),
          axis.text.y = element_text(size = general_size),
          #strip.text = element_text(size = general_size, face = "bold"),
          strip.text = element_blank(),
          plot.margin = margin(0, 0, 0, 0, unit = "pt"),
          legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
          legend.margin = margin(0, 0, 0, 0, unit = "pt"),
          panel.spacing = unit(0, "pt"),
          legend.text = element_text(size = general_size),
          title = element_text(size = general_size + 2, face = "bold"),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          axis.title = element_text(size = general_size + 1, face = "bold"))
  
  p_core_env
  
  ggsave(paste0("~/Documents/plots_project2/fig_s_core_",names(hab_plot)[j],".svg"), p_core_env, width = 180, height = 70, unit = "mm")

}
