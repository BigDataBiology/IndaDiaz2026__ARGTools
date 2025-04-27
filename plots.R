#rgi.blast
#rgi.diamond
#deeparg.norm
#fargene
#amrfinder.norm.prot
#abricate.argannot.norm
#abricate.card.norm
#abricate.megares.norm
#abricate.ncbi.norm
#abricate.resfinder.norm
#resfinder.norm

library(dplyr)
library(ggplot2)
library(ggVennDiagram)
library(gridExtra)
library(tidyverse)

new_level <- c("chloramphenicol phosphotransferase", "CHL ph.",
               "ciprofloxacin phosphotransferase", "CIP ph.",
               "viomycin phosphotransferase", "VIO ph. ",
               "major facilitator superfamily (MFS) antibiotic efflux pump", "MFS - efflux p.",
               "class C beta-lactamase", "Class C",
               "aminoglycoside bifunctional resistance protein", "Aminoglycoside b.",
               "streptothricin acetyltransferase (SAT)", "SAT",
               "AAC(2')", "AAC",
               "AAC(6')", "AAC",
               "AAC(3)", "AAC",
               "APH(3')", "APH",
               "APH(6)" , "APH",
               "class A beta-lactamase", "Class A",
               "macrolide phosphotransferase (MPH)", "MPH",
               "class B (metallo-) beta-lactamase", "Class B",
               "class D beta-lactamase", "Class D",
               "quinolone resistance protein (qnr)", "QNR",
               "Erm 23S ribosomal RNA methyltransferase", "ERM",
               "APH(2'')", "APH",
               "tetracycline inactivation enzyme", "TET - enzyme",
               "tetracycline-resistant ribosomal protection protein", "TET - RPG",
               "gene(s) or protein(s) associated with a glycopeptide resistance cluster", "GPA",
               "protein(s) and two-component regulatory system modulating antibiotic efflux", "Efflux p.",
               "fosfomycin inactivation enzyme", "FO",
               "antibiotic resistant dihydrofolate reductase", "DHFR",
               "rifampin-resistant RNA polymerase-binding protein", "RIF",
               "antibiotic resistant gene variant or mutant", "V/M",            
               "gene involved in antibiotic sequestration", "Sequestration",
               "gene modulating beta-lactam resistance", "Beta-lactam modulation",
               "rifampin inactivation enzyme", "RIF",
               "nitroimidazole reductase", "NTR",
               "cpa acetyltransferase", "CPA ac.",
               "protein(s) conferring resistance via host-dependent nutrient acquisition", "Host-dependent",
               "protein modulating permeability to antibiotic", "Permeability modulation",
               "chloramphenicol acetyltransferase (CAT)", "CAT",
               "streptogramin inactivation enzyme", "Streptogramin enzyme",
               "gene altering cell wall charge", "Cell wall charge",
               "ANT(2'')", "ANT",
               "ANT(9)", "ANT",
               "ANT(4')", "ANT",
               "APH(3'')", "APH",
               "lincosamide nucleotidyltransferase (LNU)", "LNU",
               "ANT(3'')", "ANT",
               "gene involved in self-resistance to antibiotic", "Self-resistance",
               "ANT(6)", "ANT",
               "APH(9)", "APH",
               "APH(7'')", "APH",
               "sulfonamide resistant sul", "SUL",
               "ABC-F ATP-binding cassette ribosomal protection protein", "ABC-F",
               "macrolide glycosyltransferase", "MGT",
               "macrolide esterase", "MEs",
               "fusidic acid inactivation enzyme", "FUS enzyme",
               "Bah amidohydrolase", "Bah amidohydrolase",
               "Target protecting FusB-type protein conferring resistance to Fusidic acid", "FUS protection",
               "APH(4)", "APH",
               "gene conferring resistance via absence", "Absence",
               "glycopeptide resistance gene cluster", "GPA",
               "beta-lactam resistant penicillin-binding proteins", "PBP",
               "Edeine acetyltransferase", "EdeQ",
               "rifamycin-resistant beta-subunit of RNA polymerase (rpoB)", "rpoB",
               "efflux pump complex or subunit conferring antibiotic resistance", "Efflux p.",
               "protein(s) conferring antibiotic resistance via molecular bypass", "Molecular bypass",
               "antibiotic target modifying enzyme", "Target-modifying enzyme",
               "antibiotic inactivation enzyme", "Inactivation enzyme")

odd_vals  <- new_level[seq(1, length(new_level), by = 2)]
even_vals <- new_level[seq(2, length(new_level), by = 2)]
new_level_df <- df <- data.frame(old = odd_vals, new = even_vals)
rm(new_level, even_vals, odd_vals)


df2 <- readRDS(file = "~/df2.rds")
abundance_parent <- readRDS("GitHub/arg_compare/abundance_per_tool.rds")
diversity_parent <- readRDS("GitHub/arg_compare/diversity_per_tool.rds")

#pal0 <- c("#8c510a","#bf812d","#dfc27d","#f6e8c3","#c7eae5","#80cdc1","#35978f","#01665e")
pal <- c("#543005", "#8c510a", "#bf812d", "#dfc27d", "#f6e8c3", "#f5f5f5", "#c7eae5", "#80cdc1", "#35978f", "#01665e", "#003c30")
pal2a <- c("#00B5E2", "#013B67", "#CC0C01", "#985428", "#73006D", "#E08728", "#377E60", "#FD8CC0")
pal3a <- c("#00B5E2", "#377E60", "#CC0C01", "#985428", "#73006D", "#E08728",  "#FD8CC0")
pal4 <- c("#a6611a", "#dfc27d", "#80cdc1", "#018571")
pal5 <- c("#a6611a", "#dfc27d", "grey30", "#80cdc1", "#018571")

# HABITATS
EN <- c("human gut", "human oral",  "human skin", "human nose", "human vagina", 
        "dog gut", "cat gut", "mouse gut", "pig gut", "wastewater", "marine", "freshwater",  
        "soil" , "amplicon", "isolate",  "built-environment" )

# SOURCE FOR EACH HABITAT
SO <- c(rep("Humans", 5), rep("Mammals", 4),  "Wastewater", "Marine", "Freshwater", "Soil", rep("Other", 3))
names(SO) <- EN

# changing habitats and tools to factor
abundance_parent$habitat <- factor(abundance_parent$habitat, levels = EN)
abundance_parent$habitat2 <- factor(SO[abundance_parent$habitat], levels = unique(SO))
abundance_parent$tool <- factor(abundance_parent$tool, 
                                levels = c("rgi.diamond", "rgi.diamond.prot", "rgi.blast", "deeparg", "deeparg.prot", 
                                "fargene", "fargene.prot", "resfinder", "amrfinder", "amrfinder.prot", "abricate.argannot", 
                                "abricate.card", "abricate.megares", "abricate.ncbi", "abricate.resfinder"))

diversity_parent$habitat <- factor(diversity_parent$habitat, levels = EN)
diversity_parent$habitat2 <- factor(SO[diversity_parent$habitat], levels = unique(SO))
diversity_parent$tool <- factor(diversity_parent$tool, 
                                levels = c("rgi.diamond", "rgi.diamond.prot", "rgi.blast", "deeparg", "deeparg.prot", 
                                           "fargene", "fargene.prot", "resfinder", "amrfinder", "amrfinder.prot", "abricate.argannot", 
                                           "abricate.card", "abricate.megares", "abricate.ncbi", "abricate.resfinder"))

# add ontology description
abundance_parent <- abundance_parent %>% mutate(parent_label = df2$Parent_Label[match(parent,df2$Parent_ID)])
diversity_parent <- diversity_parent %>% mutate(parent_label = df2$Parent_Label[match(parent,df2$Parent_ID)])

# add the condense ontology abbreviations
abundance_parent <- abundance_parent %>% mutate(new_level = new_level_df$new[match(parent_label, new_level_df$old)])
diversity_parent <- diversity_parent %>% mutate(new_level = new_level_df$new[match(parent_label, new_level_df$old)])

## factors for new_level, we take the highest abundance per ontology by tool
factor_new_level <- abundance_parent %>% ungroup() %>% 
  group_by(tool, new_level) %>% summarise(total = sum(normed10m)) %>%
  ungroup() %>% arrange(tool, desc(total)) %>% 
  group_by(tool, new_level) %>%  ungroup() %>% select(new_level) %>% distinct() %>% pull()

factor_new_level2 <- c(factor_new_level[seq(1, length(factor_new_level), by = 4)],
                       factor_new_level[seq(1, length(factor_new_level), by = 4) + 1],
                       factor_new_level[seq(1, length(factor_new_level), by = 4) + 2],
                       factor_new_level[seq(1, length(factor_new_level), by = 4) + 3])



abundance_parent$new_level <- factor(abundance_parent$new_level, levels = factor_new_level2)
diversity_parent$new_level <- factor(diversity_parent$new_level, levels = factor_new_level2)

# factor for habitat
abundance_parent$habitat <- factor(abundance_parent$habitat, levels = EN)
abundance_parent$habitat <- as.character(abundance_parent$habitat)

diversity_parent$habitat <- factor(diversity_parent$habitat, levels = EN)
diversity_parent$habitat <- as.character(diversity_parent$habitat)

# environments that we are not interested in
not_env <- c("built-environment", "amplicon", "isolate")
total_abundance_sample <- abundance_parent %>% group_by(sample) %>% summarise(total = sum(normed10m)) %>%
  arrange(desc(total))


# sample outliers
extreme_samples <- total_abundance_sample %>% filter(total > quantile(total, probs = .9977)) %>%
  ungroup() %>% select(sample) %>% pull()


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


lst <- readRDS("GitHub/arg_compare/results_tools.rds")

fargene_with_rgi <- read.delim("GitHub/arg_compare/check_missing_annot/fargene_predicted_fna.txt")
fargene_with_rgi <- fargene_with_rgi[,-c(1,18,19,20)]
fargene_with_rgi$query <- gsub(' ', '', fargene_with_rgi$Contig)
fargene_with_rgi$query <- gsub('.{2}$', '', fargene_with_rgi$query)
fargene_with_rgi$query <- gsub('tmp_', '', fargene_with_rgi$query)
fargene_with_rgi <- fargene_with_rgi %>%
  mutate(ARO = paste0("ARO:", ARO))
sum(fargene_with_rgi$query %in% lst$fargene$query)
sum(!fargene_with_rgi$query %in% lst$fargene$query)
sum(lst$fargene$query %in% fargene_with_rgi$query)
sum(!lst$fargene$query %in% fargene_with_rgi$query)

number_unigenes_per_tool <- sapply(lst, function(x) length(unique(x$query)))

names(lst)

lst$deeparg.norm.prot$id <- lst$deeparg.norm.prot$identity
lst$deeparg.norm$id <- lst$deeparg.norm$identity
#lst$fargene$id <- lst$fargene$hmm
lst$fargene$id <- fargene_with_rgi$Best_Identities[match(lst$fargene$query, fargene_with_rgi$query)]
lst$fargene.prot$id <- lst$fargene.prot$hmm
lst$rgi.blast$id <- lst$rgi.blast$Best_Identities
lst$rgi.diamond$id <- lst$rgi.diamond$Best_Identities
lst$rgi.diamond.prot$id <- lst$rgi.diamond.prot$Best_Identities
lst$amrfinder.norm$id <- lst$amrfinder.norm$X..Identity.to.reference
lst$amrfinder.norm.prot$id <- lst$amrfinder.norm.prot$X..Identity.to.reference
lst$abricate.argannot.norm$id <- lst$abricate.argannot.norm$V11
lst$abricate.card.norm$id <- lst$abricate.card.norm$V11
lst$abricate.megares.norm$id <- lst$abricate.megares.norm$V11
lst$abricate.ncbi.norm$id <- lst$abricate.ncbi.norm$V11
lst$abricate.resfinder.norm$id <- lst$abricate.resfinder.norm$V11
lst$resfinder.norm$id <- lst$resfinder.norm$Identity

unigenes_per_tool <- do.call(rbind, lapply(seq_along(1:length(lst)), function(j) data.frame(tool = rep(names(lst)[j], dim(lst[[j]])[1]), query = lst[[j]]$query, 
                                                                                            ARO = lst[[j]]$ARO, parent = lst[[j]]$parent, 
                                                                                            parent_description = lst[[j]]$parent_description,
                                                                                            id = lst[[j]]$id)))

unigenes_per_tool <- unigenes_per_tool %>% filter(!is.na(parent)) %>% distinct()


tool_2 <- c("deeparg.norm", "rgi.diamond", "fargene", "amrfinder.norm.prot", "abricate.argannot.norm",
            "abricate.card.norm", "abricate.megares.norm", "abricate.ncbi.norm", "abricate.resfinder.norm",
            "resfinder.norm")

tools_per_unigene <- unigenes_per_tool %>% ungroup() %>% filter(tool %in% tool_2) %>% 
  arrange(query) %>% 
  group_by(query) %>% 
  mutate(n_tools = n()) %>% 
  mutate(single = (n_tools ==1)) 


tools_per_unigene$new_level <- new_level_df$new[match(tools_per_unigene$parent_description, new_level_df$old)]
tools_per_unigene[duplicated(paste(tools_per_unigene$tool, tools_per_unigene$query)),]

########################################################################
#########################################################################
#########################################################################
# heatmaps

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


PM %>% filter(tool %in% "fargene")

PM_long <- PM %>% 
  pivot_longer(
    cols = -c(tool, new_level),
    names_to = "variable",
    values_to = "value"
  ) 

data.frame(PM_long %>% filter(tool %in% "fargene"))
PM_long %>% filter(tool %in% "fargene") %>% group_by(new_level) %>% summarise(n = sum(value))

PM_long %>% filter(variable %in% "fargene")
PM_long %>% filter(variable %in% "fargene") %>% group_by(new_level) %>% summarise(n = sum(value))

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


PM_new %>% filter(new_level %in% "GPA", tool %in% c("rgi.diamond", "deeparg.norm"), variable %in% c("rgi.diamond", "deeparg.norm"))


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
