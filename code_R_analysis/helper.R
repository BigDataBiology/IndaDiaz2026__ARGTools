# gene classes 
# gene_classes_list is a vector for labeling the plots 
# gene_classes is the list with the names from the prepared datasets and used for shiny

gene_classes = data.frame(rbind(
  c("glycopeptide resistance (van)" , "van"), 
  c("protein(s) and two-component regulatory system modulating antibiotic efflux" , "efflux pump"),
  c("gene altering cell wall charge" , "cell wall charge"), 
  c("rifamycin-resistant beta-subunit of RNA polymerase (rpoB)" , "rpoB"),
  c("class A beta-lactamase" , "class A beta-lactamase"), 
  c("class B beta-lactamase" , "class B beta-lactamase"), 
  c("class C beta-lactamase" , "class C beta-lactamase"), 
  c("class D beta-lactamase" , "class D beta-lactamase"),
  c("gene modulating beta-lactam resistance" , "beta-lactam modulation resistance"),
  c("aminoglycoside acetyltransferase (aac)" , "aac"), 
  c("aminoglycoside phosphotransferase (aph)" , "aph"), 
  c("aminoglycoside nucleotidyltransferase (ant)" , "ant"), 
  c("aminoglycoside bifunctional resistance protein" , "bifunctional aminoglycoside"),
  c("tetracycline-resistant ribosomal protection protein (tet RPG)" , "tet RPG"), 
  c("tetracycline inactivation enzyme (tet enzyme)" , "tet enzyme"),
  c("major facilitator superfamily antibiotic efflux pump (MFS efflux pump)" , "MFS efflux pump"), 
  c("erm 23S ribosomal RNA methyltransferase (erm)" , "erm"), 
  c("macrolide phosphotransferase (mph)"  , "mph"), 
  c("quinolone resistance protein (qnr)" , "qnr"),
  c("beta-lactam resistant penicillin-binding proteins (PBP)" , "PBP"), 
  c("antibiotic target modifying enzyme" , "target-modifying enzyme"),
  c("ciprofloxacin phosphotransferase (crpP)" , "crpP"), 
  c("rifampin-resistant RNA polymerase-binding protein (rifampin Rbp)" , "rifampin Rbp"),
  c("chloramphenicol phosphotransferase (cpt)" , "cpt"), 
  c("chloramphenicol acetyltransferase (cat)" , "cat"),
  c("sulfonamide resistant (sul)" , "sul"), 
  c("viomycin phosphotransferase (vph)" , "vph"),
  c("antibiotic resistant gene variant or mutant (variant or mutant)" , "variant or mutant"), 
  c("fosfomycin inactivation enzyme (fos)" , "fos"),
  c("antibiotic inactivation enzyme" , "antibiotic inactivation enzyme"), 
  c("antibiotic resistant dihydrofolate reductase (dfr)" , "dfr"),
  c("streptothricin acetyltransferase (sat)" , "sat"), 
  c("nitroimidazole reductase (nim)" , "nim"), 
  c("ABC-F ATP-binding cassette ribosomal protection protein (abcF)" , "abcF"), 
  c("protein modulating permeability to antibiotic" , "permeability modulation"),
  c("protein(s) conferring resistance via host-dependent nutrient acquisition" , "host-dependent nutrient acquisition"), 
  c("gene involved in antibiotic sequestration" , "antibiotic sequestration"),
  c("lincosamide nucleotidyltransferase (lnu)" , "lnu"), 
  c("macrolide glycosyltransferase (mgt)" , "mgt"),
  c("fusidic acid inactivation enzyme (fai)" , "fai"), 
  c("target protecting FusB-type protein conferring resistance to Fusidic acid (fusB-type)" , "fusB-type"), 
  c("macrolide esterase (mel)" , "mel"), 
  c("bah amidohydrolase (bah)" , "bah"),
  c("cpa acetyltransferase (cpa)" , "cpa"), 
  c("gene involved in self-resistance to antibiotic" , "self-resistance"),
  c("streptogramin inactivation enzym (vat)" , "vat"), 
  c("gene conferring resistance via absence" , "resistance by absence"),
  c("capreomycin phosphotransferase (cph)" , "cph"), 
  c("edeine acetyltransferase (edeQ)" , "edeQ"), 
  c("protein(s) conferring antibiotic resistance via molecular bypass" , "molecular bypass"),
  c("rifampin inactivation enzyme", "rifampin inactivation enzyme")))
colnames(gene_classes) <- c("old", "new")

#  gene_classes0 and gene_classes_list are the labels for the genes in the plots 
gene_classes0 <- gene_classes
gene_classes0 <- gene_classes0 %>% mutate(new = ifelse(new %in% "class A beta-lactamase", "A beta-lactamase", new)) %>% 
  mutate(new = ifelse(new %in% "class B beta-lactamase", "B beta-lactamase", new)) %>% 
  mutate(new = ifelse(new %in% "class C beta-lactamase", "C beta-lactamase", new)) %>%
  mutate(new = ifelse(new %in% "class D beta-lactamase", "D beta-lactamase", new)) %>%
  mutate(new = ifelse(new %in% "beta-lactam modulation resistance", "beta-lactam modulation", new))

gene_classes_list <- gene_classes0$new
rm(gene_classes0)
names(gene_classes_list) <- gene_classes$new

gene_classes <- setNames(as.list(gene_classes$new), gene_classes$old)

# class overlap and general overlap for all genes and for core resistome

# new_intersect <- function(qc_ref, q_ref, qc_comp, q_comp){
#   A <- unlist(qc_ref)
#   B <- unlist(q_ref)
#   C <- unlist(qc_comp)
#   D <- unlist(q_comp)
#   x <- intersect(A, setdiff(D, C))
#   y <- setdiff(C, intersect(C, setdiff(B, A)))
#   z <- union(x, y)
#   d <- intersect(A, z)
#   r <- ifelse(length(z) == 0, NA, length(d) / length(z))
#   return(r)
# }


new_intersect_lists <- function(qc_ref, q_ref, qc_comp, q_comp){
  A <- unlist(qc_ref)
  B <- unlist(q_ref)
  C <- unlist(qc_comp)
  D <- unlist(q_comp)
  comp_class <- C[!C %in% setdiff(B, A)] # remove from the class in the comparison tool those found in a different class in reference tool
  comp_class <- c(comp_class, D[D %in% intersect(setdiff(D, B), A)]) # complement class in comparison tool those genes found in a different class but they are in the right class in reference tool
  r <- ifelse(length(A) == 0 & length(comp_class) == 0, NA, 
              ifelse(length(A) == 0 & length(comp_class) != 0, NA, 
              ifelse(length(A) != 0 & length(comp_class) == 0, NA, 
                     length(intersect(A, comp_class)) / length(comp_class))))
  return(r)
}

new_difference_list <- function(qc_ref, q_ref, qc_comp, q_comp){
  A <-  unlist(qc_ref)
  B <- unlist(q_ref)
  C <- unlist(qc_comp)
  D <- unlist(q_comp)
  comp_class <- C[!C %in% setdiff(B, A)] # remove from the class in the comparison tool those found in a different class in reference tool
  comp_class <- c(comp_class, D[D %in% intersect(setdiff(D, B), A)]) # complement class in comparison tool those genes found in a different class but they are in the right class in reference tool
  r <- ifelse(length(A) == 0 & length(comp_class) == 0, NA, 
       ifelse(length(A) == 0 & length(comp_class) != 0, NA, 
       ifelse(length(A) != 0 & length(comp_class) == 0, NA, 
       length(setdiff(A, comp_class)) / length(A))))
  return(r)
}

create_class_overlaps <- function(unigenes){
  tools_per_unigene <- unigenes %>% ungroup()  %>% 
    arrange(query) %>% 
    group_by(query) %>% 
    mutate(n_tools = n_distinct(tool)) %>% 
    mutate(single = (n_tools ==1)) 
  
  sets0 <- tools_per_unigene %>%
    group_by(tool) %>%
    summarise(query = list(query), .groups = "drop")
  
  sets1 <- tools_per_unigene %>%
    group_by(new_level, tool) %>%
    summarise(query = list(query), .groups = "drop")  
  
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
  
  JI_class_other <- JI_class_other %>% 
    rowwise() %>%
    mutate(recall = new_intersect_lists(qc_ref, q_ref, qc_comp, q_comp)) %>% 
    rowwise() %>% 
    mutate(fnr = new_difference_list(qc_ref, q_ref, qc_comp, q_comp)) %>% 
    mutate(ref_n_class = length(qc_ref), comp_n_class = length(qc_comp), ref_n_all = length(q_ref), comp_n_all = length(q_comp)) %>% 
    ungroup() %>% 
    filter(tool_ref != tool_comp) %>% 
    select(-c(qc_ref, qc_comp, q_ref, q_comp))
  return(JI_class_other)
}

return_overlap_tools <- function(unigenes) { 
  sets <- unigenes %>% ungroup()  %>% 
    arrange(query) %>% 
    group_by(query) %>% 
    mutate(n_tools = n_distinct(tool)) %>% 
    mutate(single = (n_tools ==1))  %>%  
    group_by(tool) %>%
    summarise(query = list(query), .groups = "drop") # put every query in a list
  
  JI_all <- expand_grid(tool_ref = sets$tool, tool_comp = sets$tool)  %>%
    left_join(sets, by = c("tool_ref" = "tool")) %>%
    rename(values1 = query) %>%
    left_join(sets, by = c("tool_comp" = "tool")) %>%
    rename(values2 = query) %>%
    mutate(jaccard = map2_dbl(values1, values2, ~ length(intersect(.x, .y)) / length(union(.x, .y)))) %>%
    mutate(recall = map2_dbl(values1, values2, ~ length(intersect(.x, .y)) / length( .y))) %>%
    mutate(fnr = map2_dbl(values1, values2, ~ length(setdiff(.x, .y)) / length( .x))) %>%
    select(tool_ref, tool_comp, jaccard, recall, fnr) 
  
  return(JI_all)
  
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


create_class_overlaps_core <- function(core, cut_threshold, count_threshold, habitat){
  sets0 <- core %>% 
    filter(cut %in% cut_threshold & cnt > count_threshold, habitat %in% habitat) %>%
    group_by(tool) %>%
    summarise(query = list(X), .groups = "drop")
  
  sets1 <- core %>% 
    filter(cut %in% cut_threshold & cnt > count_threshold, habitat %in% habitat) %>%
    group_by(new_level, tool) %>%
    summarise(query = list(X), .groups = "drop")
  
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
  
  JI_core_class <- JI_core_class %>% 
    rowwise() %>% 
    mutate(recall = new_intersect_lists(qc_ref, q_ref, qc_comp, q_comp)) %>% 
    mutate(fnr = new_difference_list(qc_ref, q_ref, qc_comp, q_comp)) %>%
    mutate(union = new_union(qc_ref, q_ref, qc_comp, q_comp)) %>%
    mutate(intersect = new_intersect2(qc_ref, q_ref, qc_comp, q_comp)) %>%
    mutate(ref_n_class = length(qc_ref), comp_n_class = length(qc_comp), ref_n_all = length(q_ref), comp_n_all = length(q_comp)) %>%
    mutate(jaccard = intersect / union) %>% 
    ungroup() %>% 
    filter(tool_ref != tool_comp) %>% 
    select(-c(qc_ref, qc_comp, q_ref, q_comp))
  
  return(JI_core_class)
}


########################################################################
########################################################################
# Figure 1

plot_count_genes_tool <- function(unigenes, tools_for_figure, general_size, pal_10_q, tool_label, tools_levels, texture){
  
  values_plot = pal_10_q[tools_levels %in% unigenes$tool & tools_levels %in% tools_for_figure]
  labels_plot  = tool_label[tools_levels %in% unigenes$tool & tools_levels %in% tools_for_figure]
  unigenes <- unigenes %>% mutate(texture = ifelse(tool %in% texture, "yes", "no"))
  
  p <- unigenes  %>% filter(tool %in% tools_for_figure) %>% 
    ggplot(aes( x = tool)) +
    geom_bar_pattern(aes(fill = tool, pattern = texture), color = "black", linewidth = 0.2,
                     pattern_color = "white",
                     pattern_density = 0.1, 
                     pattern_spacing = 0.025, 
                     pattern_fill = "white",
                     pattern_key_scale_factor = 0.6) +
    scale_pattern_manual(values = c('no' = 'none', 'yes' = 'stripe')) +
    scale_fill_manual(values = values_plot, labels = labels_plot) +
    scale_x_discrete( labels =  lab_fn(labels_plot)) +
    theme_minimal() +
    ylab("Number of ARGs") +
    xlab("") +
    ggtitle("") +
    labs(fill = "") +
    scale_y_continuous(expand = c(0.01, 0.01), breaks = pretty_breaks(n = 5),
                       labels = scales::label_number()) + 
    theme(
      legend.position = "bottom",
      legend.text = element_text(size = general_size),
      panel.border = element_blank(),
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
  return(p)
}


plot_count_genes_tool_legends <- function(unigenes, tools_for_figure, general_size, pal_10_q, tool_label, tools_levels){
  
  values_plot = pal_10_q[tools_levels %in% unigenes$tool & tools_levels %in% tools_for_figure]
  labels_plot  = tool_label[tools_levels %in% unigenes$tool & tools_levels %in% tools_for_figure]
  
  p <- unigenes  %>% filter(tool %in% tools_for_figure) %>% 
    ggplot(aes( x = tool)) +
    geom_bar(aes(fill = tool), color = "black", linewidth = 0.2) +
    scale_fill_manual(values = values_plot, labels = labels_plot) +
    scale_x_discrete( labels =  labels_plot) +
    theme_minimal() +
    ylab("Number of ARGs") +
    xlab("") +
    ggtitle("") +
    labs(fill = "") +
    scale_y_continuous(expand = c(0, 0)) + 
    theme(
      legend.position = "bottom",
      legend.text = element_text(size = general_size),
      panel.border = element_blank(),
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
  #return(plotly::ggplotly(p))
  return(p)
}


summarise_abundance_diversity <- function(df, type_data){
  if(type_data == "abundance"){
    df <- df %>% 
      group_by(habitat2, location, tool, sample, aggregation) %>% 
      summarise(total = sum(normed10m)) %>%
      ungroup() %>% group_by(habitat2, location, tool, aggregation) %>%
      summarise(md = median(total) + 1e-20, q25 = quantile(total, 0.25) + 1e-20, q75 = quantile(total, 0.75) + 1e-20)
  } else{
    if(type_data == "diversity"){
      df <- df  %>% 
        group_by(habitat2, location, tool, sample, aggregation) %>% 
        summarise(total = sum(unigenes)) %>%
        ungroup() %>% group_by(habitat2, location, tool, aggregation) %>%
        summarise(md = median(total) + 1e-20, q25 = quantile(total, 0.25) + 1e-20, q75 = quantile(total, 0.75) + 1e-20)
    } else {
      print("error")
    }
  }
  return(df)
}  





plot_total_abundance_diversity_tool <- function(dataset, tools_to_plot, environments_plot, general_size, text_yaxis, pal_10_q, tool_label, tools_levels){
  n <- length(unique(dataset$habitat2))
  
  df <- data.frame(x = 1:n, y = c(1:n))
  max_y <- max(dataset$q75)
  rects <- data.frame(
    xmin = df$x - .1, xmax = df$x + .9, ymin = 0.1, ymax = max_y,
    alpha = rep( c(0, .1), length.out = nrow(df)))
  
  
  
  dataset2 <- dataset  %>% 
    filter(tool %in% tools_to_plot) %>% 
    mutate(q25 = ifelse(q75 >= 0.1, max(0.1, q25), 0.1),
           md = ifelse(q75 >= 0.1, max(0.1, md), 0.1)) %>%
    mutate(q75 = ifelse(q75 >= 0.1, q75, 0.1))
  
  labels_x_axis <- levels(dataset2$habitat2)[levels(dataset2$habitat2) %in% environments_plot]
  
  values_plot = pal_10_q[tools_levels %in% tools_to_plot & tools_levels %in% dataset2$tool ]
  labels_plot = tool_label[tools_levels %in% tools_to_plot & tools_levels %in% dataset2$tool ]
  
  p <-  dataset2  %>%
    ggplot(aes( x = habitat2_factor, y = md, fill = tool, color = tool, group = interaction(tool, location))) +
    geom_rect(data = rects, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, alpha = alpha), fill = "grey", inherit.aes = FALSE) +
    geom_line(linetype = 3, linewidth = 0.7, show.legend = FALSE) +
    geom_rect(aes(
      xmin = as.numeric(habitat2_factor) - 0.03,
      xmax = as.numeric(habitat2_factor) + 0.03,
      ymin = q25,
      ymax = q75), color = "black", linewidth = 0.2) + 
    geom_point(shape = 15, color = "black") +
    scale_color_manual(values = values_plot, labels = labels_plot) +
    scale_fill_manual(values = values_plot, labels = labels_plot) +
    theme_minimal() +
    scale_y_log10(labels = scales::label_math(),
                  limits = c(1e-1, max_y), expand = c(0, 0)) + 
    scale_x_continuous(limits = c(0.9, n + 1), expand = c(0, 0), breaks = c(1:n) + 0.5, labels = labels_x_axis) + 
    scale_alpha(range = c(0, 0.3), guide = "none") +
    ylab(text_yaxis) +
    xlab("") +
    labs(fill = "") +
    ggtitle("") +
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
      axis.text.y = element_text(size = general_size))
  
  return(p)
}


################################################################################################
################################################################################################
# Figure 2

round_up_pow10 <- function(x) {
  if (x <= 0) return(0)
  10^ceiling(log10(x))
}

vector_to_next_pow10 <- function(x) {
  if (x <= 0) return(numeric(0))  
  start_exp <- -1
  end_exp <- ceiling(log10(x))
  10^seq(start_exp, end_exp)
}

# pan resistome 
pan_resistome_plot <- function(dataset, tools_levels_fig2, environments_plot_fig2, general_size, pal_10_q, tool_label, h2, tools_levels){
  
  dataset2 <- dataset %>% 
    filter(habitat %in% environments_plot_fig2, tool %in% tools_levels_fig2) %>% 
    filter(aggregation %in% "new_level_centroid") %>% 
    mutate(log10_mn = ifelse(mn > 0, -log10(mn), - 1e-20),
           log10mn_sm_low = ifelse(mn + sd > 0, -log10(mn + sd), - 1e-20),
           log10mn_sm_top = ifelse(mn - sd > 0, -log10(mn - sd), - 1e-20))
  
  max_pan <- -10^round_up_pow10(max(dataset2$mn + dataset2$sd, na.rm = T))
  
  exponents <- round(seq(from = 0, to = max(dataset2$mn + dataset2$sd), length.out = 5), -3)
  values_plot <- pal_10_q[tools_levels %in% tools_levels_fig2 & tools_levels %in% dataset2$tool] 
  labels_plot <- tool_label[tools_levels %in% tools_levels_fig2 & tools_levels %in% dataset2$tool]
  
  pan_plot <- dataset2 %>% 
    ggplot(aes(x = -mn, y = habitat, fill = tool)) +
    geom_col(position = position_dodge2(preserve = "single", width = 1), width = 1, color = "black", linewidth = 0.2) +
    theme_minimal() +
    labs(fill = "") +
    facet_grid(habitat ~ ., scales = "free_y") + 
    ylab("") +
    xlab("Number of ARGs") +
    ggtitle("Pan-resistome") + 
    scale_x_continuous(breaks = -exponents, labels = exponents) +
    scale_fill_manual(values = values_plot, labels = labels_plot) +
    theme(legend.position = "bottom",
          panel.border = element_blank(),
          axis.text.x = element_text(angle = 0, hjust = 1, vjust = 0.5, size = general_size),
          axis.text.y = element_text(size = general_size),
          strip.text = element_blank(),
          plot.margin = margin(0, 0, 0, 0, unit = "pt"),
          legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
          legend.margin = margin(0, 0, 0, 0, unit = "pt"),
          title = element_text(size = general_size + 2, face = "bold"),
          panel.spacing = unit(0, "pt"),
          legend.text = element_text(size = general_size),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          axis.title = element_text(size = general_size + 1, face = "bold"),
          panel.background = element_rect(colour = "black", fill = NA)) 
  
  return(pan_plot)
}


g_legend <- function(a.gplot){
  tmp <- ggplotGrob(a.gplot)
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}


sum_core_adjust <- function(core, cnt_subset = 900, threshold_samples = 0.5){
  return(core %>% 
           filter(cut %in% threshold_samples & cnt > cnt_subset) %>% 
           ungroup() %>% 
           group_by(new_level, tool, habitat) %>% 
           summarise(unigenes = n_distinct(X)))
}



core_resistome_plot <- function(dataset, tools_levels_fig2, environments_plot_fig2, general_size, pal_10_q, tool_label, h2, tools_levels){
  dataset2 <- dataset  %>% ungroup() %>% 
    group_by(habitat, tool) %>% 
    summarise(unigenes = sum(unigenes)) %>% 
    filter(habitat %in% environments_plot_fig2) %>% 
    filter(tool %in% tools_levels_fig2)
  
  dataset2 <- dataset2 %>% ungroup() %>% complete( habitat, tool, 
                                                   fill = list(unigenes = 0)) %>%
    filter(habitat %in% environments_plot_fig2)
  
  values_plot_core = pal_10_q[tools_levels %in% tools_levels_fig2 & tools_levels %in% dataset2$tool] 
  labels_plot_core = tool_label[tools_levels %in% tools_levels_fig2 & tools_levels %in% dataset2$tool]
  
  core_environment_supplement <- dataset2 %>% 
    ggplot(aes(x = unigenes, y = habitat, fill = tool)) +
    geom_col(position = position_dodge2(preserve = "single", width = 1), width = 1, color = "black", linewidth = 0.2) +
    theme_minimal() +
    labs(fill = "") +
    facet_grid(habitat ~ ., scales = "free_y") +
    ylab("") +
    xlab("Number of ARGs") +
    ggtitle("Core-resistome") +
    scale_fill_manual(values = values_plot_core, labels = labels_plot_core) +
    theme(legend.position = "bottom",
          legend.text = element_text(size = general_size),
          panel.border = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          plot.margin = margin(0, 0, 0, 0, unit = "pt"),
          legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
          legend.margin = margin(0, 0, 0, 0, unit = "pt"),
          panel.spacing = unit(0, "pt"),
          title = element_text(size = general_size + 2, face = "bold"),
          strip.text = element_blank(),
          axis.title = element_text(size = general_size + 1, face = "bold"),
          axis.text.x = element_text(angle = 0, size = general_size),
          axis.text.y = element_text(size = general_size, hjust = 1.0),
          panel.background = element_rect(colour = "black", fill = NA)) 
  
  return(core_environment_supplement)
}


plot_fig2 <- function(pan_resistome_plot, core_resistome_plot){
  p <- grid.arrange(pan_resistome_plot + theme(axis.text.y = element_blank(), legend.position = "none"), 
                    core_resistome_plot + theme(legend.position = "none"),
                    g_legend(core_resistome_plot),
                    layout_matrix = rbind(c(1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2),
                                          c(1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2),
                                          c(1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2),
                                          c(1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2),
                                          c(1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2),
                                          c(1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2),
                                          c(1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2),
                                          c(1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2),
                                          c(3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3)))
  return(p)
}



##########################################################################################
##########################################################################################
## Figure 3


plot_core_env_tool_class <- function(dataset, tools_levels_fig3, environments_plot_fig3, classes_plot_fig3, general_size, pal_10_q, tool_label, tools_levels){
  
  env_core <- dataset %>% 
    filter(habitat %in% environments_plot_fig3, tool %in% tools_levels_fig3) %>% 
    ungroup() %>% 
    mutate(unigenes = as.numeric(unigenes), new_level = as.character(new_level)) %>% 
    mutate(new_level = ifelse(new_level %in% classes_plot_fig3, new_level, "Other")) %>% 
    mutate(new_level = factor(new_level, levels = c(classes_plot_fig3, "Other"))) %>% 
    group_by(habitat, tool, new_level) %>%
    summarise(partial_sum = sum(unigenes)) %>%
    ungroup() %>% 
    arrange(habitat, tool, desc(partial_sum)) %>% 
    ungroup() %>%
    complete(habitat, tool, new_level,  
             fill = list(partial_sum = 0, partial_order = 0)) %>%
    group_by(habitat, tool) %>%
    mutate(partial_order = 1:n())
  
  env_core <- droplevels(env_core)
  env_core <- env_core %>% filter(!is.na(partial_sum))
  
  env_core_env_non_zero <- env_core %>% ungroup() %>% group_by(habitat) %>% summarise(n = sum(partial_sum)) %>%
    filter(n > 0 )
  
  env_core2 <- env_core %>% filter(habitat %in% environments_plot_fig3, tool %in% tools_levels_fig3, new_level %in% c(classes_plot_fig3, "Other")) %>%
    ungroup() %>% group_by(new_level) %>%
    filter(habitat %in% env_core_env_non_zero$habitat) %>%
    ungroup() %>% 
    group_by(tool, habitat, new_level)
  
  
  p_core_env <- env_core2 %>% 
    ggplot(aes(y = partial_sum, x = new_level, fill = tool)) +
    geom_col(position = position_dodge2(preserve = "single"), color = "black", linewidth = 0.1) +
    theme_minimal() +
    labs(color = "") +
    facet_grid(habitat ~ new_level, scales = "free") + 
    xlab("") +
    ylab("") +
    ggtitle("") +
    labs(fill = "") +
    scale_fill_manual(values = pal_10_q[tools_levels %in% tools_levels_fig3 & tools_levels %in% env_core2$tool], 
                      labels = tool_label[tools_levels %in% tools_levels_fig3 & tools_levels %in% env_core2$tool]) +
    scale_y_continuous(labels = function(x) sprintf("%6d", as.integer(x))) + 
    theme(legend.position = "bottom",
          panel.border = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = general_size),
          axis.text.y = element_text(size = general_size),
          strip.text = element_text(size = general_size, face = "bold"),
          strip.text.x  = element_blank(),
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
  
  return(p_core_env)
  
}


abundance_medians_env_class <- function(abundance_class, environments_plot_fig3, tools_levels_fig3, classes_plot_fig3){
  abundance_medians_env <- abundance_class %>% ungroup() %>%
    filter(habitat %in% environments_plot_fig3, tool %in% tools_levels_fig3) %>%
    mutate(gene = as.character(gene)) %>%
    mutate(gene = ifelse(gene %in% classes_plot_fig3, gene, "Other")) %>% 
    mutate(gene = factor(gene, levels = c(classes_plot_fig3, "Other"))) %>% 
    group_by(tool, habitat, gene, sample) %>% 
    summarise(total = sum(normed10m)) %>%
    ungroup() 
  
  return(abundance_medians_env)
}


calc_boxplot_stat <- function(x) {
  coef <- 1e-20
  n <- sum(!is.na(x))
  stats <- quantile(x, probs = c(0.0, 0.25, 0.5, 0.75, 1.0))
  names(stats) <- c("ymin", "lower", "middle", "upper", "ymax")
  iqr <- diff(stats[c(2, 4)])
  outliers <- x < (stats[2] - coef * iqr) | x > (stats[4] + coef * iqr)
  if (any(outliers)) {
    stats[c(1, 5)] <- range(c(stats[2:4], x[!outliers]), na.rm = TRUE)
  }
  return(stats)
}



plot_abundance_env_tool_class <- function(dataset, tools_levels_fig3, environments_plot_fig3, classes_plot_fig3, general_size, pal_10_q, tool_label, tools_levels){
  
  dataset2 <- dataset %>% 
    ungroup() %>%
    group_by(tool, habitat, gene) %>%
    summarise(q75 = quantile(total, 0.75)) %>%
    ungroup() %>%
    filter(q75 > 0) %>%
    mutate(id = paste(tool, habitat, gene))
  
  
  n <- length(unique(dataset2$gene))
  n <- length(unique(classes_plot_fig3)) + 1
  vlines <- data.frame(xint = seq(1.5, n - 0.5, by = 1))
  
  
  abundance_class_plot_env <- dataset %>% 
    complete(habitat, tool, gene, fill = list(total = 0)) %>%
    ggplot(aes( x = gene, y = total , fill = tool)) +
    #geom_vline(data = vlines, aes(xintercept = xint), color = "black", linewidth = 0.5) +
    stat_summary(fun.data = calc_boxplot_stat, geom="boxplot", position = position_dodge2(preserve = "single"), color = "black", linewidth = 0.2) + 
    scale_fill_manual(values = pal_10_q[tools_levels %in% tools_levels_fig3 ], 
                      labels = tool_label[tools_levels %in% tools_levels_fig3]) +
    scale_y_continuous(labels = function(x) sprintf("%6d", as.integer(x))) + 
    #facet_wrap( ~ habitat, scales="free_y", ncol =  1) + 
    facet_grid( ~ gene, scales = "free") + 
    theme_minimal() +
    ylab("Relative abundance") +
    xlab("") +
    labs(fill = "") +
    ggtitle("") +
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
      strip.text = element_text(size = general_size, face = "bold"),
      title = element_text(size = general_size + 2, face = "bold"),
      axis.title = element_text(size = general_size + 1, face = "bold"),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = general_size),
      axis.text.y = element_text(size = general_size),
      panel.background = element_rect(colour = "black", fill = NA)) 
  
  return(abundance_class_plot_env)
  
}



fix_abundance_class <- function(abundance, not_env, tools_levels, metadata){
  df <- abundance %>% 
    ungroup() %>% filter(tool %in% tools_levels, !habitat %in% not_env, aggregation %in% "new_level") %>% 
    mutate(gene = factor(gene)) %>%
    mutate(sample = factor(sample)) %>%
    mutate(tool = factor(tool, levels = tools_levels)) %>%
    select(!c(raw, raw_unique, scaled, unigenes_rarefied, unigenes_raw, unigenes_raw_unique)) %>% 
    complete( sample, gene, tool, 
              fill = list(normed10m = 0,  unigenes = 0)) %>%
    mutate(aggregation = "new_level") %>%
    mutate(gene = as.character(gene))
  
  df <- df %>% mutate(habitat = metadata$habitat[match(sample, metadata$sample_id)])
  df <- df %>% mutate(habitat2 = metadata$habitat2[match(sample, metadata$sample_id)])
  df <- df %>% mutate(location = metadata$location[match(sample, metadata$sample_id)])
  df <- df %>% mutate(aggregation = "new_level")
  return(df)
} 



get_unigenes_class <- function(dataset, tools_levels_fig3, classes_plot_fig3){
  df <- dataset %>% ungroup() %>%  
    filter(tool %in% tools_levels_fig3) %>% 
    mutate(new_level = ifelse(new_level %in% classes_plot_fig3, new_level, "Other")) %>%
    mutate(new_level = factor(new_level, levels = c(classes_plot_fig3, "Other"))) %>% 
    group_by(tool, new_level) %>% summarise( n = n_distinct(query)) %>% 
    arrange(tool, desc(n))
  return(df)
}



plot_unigenes_env_tool_class <- function(dataset, general_size, pal_10_q, tool_label, tools_levels, tools_levels_fig3, min_value, max_value, n_breaks, expand_y){
  
  
  dataset2  <- dataset %>%
    ungroup() %>%
    complete(tool, new_level, fill = list(n = 0)) 
  
  n <- length(unique(dataset2$new_level)) 
  vlines <- data.frame(xint = seq(1.5, n - 0.5, by = 1))
  
  unigenes_class_plot <- dataset2 %>%
    ggplot(aes( x = new_level, y = n, fill = tool, color = tool)) +
    geom_col(position = position_dodge2(preserve = "single"), color = "black", linewidth = 0.2) +
    coord_cartesian(ylim = c(min_value, max_value), expand = expand_y) +
    scale_y_continuous(breaks = scales::breaks_pretty(n = n_breaks), labels = function(x) sprintf("%5d", as.integer(x))) + 
    facet_grid(. ~ new_level, scales  = "free") +
    scale_fill_manual(values = pal_10_q[tools_levels %in% tools_levels_fig3 & tools_levels %in% dataset$tool], 
                      labels = tool_label[tools_levels %in% tools_levels_fig3 & tools_levels %in% dataset$tool]) +
    scale_color_manual( values = pal_10_q[tools_levels %in% tools_levels_fig3 & tools_levels %in% dataset$tool], 
                        labels = tool_label[tools_levels %in% tools_levels_fig3 & tools_levels %in% dataset$tool]) +
    xlab("") +
    #ggtitle("A") + 
    ylab("Number of ARGs") +
    theme_minimal() +
    labs(fill = "", color = "") +
    theme(legend.position = "none",
          panel.border = element_blank(),
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
          axis.title = element_text(size = general_size + 1, face = "bold"),
          panel.background = element_rect(colour = "black", fill = NA))
  
  return(unigenes_class_plot)
  
}



plot_unigenes_env_tool_class_log10 <- function(dataset, general_size, pal_10_q, tool_label, tools_levels, tools_levels_fig3){
  
  
  dataset2  <- dataset %>%
    ungroup() %>%
    complete(tool, new_level, fill = list(n = 0)) 
  
  n <- length(unique(dataset2$new_level)) 
  vlines <- data.frame(xint = seq(1.5, n - 0.5, by = 1))
  
  unigenes_class_plot <- dataset2 %>%
    ggplot(aes( x = new_level, y = n+1e-20, fill = tool, color = tool)) +
    geom_col(position = position_dodge2(preserve = "single"), color = "black", linewidth = 0.2) +
    scale_y_log10(labels = scales::math_format(10^.x)(0:5), expand = c(0, 0),
                  breaks = 10^(0:5)) +
    coord_cartesian(ylim = c(1, 8e4)) +
    facet_grid(. ~ new_level, scales  = "free") +
    scale_fill_manual(values = pal_10_q[tools_levels %in% tools_levels_fig3 & tools_levels %in% dataset$tool], 
                      labels = tool_label[tools_levels %in% tools_levels_fig3 & tools_levels %in% dataset$tool]) +
    scale_color_manual( values = pal_10_q[tools_levels %in% tools_levels_fig3 & tools_levels %in% dataset$tool], 
                        labels = tool_label[tools_levels %in% tools_levels_fig3 & tools_levels %in% dataset$tool]) +
    xlab("") +
    #ggtitle("A") + 
    ylab("Number of ARGs") +
    theme_minimal() +
    labs(fill = "", color = "") +
    theme(legend.position = "none",
          panel.border = element_blank(),
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
          axis.title = element_text(size = general_size + 1, face = "bold"),
          panel.background = element_rect(colour = "black", fill = NA))
  
  return(unigenes_class_plot)
  
}

##########################################################################################
##########################################################################################
## Figure 4 


###
# overlap all genes between tools and classes
# create lists


#select_class_tool_fnr <- function(dataset, tool_fig4, class_fig4){
#  return(dataset %>% filter(!is.na(recall)) %>%
#           filter(tool_ref %in% tool_fig4, tool_comp %in% tool_fig4, new_level %in% class_fig4))
#}




plot_recall_fnr <- function(dataset, tool_fig4, class_fig4, tool_label, pal_10_q, general_size, tools_levels){
  
  filter_dataset <- dataset %>% filter(!is.na(recall)) %>%
    filter(tool_ref %in% tool_fig4, tool_comp %in% tool_fig4, new_level %in% class_fig4)
  
  recall_plot <- filter_dataset %>%
    ggplot(aes(x = tool_ref, y = recall, fill = tool_ref)) + 
    geom_boxplot( linewidth = 0.2) +
    scale_fill_manual(values = pal_10_q[tools_levels %in% tool_fig4 & tools_levels %in% filter_dataset$tool_ref], 
                      labels = tool_label[tools_levels %in% tool_fig4 & tools_levels %in% filter_dataset$tool_ref]) +
    scale_x_discrete( labels = tool_label[tools_levels %in% tool_fig4 & tools_levels %in% filter_dataset$tool_ref]) +
    theme_minimal() +
    coord_cartesian(ylim = c(0,1)) +
    ylab("CSTC") +
    xlab("") +
    labs(fill = "") +
    ggtitle("A") +
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
      panel.background = element_rect(colour = "black", fill = NA))
  
  filter_dataset2 <- dataset %>% filter(!is.na(fnr)) %>%
    filter(tool_ref %in% tool_fig4, tool_comp %in% tool_fig4, new_level %in% class_fig4)  
  
  fnr_plot <- filter_dataset2 %>% 
    ggplot(aes(x = tool_ref, y = fnr, fill = tool_ref)) + 
    geom_boxplot(linewidth = 0.2) +
    scale_fill_manual(values = pal_10_q[tools_levels %in% tool_fig4 & tools_levels %in% filter_dataset2$tool_ref], 
                      labels = tool_label[tools_levels %in% tool_fig4 & tools_levels %in% filter_dataset2$tool_ref]) +
    scale_x_discrete( labels = tool_label[tools_levels %in% tool_fig4 & tools_levels %in% filter_dataset2$tool_ref]) +
    theme_minimal() +
    coord_cartesian(ylim = c(0,1)) +
    ylab("CSNO") +
    xlab("") +
    labs(fill = "") +
    ggtitle("B") +
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
      panel.background = element_rect(colour = "black", fill = NA))
  
  return(grid.arrange(recall_plot, fnr_plot, ncol = 2))
  
}


plot_recall_fnr_texture <- function(dataset, tool_fig4, class_fig4, tool_label, pal_10_q, general_size, tools_levels, texture){
  
  dataset <- dataset %>% mutate(texture = ifelse(tool_ref %in% texture, "yes", "no"))
  
  filter_dataset <- dataset %>% filter(!is.na(recall)) %>%
    filter(tool_ref %in% tool_fig4, tool_comp %in% tool_fig4, new_level %in% class_fig4)
  
  recall_plot <- filter_dataset %>%
    ggplot(aes(x = tool_ref, y = recall, fill = tool_ref, pattern = texture)) + 
    geom_boxplot_pattern(linewidth = 0.2,
                         pattern_color = "white",
                         pattern_density = 0.1, 
                         pattern_spacing = 0.025, 
                         pattern_fill = "white",
                         pattern_key_scale_factor = 0.6) +
    scale_pattern_manual(values = c('no' = 'none', 'yes' = 'stripe')) +
    scale_fill_manual(values = pal_10_q[tools_levels %in% tool_fig4 & tools_levels %in% filter_dataset$tool_ref], 
                      labels = tool_label[tools_levels %in% tool_fig4 & tools_levels %in% filter_dataset$tool_ref]) +
    scale_x_discrete( labels = tool_label[tools_levels %in% tool_fig4 & tools_levels %in% filter_dataset$tool_ref]) +
    theme_minimal() +
    coord_cartesian(ylim = c(0,1)) +
    ylab("CSTC") +
    xlab("") +
    labs(fill = "") +
    ggtitle("A") +
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
      panel.background = element_rect(colour = "black", fill = NA))
  
  filter_dataset2 <- dataset %>% filter(!is.na(fnr)) %>%
    filter(tool_ref %in% tool_fig4, tool_comp %in% tool_fig4, new_level %in% class_fig4)  
  
  fnr_plot <- filter_dataset2 %>% 
    ggplot(aes(x = tool_ref, y = fnr, fill = tool_ref, pattern = texture)) + 
    geom_boxplot_pattern(linewidth = 0.2,
                         pattern_color = "white",
                         pattern_density = 0.1, 
                         pattern_spacing = 0.025, 
                         pattern_fill = "white",
                         pattern_key_scale_factor = 0.6) +
    scale_pattern_manual(values = c('no' = 'none', 'yes' = 'stripe')) +
    scale_fill_manual(values = pal_10_q[tools_levels %in% tool_fig4 & tools_levels %in% filter_dataset2$tool_ref], 
                      labels = tool_label[tools_levels %in% tool_fig4 & tools_levels %in% filter_dataset2$tool_ref]) +
    scale_x_discrete( labels = tool_label[tools_levels %in% tool_fig4 & tools_levels %in% filter_dataset2$tool_ref]) +
    theme_minimal() +
    coord_cartesian(ylim = c(0,1)) +
    ylab("CSNO") +
    xlab("") +
    labs(fill = "") +
    ggtitle("B") +
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
      panel.background = element_rect(colour = "black", fill = NA))
  
  return(grid.arrange(recall_plot, fnr_plot, ncol = 2))
  
}



plot_recall_detailed <- function(dataset, tool_fig4, class_fig4, tool_label, pal_10_q, general_size, tools_levels){
  
  filter_dataset <- dataset %>% filter(!is.na(recall)) %>%
    filter(tool_ref %in% tool_fig4, new_level %in% class_fig4)
  
  recall_plot <- filter_dataset %>%
    ggplot(aes(x = new_level, y = recall, fill = tool_ref)) + 
    geom_boxplot( linewidth = 0.2) +
    scale_fill_manual(values = pal_10_q[tools_levels %in% tool_fig4 & tools_levels %in% filter_dataset$tool_ref], 
                      labels = tool_label[tools_levels %in% tool_fig4 & tools_levels %in% filter_dataset$tool_ref]) +
    facet_grid(tool_ref ~ ., scales = "free") +
    theme_minimal() +
    coord_cartesian(ylim = c(0,1)) +
    ylab("CSTC") +
    xlab("") +
    labs(fill = "") +
    ggtitle("") +
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
      strip.text = element_text(size = general_size, face = "bold"),
      strip.text.x = element_blank(),
      panel.background = element_rect(colour = "black", fill = NA))
  
  
  return(recall_plot)
  
} 




plot_fnr_detailed <- function(dataset, tool_fig4, class_fig4, tool_label, pal_10_q, general_size, tools_levels){
  
  
  filter_dataset <- dataset %>% filter(!is.na(fnr)) %>%
    filter(tool_ref %in% tool_fig4, new_level %in% class_fig4)
  
  fnr_plot <- filter_dataset %>%
    ggplot(aes(x = new_level, y = fnr, fill = tool_ref)) + 
    geom_boxplot(linewidth = 0.2) +
    scale_fill_manual(values = pal_10_q[tools_levels %in% tool_fig4 & tools_levels %in% filter_dataset$tool_ref], 
                      labels = tool_label[tools_levels %in% tool_fig4 & tools_levels %in% filter_dataset$tool_ref]) +
    theme_minimal() +
    facet_grid(tool_ref ~ ., scales = "free") +
    coord_cartesian(ylim = c(0,1)) +
    ylab("CSNO") +
    xlab("") +
    labs(fill = "") +
    ggtitle("") +
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
      strip.text = element_text(size = general_size, face = "bold"),
      strip.text.x = element_blank(),
      panel.background = element_rect(colour = "black", fill = NA))
  
  return(fnr_plot)
}




#### ##################
one_tool_recall_fnr_id <- function(unigenes, JI_class_other_filter, tool1, pal_10_q, general_size, tool_label, tools_levels) { 
  data_plot_recall_fnr_id <- bind_rows(
    unigenes %>% 
      filter(tool %in% c(tool1)) %>%
      mutate(id = id/100, dt = "unigenes") %>%
      ungroup() %>% group_by(tool, new_level) %>% 
      summarise(md = median(id, na.rm = T), q25 = quantile(id , 0.25, na.rm = T), q75 = quantile(id , 0.75, na.rm = T)) %>% 
      mutate(dt = "Identity"),
    JI_class_other_filter %>% 
      filter(tool_ref %in% c(tool1)) %>%
      select(tool_ref, recall, new_level) %>% 
      ungroup() %>% group_by(tool_ref, new_level) %>% 
      summarise(md = median(recall, na.rm = T), q25 = quantile(recall , 0.25, na.rm = T), q75 = quantile(recall , 0.75, na.rm = T)) %>% 
      rename(tool = tool_ref) %>%
      mutate(dt  = "CSTC"),
    JI_class_other_filter %>% 
      filter(tool_ref %in% c(tool1)) %>%
      select(tool_ref, fnr, new_level) %>% 
      ungroup() %>% group_by(tool_ref, new_level) %>% 
      summarise(md = median(fnr, na.rm = T), q25 = quantile(fnr , 0.25, na.rm = T), q75 = quantile(fnr , 0.75, na.rm = T)) %>% 
      rename(tool = tool_ref) %>%
      mutate(dt  = "CSNO"))
  
  all_id_recall_fnr <- data_plot_recall_fnr_id   %>% 
    filter(tool %in% tool1) %>% 
    mutate(dt = factor(dt, levels = c("Identity", "CSTC", "CSNO"))) %>% 
    ungroup() %>% group_by(dt, tool) %>%
    arrange(dt, md) %>% 
    mutate(cl1 = factor(new_level, levels = unique(new_level))) %>%
    mutate(cl = as.numeric(cl1)) %>%
    mutate(cl = ifelse(dt %in% c("Identity"), cl - 0.15, cl)) %>%
    mutate(cl = ifelse(dt %in% c("CSNO"), cl + 0.15, cl))
  
  
  p <- all_id_recall_fnr %>% 
    ggplot(aes( y = cl, x = md, fill = dt, color = dt)) +
    geom_rect(aes(
      ymin = as.numeric(cl) - 0.09,
      ymax = as.numeric(cl) + 0.09,
      xmin = q25,
      xmax = q75), color = "black", linewidth = 0.1) + 
    geom_point( aes(color = dt), alpha = 0.75, size  = 3) +
    scale_color_manual(values = pal_10_q) +
    scale_fill_manual(values = pal_10_q) +
    theme_minimal() +
    scale_y_continuous(limits = c(0, length(levels(all_id_recall_fnr$cl1)) + 0.4), 
                       expand = c(0, 0), 
                       breaks = 1:length(levels(all_id_recall_fnr$cl1)), 
                       labels = levels(all_id_recall_fnr$cl1)) + 
    facet_grid(. ~ tool, scales = "free_y") +
    ylab("") +
    xlab("") +
    labs(fill = "", color = "") +
    ggtitle("A") +
    theme(
      legend.position = "bottom",
      legend.text = element_text(size = general_size ),
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      plot.margin = margin(0, 0, 0, 0, unit = "pt"),
      strip.text = element_blank(),
      legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
      legend.margin = margin(0, 0, 0, 0, unit = "pt"),
      panel.spacing = unit(0, "pt"),
      title = element_text(size = general_size + 2, face = "bold"),
      axis.title = element_text(size = general_size + 1, face = "bold"),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = general_size),
      axis.text.y = element_text(size = general_size),
      panel.background = element_rect(colour = "black", fill = NA))
  return(p)
  
}

merge_recall_fnr_id <- function(unigenes, JI_class_other_filter, tool1, tool2, pal_10_q, general_size, tool_label, tools_levels){
  
  data_plot_recall_fnr_id <- bind_rows(
    unigenes %>% 
      filter(tool %in% c(tool1, tool2)) %>%
      mutate(id = id/100, dt = "unigenes") %>%
      ungroup() %>% group_by(tool, new_level) %>% 
      summarise(md = median(id, na.rm = T), q25 = quantile(id , 0.25, na.rm = T), q75 = quantile(id , 0.75, na.rm = T)) %>% 
      mutate(dt = "Identity"),
    JI_class_other_filter %>% 
      filter(tool_ref %in% c(tool1, tool2)) %>%
      select(tool_ref, recall, new_level) %>% 
      ungroup() %>% group_by(tool_ref, new_level) %>% 
      summarise(md = median(recall, na.rm = T), q25 = quantile(recall , 0.25, na.rm = T), q75 = quantile(recall , 0.75, na.rm = T)) %>% 
      rename(tool = tool_ref) %>%
      mutate(dt  = "CSTC"),
    JI_class_other_filter %>% 
      filter(tool_ref %in% c(tool1, tool2)) %>%
      select(tool_ref, fnr, new_level) %>% 
      ungroup() %>% group_by(tool_ref, new_level) %>% 
      summarise(md = median(fnr, na.rm = T), q25 = quantile(fnr , 0.25, na.rm = T), q75 = quantile(fnr , 0.75, na.rm = T)) %>% 
      rename(tool = tool_ref) %>%
      mutate(dt  = "CSNO"))
  
  all_id_recall_fnr <- data_plot_recall_fnr_id   %>% 
    filter(tool %in% tool1) %>% 
    mutate(dt = factor(dt, levels = c("Identity", "CSTC", "CSNO"))) %>% 
    ungroup() %>% group_by(dt, tool) %>%
    arrange(dt, md) %>% 
    mutate(cl1 = factor(new_level, levels = unique(new_level))) %>%
    mutate(cl = as.numeric(cl1)) %>%
    mutate(cl = ifelse(dt %in% c("Identity"), cl - 0.15, cl)) %>%
    mutate(cl = ifelse(dt %in% c("CSNO"), cl + 0.15, cl))
  
  all_id_recall_fnr_2 <- data_plot_recall_fnr_id   %>% 
    filter(tool %in% tool2) %>% 
    mutate(dt = factor(dt, levels = c("Identity", "CSTC", "CSNO"))) %>% 
    ungroup() %>% group_by(dt, tool) %>%
    arrange(dt, md) %>% 
    mutate(cl1 = factor(new_level, levels = unique(new_level))) %>%
    mutate(cl = as.numeric(cl1)) %>%
    mutate(cl = ifelse(dt %in% c("Identity"), cl - 0.15, cl)) %>%
    mutate(cl = ifelse(dt %in% c("CSNO"), cl + 0.15, cl))
  
  p <- all_id_recall_fnr %>% 
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
    scale_y_continuous(limits = c(0, length(levels(all_id_recall_fnr$cl1)) + 0.4), 
                       expand = c(0, 0), 
                       breaks = 1:length(levels(all_id_recall_fnr$cl1)), 
                       labels = levels(all_id_recall_fnr$cl1)) + 
    facet_grid(. ~ tool, scales = "free_y") +
    ylab("") +
    xlab("") +
    labs(fill = "", color = "") +
    ggtitle("") +
    theme(
      legend.position = "bottom",
      legend.text = element_text(size = general_size ),
      panel.border = element_blank(),
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
      axis.text.y = element_text(size = general_size),
      panel.background = element_rect(colour = "black", fill = NA))
  
  
  p2 <- all_id_recall_fnr_2 %>% 
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
    scale_y_continuous(limits = c(0, length(levels(all_id_recall_fnr_2$cl1)) + 0.4), 
                       expand = c(0, 0), 
                       breaks = 1:length(levels(all_id_recall_fnr_2$cl1)), 
                       labels = levels(all_id_recall_fnr_2$cl1)) + 
    facet_grid(. ~ tool, scales = "free_y") +
    ylab("") +
    xlab("") +
    labs(fill = "", color = "") +
    ggtitle("") +
    theme(
      legend.position = "bottom",
      legend.text = element_text(size = general_size ),
      panel.border = element_blank(),
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
      axis.text.y = element_text(size = general_size),
      panel.background = element_rect(colour = "black", fill = NA))
  
  return(grid.arrange(p, p2, nrow = 1))
}




plot_id_levels <- function(unigenes, tools_levels, pal_10_q, tool_label, general_size) {
  
  df0 <- unigenes %>% ungroup() %>% filter(tool %in% tools_levels) %>% 
    mutate(bin = cut_width(round(id), width = 2, boundary = min(round(id), na.rm = T))) %>%
    group_by(tool, bin) %>%
    summarise(n = n(), .groups = "drop_last") %>%
    mutate(frac = n / sum(n)) %>%
    ungroup() 
  
  df_min <- unigenes %>% group_by(tool) %>% summarise(mi = min(id, na.rm = T), ma = max(id, na.rm = T))
  
  df0 <- df0 %>% left_join(df_min)
  
  df <- df0 %>%
    ungroup() %>%
    mutate(id2 =  as.numeric(sub(".*,(.*)\\]", "\\1", bin))) %>%
    mutate(id1 =  as.numeric(sub("^[\\(\\[]([0-9]+),.*\\]$", "\\1", bin))) %>% 
    group_by(tool) %>%
    mutate(xmin = ifelse(id1 < mi, mi, id1)) %>%
    mutate(xmax = ifelse(id2 > ma, ma, id2)) %>%
    mutate(frac0 = frac) %>% 
    mutate(frac = ifelse(tool %in% c("ResFinder", "ABRicate-ResFinder"), 0.6 + frac,
                         ifelse(tool %in% c("AMRFinderPlus", "ABRicate-NCBI"), frac + 0.45,
                                ifelse(tool %in% c("fARGene", "ABRicate-ARGANNOT"), frac + 0.3,
                                       ifelse(tool %in% c("DeepARG", "ABRicate-MEGARes"), frac + 0.15, frac ))))) %>%
    mutate(miny = ifelse(tool %in% c("ResFinder", "ABRicate-ResFinder"), 0.6,
                         ifelse(tool %in% c("AMRFinderPlus", "ABRicate-NCBI"),  0.45,
                                ifelse(tool %in% c("fARGene", "ABRicate-ARGANNOT"),  0.3,
                                       ifelse(tool %in% c("DeepARG", "ABRicate-MEGARes"),  0.15,  0)))))
  

  p <-  df %>%
    ggplot(aes(x = id2 , color = tool, y = frac, fill = tool)) +
    geom_rect(aes(xmin = xmin, xmax = xmax, 
                  ymin = miny, ymax = frac), linewidth = 0) + 
    theme_minimal() + 
    scale_fill_manual(values = pal_10_q[tools_levels %in% df$tool], labels = tool_label[tools_levels %in% df$tool], guide = guide_legend(nrow = 2)) +
    scale_color_manual(values = pal_10_q[tools_levels %in% df$tool], labels = tool_label[tools_levels %in% df$tool], guide = guide_legend(nrow = 2)) +
    scale_x_continuous(limits = c(20, 100), breaks = seq(from = 0, to = 100, by = 10)) + 
    ggtitle("") + 
    ylab("Proportion of genes") + 
    xlab("% Identity level") + 
    labs(color = "", fill = "") +
    theme(legend.position = "bottom",
          panel.border = element_blank(),
          axis.text.x = element_text(angle = 0, hjust = 1, vjust = 0.5, size = general_size),
          axis.text.y = element_blank(),
          strip.text = element_blank(),
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
  
  return(p)
  
}




return_heatmap_overalp  <- function(JI_all, tools_selected, general_size, tool_label, tools_levels) {
  df <- JI_all %>% 
    filter(tool_ref %in% tools_selected, tool_comp %in% tools_selected)
  
  labels_plot  = tool_label[tools_levels %in% df$tool_ref & tools_levels %in% df$tool_comp  & tools_levels %in% tools_selected]
  
  heatmap_proportion_unigenes_found_in_other <- df %>% 
    ggplot(aes(x = tool_comp, y = tool_ref, fill = recall)) +
    geom_tile() +
    scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(9, "YlOrBr"))) +
    scale_x_discrete(labels = labels_plot) +
    scale_y_discrete(labels = labels_plot) +
    theme_minimal() +
    ggtitle("") + 
    labs(fill = "") +
    xlab("") +
    theme(axis.title.y = element_blank(),
          legend.position = "right",
          panel.border = element_blank(),
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
          axis.title = element_text(size = general_size + 1, face = "bold"),
          panel.background = element_rect(colour = "black", fill = NA))
  
  return(heatmap_proportion_unigenes_found_in_other)
}







plot_total_abundance_diversity_new_version <- function(
    dataset, tools_labels, tools_to_plot, 
    environments_plot, general_size,  pal_10_q, 
    metric = "abundance", sd = 2025, obs = 300, texture,
    pattern_density = 0.005,
    pattern_spacing = 0.025,
    pattern_fill = "white",
    pattern_size = 0.3){
  set.seed(sd)
  dataset2 <- dataset %>% filter(tool %in% tools_to_plot, habitat %in% environments_plot) %>%
    mutate(texture = ifelse(tool %in% texture, "yes", "no"))
  
  labels_plot  = tools_labels[tools_to_plot %in% dataset2$tool]
  values_plot = pal_10_q[tools_to_plot %in%  dataset2$tool]
  
  if(metric == "abundance"){
    
    dataset_abundance <- dataset2 %>% filter(normed10m > 0) %>%
      ungroup() %>%
      select(sample, habitat) %>%
      distinct() %>%
      arrange(habitat, sample) %>% 
      group_by(habitat) %>% 
      group_modify(~ {
        n_rows <- nrow(.x)
        .x[sample(n_rows, min(n_rows, obs)), ]
      }) %>% 
      ungroup()
    
    dataset2 <- dataset2 %>% mutate(p = NA)
    dataset2 <- dataset2 %>% mutate(p = ifelse(paste(habitat, sample) %in% paste(dataset_abundance$habitat, dataset_abundance$sample), normed10m, NA))
    
    
    
    p <- dataset2 %>%
      ggplot(aes(x = tool, fill = tool)) +
      geom_boxplot_pattern(aes(y = normed10m + 1e-20, pattern = texture), position = position_dodge2(preserve = "single"),  
                           linewidth = 0.2, outlier.shape = NA, color = "black",
        pattern_color = "white",
        pattern_density = pattern_density, 
        pattern_spacing = pattern_spacing, 
        pattern_fill = pattern_fill,
        pattern_size = pattern_size,
        pattern_key_scale_factor = 0.6) +
      scale_pattern_manual(values = c('no' = 'none', 'yes' = 'stripe')) +
      geom_jitter(aes(y = p+1e-20), alpha =0.15, size = 0.3, width = 0.25, height = 0, color = "black") +
      scale_fill_manual(values = values_plot, labels = labels_plot) +
      facet_wrap(~ habitat, scales = "free_x", nrow = 2) +
      scale_y_log10(labels = scales::math_format(10^.x)(0:5), expand = c(0, 0),
                    breaks = 10^(0:5)) +
      coord_cartesian(ylim = c(9e-1, 5e4)) +
      xlab("") +
      labs(fill = "") +
      ylab("Relative abundance") +
      ggtitle("") +
      theme_minimal() +
      theme(
        legend.position = "bottom",
        strip.text = element_text(size = general_size),
        legend.text = element_text(size = general_size ),
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
  } else{
    if( metric == "diversity"){
      dataset_abundance <- dataset2 %>% filter(unigenes > 0) %>%
        ungroup() %>%
        select(sample, habitat) %>%
        distinct() %>%
        arrange(habitat, sample) %>% 
        group_by(habitat) %>% 
        group_modify(~ {
          n_rows <- nrow(.x)
          .x[sample(n_rows, min(n_rows, obs)), ]
        }) %>% 
        ungroup()
      
      dataset2 <- dataset2 %>% mutate(p = NA)
      dataset2 <- dataset2 %>% mutate(p = ifelse(paste(habitat, sample) %in% paste(dataset_abundance$habitat, dataset_abundance$sample), unigenes, NA))
      

      p <- dataset2 %>%
        ggplot(aes(x = tool, fill = tool)) +
        geom_boxplot_pattern(aes(y = unigenes + 1e-20, pattern = texture), position = position_dodge2(preserve = "single"),  
                     linewidth = 0.2, outlier.shape = NA, color = "black", 
                     pattern_color = "white",
                     pattern_density = pattern_density, 
                     pattern_spacing = pattern_spacing, 
                     pattern_fill = pattern_fill,
                     pattern_size = pattern_size,
                     pattern_key_scale_factor = 0.6) +
        scale_pattern_manual(values = c('no' = 'none', 'yes' = 'stripe')) +
        geom_jitter(aes(y = p+1e-20), alpha =0.15, size = 0.3, width = 0.25, height = 0, color = "black") +
        scale_fill_manual(values = values_plot, labels = labels_plot) +
        facet_wrap(~ habitat, scales = "free_x", nrow = 2) +
        scale_y_log10(labels = scales::math_format(10^.x)(0:4), expand = c(0, 0),
                      breaks = 10^(0:4)) +
        coord_cartesian(ylim = c(9e-1, 5e3)) +
        xlab("") +
        ylab("Diversity") +
        labs(fill = "") +
        ggtitle("") +
        theme_minimal() +
        theme(
          legend.position = "bottom",
          strip.text = element_text(size = general_size),
          legend.text = element_text(size = general_size ),
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
    } else{
      p <- NA
    }
  }
  return(p)
}



plot_total_abundance_diversity_new_version_both_together <- function(dataset, tools_to_plot, environments_plot, general_size, text_yaxis, pal_10_q, tool_label, tools_levels){
  n <- length(unique(dataset$habitat2))
  
  df <- data.frame(x = 1:n, y = c(1:n))
  max_y <- max(dataset$q75[dataset$metric == "abundance"])
  max_x <- max(dataset$q75[dataset$metric == "diversity"])
  rects <- data.frame(
    xmin = df$x - .1, xmax = df$x + .9, ymin = 0.1, ymax = max_y,
    alpha = rep( c(0, .1), length.out = nrow(df)))
  
  
  
  dataset2 <- dataset  %>% ungroup() %>%
    filter(tool %in% tools_to_plot) %>% 
    group_by(tool, habitat2, location, aggregation, metric) %>% 
    mutate(q25 = ifelse(q75 >= 0.1, max(0.1, q25), 0.1),
           md = ifelse(q75 >= 0.1, max(0.1, md), 0.1)) %>%
    mutate(q75 = ifelse(q75 >= 0.1, q75, 0.1))
  
  dataset2 <- dataset2 %>% pivot_wider(names_from = metric, values_from = c(md, q25, q75))

  values_plot = pal_10_q[tools_levels %in% tools_to_plot & tools_levels %in% dataset2$tool ]
  labels_plot = tool_label[tools_levels %in% tools_to_plot & tools_levels %in% dataset2$tool ]
  
  p <-  dataset2  %>%
    ggplot(aes( x = q75_diversity, y = q75_abundance, fill = tool, color = tool)) +
    geom_rect(aes(
      xmin = md_diversity - log10(md_diversity) ,
      xmax = md_diversity + log10(md_diversity),
      ymin = q25_abundance,
      ymax = q75_abundance), color = "black", linewidth = 0.2) + 
    geom_rect(aes(
      ymin = md_abundance - log10(md_abundance),
      ymax = md_abundance + log10(md_abundance),
      xmin = q25_diversity,
      xmax = q75_diversity), color = "black", linewidth = 0.2) + 
    scale_color_manual(values = values_plot, labels = labels_plot) +
    scale_fill_manual(values = values_plot, labels = labels_plot) +
    theme_minimal() +
    scale_x_log10(labels = scales::label_math(),
                  limits = c(1e-4, max(dataset2$q75_diversity))) + 
    scale_y_log10(labels = scales::label_math(),
                  limits = c(1e-4, max(dataset2$q75_abundance))) + 
    facet_wrap(~ habitat2, scales = "free") +
    xlab("") +
    labs(fill = "") +
    ggtitle("") +
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
      axis.text.x = element_text(hjust = 1, vjust = 0.5, size = general_size),
      axis.text.y = element_text(size = general_size))
  
  return(p)
}


dot_plot_pan_core <- function(pan_core, pan_core_env, tool_label, tools_to_plot, tools_levels, tools_texture, general_size, pal_10_q){
  
  df <- pan_core %>% 
    filter(habitat %in% pan_core_env, tool %in% tools_to_plot) 
  labels_plot  = tool_label[tools_to_plot %in% tools_levels]
  values_plot = pal_10_q[tools_to_plot %in%  tools_levels]
  
  p <- df %>% 
    ggplot(aes(x = md, y = core +1, fill = tool)) +
    geom_point(aes(shape = ifelse(tool %in% tools_texture, "yes", "no")), size = 3, color = "black") +
    scale_fill_manual(values = values_plot, labels = labels_plot) +
    scale_shape_manual(values = c(21,24)) +
    facet_wrap( ~ habitat, scales = "free", nrow = 2)  +
    theme_minimal() +
    ylab("Core-resistome") +
    xlab("Pan-resistome") +
    ggtitle("") +
    labs(fill = "") +
    theme(
      legend.position = "bottom",
      legend.text = element_text(size = general_size),
      plot.margin = margin(0, 0, 0, 0, unit = "pt"),
      legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
      legend.margin = margin(0, 0, 0, 0, unit = "pt"),
      panel.spacing = unit(0, "pt"),
      title = element_text(size = general_size + 2, face = "bold"),
      axis.title = element_text(size = general_size + 1, face = "bold"),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = general_size),
      axis.text.y = element_text(size = general_size),
      panel.background = element_rect(colour = "black", fill = NA))
  
  return(p)
}



plot_pan_core_class <- function(pan, sumcore, top20, tools_labels, tools_to_plot, tools_levels, environments_plot, general_size, pal_10_q, tools_texture){
  
  labels_plot  = tools_labels[tools_to_plot %in% tools_levels]
  values_plot = pal_10_q[tools_to_plot %in%  tools_levels]
  
  sum_pan_environment <- pan %>% 
    ungroup() %>% 
    filter(habitat %in% environments_plot, 
           tool %in% tools_to_plot) %>% 
    group_by(tool, habitat, aggregation, epoch, gene_class) %>% # sum by epoch
    summarise(s = sum(unigenes)) %>%
    ungroup() %>% 
    group_by(tool, habitat, aggregation, gene_class) %>% 
    summarise(md = median(s), mn = mean(s), sd = sd(s)) %>%  # stats by class and tool
    ungroup() %>% 
    select(tool, habitat, gene_class, md) %>% 
    rename(new_level = gene_class, unigenes = md) %>%
    ungroup() %>% 
    group_by(tool) %>% 
    mutate(N = sum(unigenes)) %>%  # total pan resistome by tool
    mutate(new_level = ifelse(new_level %in% top20, new_level, "Other")) %>% # change factor of classes
    mutate(new_level = factor(new_level, levels = c(top20, "Other"))) %>%
    ungroup() %>% 
    group_by(tool, new_level, habitat) %>% 
    summarise(unigenes = sum(unigenes)) 
  
  sum_pan_environment <- sum_pan_environment %>% 
    mutate(habitat = factor(as.character(habitat), 
                            levels = levels(habitat)[levels(habitat) %in% as.character(habitat)])) %>%
    ungroup() %>%
    complete(habitat, tool, new_level,  
             fill = list(unigenes = NA))
  
  sumcore_environment <- sumcore %>% 
    filter(habitat %in% environments_plot,
           tool %in% tools_to_plot) %>% 
    ungroup() %>% 
    group_by(tool) %>% 
    mutate(N = sum(unigenes)) %>% 
    mutate(new_level = ifelse(new_level %in% top20, new_level, "Other")) %>% 
    mutate(new_level = factor(new_level, levels = c(top20, "Other"))) %>% 
    ungroup() %>% 
    group_by(tool, new_level, habitat) %>% 
    summarise(unigenes = sum(unigenes)) 
  
  sumcore_environment <- sumcore_environment %>% 
    mutate(habitat = factor(as.character(habitat), 
                            levels = levels(habitat)[levels(habitat) %in% as.character(habitat)])) %>% 
    ungroup()  %>% 
    complete(habitat, tool, new_level,  
             fill = list(unigenes = NA))
  
  p_environment_core <- sumcore_environment %>% 
    ggplot(aes(x = new_level, y = unigenes, fill = tool)) +
    geom_jitter(aes(shape = ifelse(tool %in% tools_texture, "yes", "no")), 
                size = 3, color = "black", width = 0.4, height = 0.4, na.rm = TRUE) +
    scale_fill_manual(values = values_plot, labels = labels_plot) +
    scale_shape_manual(values = c(21,24)) +
    scale_y_continuous(expand = c(0.1, 0), breaks = scales::pretty_breaks(n = 5)) + 
    theme_minimal() +
    facet_grid(habitat ~ new_level, scales = "free_x", drop = FALSE) + 
    ylab("Core-resistome") +
    xlab("") +
    ggtitle("") +
    labs(fill = "") +
    theme(
      legend.position = "none",
      legend.text = element_text(size = general_size),
      plot.margin = margin(0, 0, 0, 0, unit = "pt"),
      legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
      legend.margin = margin(0, 0, 0, 0, unit = "pt"),
      panel.spacing = unit(0, "pt"),
      title = element_text(size = general_size + 2, face = "bold"),
      axis.title = element_text(size = general_size + 1, face = "bold"),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = general_size),
      axis.text.y = element_text(size = general_size),
      strip.text = element_blank(),
      panel.background = element_rect(colour = "black", fill = NA),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank())
  
  p_environment_pan <- sum_pan_environment %>% 
    ggplot(aes(x = new_level , y = unigenes, fill = tool)) +
    geom_jitter(aes(shape = ifelse(tool %in% tools_texture, "yes", "no")), 
                size = 3, color = "black", width = 0.4, height = 0.4, na.rm = TRUE) +
    scale_fill_manual(values = values_plot, labels = labels_plot) +
    scale_shape_manual(values = c(21,24)) +
    scale_y_log10(labels = scales::math_format(10^.x)(0:5), expand = c(.1, .1),
                  breaks = 10^(0:5)) +
    facet_grid(habitat ~ new_level, scales = "free_x", drop = FALSE) + 
    theme_minimal() +
    ylab("Pan-resistome") +
    xlab("") +
    ggtitle("") +
    labs(fill = "") +
    theme(
      legend.position = "none",
      legend.text = element_text(size = general_size),
      plot.margin = margin(0, 0, 0, 0, unit = "pt"),
      legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
      legend.margin = margin(0, 0, 0, 0, unit = "pt"),
      panel.spacing = unit(0, "pt"),
      title = element_text(size = general_size + 2, face = "bold"),
      axis.title = element_text(size = general_size + 1, face = "bold"),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = general_size),
      axis.text.y = element_text(size = general_size),
      strip.text = element_blank(),
      panel.background = element_rect(colour = "black", fill = NA),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank())
  
  return(list(p_environment_core, p_environment_pan))
  
}






plot_alluvial_classes <- function(unigenes = unigenes, 
                                  levels_unigenes = levels_unigenes, 
                                  threshold_plot = 0.99,
                                  remove_class_threshold = 0.01,
                                  tools_to_plot = tools_levels, 
                                  tools_labels = tools_labels, 
                                  tools_factors = tools_levels, 
                                  pal_10_q = pal_10_q, general_size = general_size, 
                                  gene_classes_list = gene_classes_list){

  unigenes_class <- get_unigenes_class(unigenes, tools_to_plot,  unique(unigenes$new_level)) %>% 
    ungroup() %>%
    filter(tool %in% tools_to_plot) %>%
    group_by(tool) %>% 
    mutate(proportion = n / sum(n)) %>% 
    arrange(tool, desc(proportion)) %>% 
    mutate(cum_p = cumsum(proportion)) %>% 
    ungroup() %>%
    mutate( new_level = factor(as.character(new_level), 
                               levels = rev(levels_unigenes))) %>% 
    mutate(tool = factor(as.character(tool), levels = tools_factors[tools_factors %in% tools_to_plot]))
  
  unigenes_class <- unigenes_class %>%
    group_by(tool) %>%
    arrange(proportion) %>% 
    { 
      if (any(.$cum_p == threshold_plot)) {
        filter(., cum_p <= threshold_plot)
      } else {
        lower_part <- filter(., cum_p < threshold_plot)
        upper_part <- filter(., cum_p > threshold_plot) %>% slice_tail(n = 1)
        bind_rows(lower_part, upper_part)
      }
    } %>% 
    ungroup() %>% 
    arrange(tool, cum_p)
  
  
  unigenes_class <- unigenes_class %>% 
    filter(proportion > remove_class_threshold, !new_level %in% "Other") %>% 
    ungroup() %>% 
    complete(new_level, tool, 
             fill = list(n = 0,  proportion = 0.000, cum_p = 0.000))  %>%  
    mutate(gene_name = gene_classes_list[as.character(new_level)]) %>%  
    mutate(gene_name = factor(gene_name, 
                              levels = gene_classes_list[levels(unigenes_class$new_level)]))
  
  labels_plot  <- tools_labels[tools_factors %in% tools_to_plot]
  
    
  p_alluvial <- ggplot(unigenes_class,
                       aes(x = tool,
                           stratum = gene_name,
                           alluvium = gene_name,
                           y = proportion,
                           fill = gene_name,
                           label = gene_name)) +
    geom_flow(alpha = 0.5) +
    xlab("") + 
    geom_stratum(color = "black") +
    geom_text(data = unigenes_class[unigenes_class$proportion > remove_class_threshold, ],  
              stat = "stratum",
              size = 2,color = "black",hjust = 0.5,
              position = position_jitter(width = 0, height = 0)) +
    scale_y_continuous(expand = c(.01, .01), 
                       name = "Proportion", 
                       limits = c(0, 1), 
                       breaks = pretty_breaks(n = 5)) +
    scale_fill_manual(values = c(rep(pal_10_complete,10)))+
    scale_x_discrete( labels =  lab_fn(labels_plot)) +
    theme_minimal() +
    theme(
      legend.position = "none",
      legend.text = element_text(size = general_size),
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
      axis.text.y = element_text(size = general_size))
  
  return(p_alluvial)
  
}


alluvial_pan_core_env <- function(sumcore, pan, h, tools_included, pal_10_complete, general_size){
  pan_core_env <- bind_rows(
    sumcore %>% filter(habitat %in% h) %>% ungroup() %>% mutate(dt = "core"),
    pan %>% 
      filter(habitat %in% h, aggregation %in% "new_level_centroid") %>% 
      ungroup() %>% 
      group_by(tool, habitat, epoch, gene_class) %>% 
      summarise(s = sum(unigenes)) %>%
      ungroup() %>% 
      group_by(tool, habitat, gene_class) %>% 
      summarise(unigenes = ceiling(mean(s))) %>% 
      rename(new_level  = gene_class) %>% mutate(dt = "pan")) %>%
    mutate(dt = factor(dt, levels = c("pan", "core")))
  
  pan_core_env <- pan_core_env %>% ungroup() %>%  mutate(new_level = as.character(new_level)) %>% 
    group_by(tool, new_level, dt) %>% summarise( n = sum(unigenes)) %>% 
    arrange(tool, desc(n))
  
  pan_core_env <- pan_core_env %>% 
    ungroup() %>%
    group_by(tool, dt) %>% 
    mutate(proportion = n / sum(n)) %>% 
    arrange(tool, desc(proportion)) %>% 
    mutate(cum_p = cumsum(proportion)) %>% 
    ungroup() %>%
    mutate( new_level = factor(as.character(new_level), 
                               levels = c("Other", rev(levels_unigenes))))
  
  
  pan_core_env <- pan_core_env %>%
    group_by(tool, dt) %>%
    arrange(proportion) %>% 
    { 
      if (any(.$cum_p == 0.95)) {
        filter(., cum_p <= 0.95)
      } else {
        lower_part <- filter(., cum_p < 0.99999)
        upper_part <- filter(., cum_p > 0.99999) %>% slice_tail(n = 1)
        bind_rows(lower_part, upper_part)
      }
    } %>% 
    ungroup() %>% 
    arrange(tool, cum_p)
  
  
  pan_core_env <- pan_core_env %>% 
    filter(proportion > 0.01) %>% 
    ungroup() %>% 
    complete(new_level, tool, dt,
             fill = list(n = 0,  proportion = 0.000, cum_p = 0.000))  %>%  
    mutate(gene_name = gene_classes_list[as.character(new_level)]) %>%  
    mutate(gene_name = ifelse(new_level %in% "Other", "Other", gene_name)) %>% 
    mutate(gene_name = factor(gene_name, 
                              levels = c("Other",gene_classes_list[levels(pan_core_env$new_level)])))
  
  
  pan_core_env <- pan_core_env %>%  
    filter(tool %in% tools_included) %>%  
    mutate(tool = factor(as.character(tool), levels = tools_included)) 
  
  p_alluvial_pan_core_env <- pan_core_env %>% 
    ggplot(
      aes(x = dt,
          stratum = gene_name,
          alluvium = gene_name,
          y = proportion,
          fill = gene_name,
          label = gene_name)) +
    geom_flow(alpha = 0.5) +
    xlab("") + 
    geom_stratum(color = "black") +
    geom_text(data = pan_core_env[pan_core_env$proportion > 0.01 , ],  
              stat = "stratum",
              size = 2,color = "black",hjust = 0.5,
              position = position_jitter(width = 0, height = 0)) +
    scale_y_continuous(expand = c(.01, .01), 
                       name = "Proportion", 
                       limits = c(0, 1), 
                       breaks = seq(from = 0, to = 1, by = 0.2)) +
    facet_wrap(tool ~ ., nrow = 2) +
    scale_fill_manual(values = c(rep(pal_10_complete,10))) +
    xlab("Pan- vs core-resistome") + 
    #scale_x_discrete( labels =  labels_plot) +
    theme_minimal() +
    theme(
      legend.position = "none",
      legend.text = element_text(size = general_size),
      panel.border = element_blank(),
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
      strip.text = element_text(size = general_size, face = "bold"),
      panel.background = element_rect(colour = "black", fill = NA)) +
    ggtitle(h)
  
  return(p_alluvial_pan_core_env)
}


plot_overlaps <- function(x, y, nxy, nx, ny){
  db_query <- data.frame(gene = unique(c(x, y)))
  db_query <- db_query %>% mutate(dataset = ifelse(gene %in% intersect(x, y), nxy,
                                                   ifelse(gene %in% x, nx, ny))) %>% 
    mutate(dataset = factor(dataset, levels = c(nxy, nx, ny)))
  
  
  db_query_plot <- db_query %>% 
    group_by(dataset) %>% summarise(n = n()) %>% mutate(p = n/sum(n)) %>%
    ggplot(aes(x = dataset, y = n, fill = dataset)) + 
    geom_col() + 
    scale_fill_manual(values = pal_10_q[3:10])  +
    geom_text(aes(label = scales::percent(p, accuracy = 0.1)), vjust = -0.3, size = general_size / .pt) +
    ylab("ARGs") + 
    xlab("") +
    scale_x_discrete(labels = function(x) {
      x <- gsub("-", "-\n", x)
      x <- gsub(" ", "\n", x)
      x}) +
    theme_overlap
  
  return(db_query_plot)
}


plot_abundance_class_environment <- function(abundance_class, e, general_size, pal_10_q){
  
  df_abundance_class_env <- abundance_class %>% 
    filter(habitat %in% e) %>% 
    mutate(total = normed10m) %>%
    mutate(gene = ifelse(gene %in% top20, gene, "Other")) %>% 
    mutate(texture = ifelse(tool %in% tools_texture, "yes", "no")) %>%
    mutate(gene = factor(gene, levels = c(top20, "Other"))) %>%
    ungroup() %>% 
    group_by(sample, tool, gene) %>%
    summarise(total = sum(total), texture = texture[1])
  
  max_abun_class_en <- df_abundance_class_env %>% group_by(tool, gene) %>% summarise(q75 = quantile(total, .75)) %>% ungroup() %>% summarise(m = max(q75)) %>% pull()
  
  p <- df_abundance_class_env %>%
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
    scale_y_continuous(expand = c(0, 0), breaks = seq(0, floor(max_abun_class_en/1000) * 1000, length.out = 4),
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
    ggtitle(e) +
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
  return(p)
} 


# function to make breaks in the axis names  of the plots
lab_fn <- function(x) {
  x <- gsub("-", "-\n", x)
  x <- gsub(" ", "\n", x)
  x <- gsub("/", "/\n", x)
  x
}
