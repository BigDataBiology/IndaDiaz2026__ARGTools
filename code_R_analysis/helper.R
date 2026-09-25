
# function to extract the legend from ggplots
g_legend <- function(a.gplot){
  tmp <- ggplotGrob(a.gplot)
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

# function to extract the core resistome based on:
# minimum proportion of metagenomic samples in a subsample needed to count an ARG as subsample core, and
# minimum number of subsample cores an ARG needs to be in to be. counted in the core-resistome

sum_core_adjust <- function(core, cnt_subset = 900, threshold_samples = 0.5){
  return(core %>% 
           filter(cut %in% threshold_samples & cnt > cnt_subset) %>% 
           ungroup() %>% 
           group_by(gene_class, tool, habitat) %>% 
           summarise(unigenes = n_distinct(X)))
}

# gene classes

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
# conversion of gene classes
gene_classes0 <- gene_classes
gene_classes0 <- gene_classes0 %>% mutate(new = ifelse(new %in% "class A beta-lactamase", "A beta-lactamase", new)) %>% 
  mutate(new = ifelse(new %in% "class B beta-lactamase", "B beta-lactamase", new)) %>% 
  mutate(new = ifelse(new %in% "class C beta-lactamase", "C beta-lactamase", new)) %>%
  mutate(new = ifelse(new %in% "class D beta-lactamase", "D beta-lactamase", new)) %>%
  mutate(new = ifelse(new %in% "beta-lactam modulation resistance", "beta-lactam modulation", new))


gene_classes_list <- gene_classes0$new
rm(gene_classes0)
names(gene_classes_list) <- gene_classes$new

# gene class vector
gene_classes <- setNames(as.list(gene_classes$new), gene_classes$old)

# function to calculate the CSC 
# qc_ref genes reported in the gene class by the reference tool 
# q_ref all genes reported by the reference tool 
# qc_comp genes reported in the gene class by the tool we are comparing against
# q_comp all genes reported by the tool we are comparing against 
new_intersect_lists <- function(qc_ref, q_ref, qc_comp, q_comp){
  A <- unlist(qc_ref)
  B <- unlist(q_ref)
  A_complement <- setdiff(B, A)
  C <- unlist(qc_comp)
  D <- unlist(q_comp)
  C_complement <- setdiff(D, C)
  x1 <- intersect(C_complement, A) # genes in a different class in comparison group but right class in reference group
  x2 <- setdiff(C, A_complement) # genes in the right class in comparison group and not in any other class in reference group
  comp_class <- unique(union(x1, x2))
  r <- ifelse(length(A) == 0 & length(comp_class) == 0, NA, 
              ifelse(length(A) == 0 & length(comp_class) != 0, NA, 
              ifelse(length(A) != 0 & length(comp_class) == 0, NA, 
                     length(intersect(A, comp_class)) / length(comp_class))))
  return(r)
}

# Calculate the CSC  for all tools and gene classes in unigenes 
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
    group_by(gene_class, tool) %>%
    summarise(query = list(query), .groups = "drop")  
  
  pairwise <- sets1 %>%
    group_by(gene_class) %>%
    summarise(pairs = list(expand_grid(tool_ref = tool, tool_comp = tool)), .groups = "drop") %>%
    unnest(pairs)
  
  JI_class_other <- pairwise %>%
    left_join(sets1, by = c("gene_class", "tool_ref" = "tool")) %>%
    rename(qc_ref = query) %>%
    left_join(sets1, by = c("gene_class", "tool_comp" = "tool")) %>%
    rename(qc_comp = query) %>%
    left_join(sets0, by = c( "tool_ref" = "tool")) %>%
    rename(q_ref = query) %>% 
    left_join(sets0, by = c( "tool_comp" = "tool")) %>%
    rename(q_comp = query)
  
  JI_class_other <- JI_class_other %>% 
    rowwise() %>%
    mutate(csc = new_intersect_lists(qc_ref, q_ref, qc_comp, q_comp)) %>% 
    rowwise() %>% 
    mutate(ref_n_class = length(qc_ref), comp_n_class = length(qc_comp), ref_n_all = length(q_ref), comp_n_all = length(q_comp)) %>% 
    ungroup() %>% 
    filter(tool_ref != tool_comp) %>% 
    select(-c(qc_ref, qc_comp, q_ref, q_comp))
  return(JI_class_other)
}

# calculates the Jaccard Index for all tools in unigenes
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
    mutate(csc = map2_dbl(values1, values2, ~ length(intersect(.x, .y)) / length( .y))) %>%
    select(tool_ref, tool_comp, jaccard, csc) 
  
  return(JI_all)
  
}


# calculates the CSC without assigning gene classes from reference tool to the tool we compare against 
create_class_overlaps_no_shuffling <- function(unigenes){
  tools_per_unigene <- unigenes %>% ungroup()  %>% 
    arrange(query) %>% 
    group_by(query) %>% 
    mutate(n_tools = n_distinct(tool)) %>% 
    mutate(single = (n_tools ==1)) 
  
  sets0 <- tools_per_unigene %>%
    group_by(tool) %>%
    summarise(query = list(query), .groups = "drop")
  
  sets1 <- tools_per_unigene %>%
    group_by(gene_class, tool) %>%
    summarise(query = list(query), .groups = "drop")  
  
  pairwise <- sets1 %>%
    group_by(gene_class) %>%
    summarise(pairs = list(expand_grid(tool_ref = tool, tool_comp = tool)), .groups = "drop") %>%
    unnest(pairs)
  
  JI_class_other <- pairwise %>%
    left_join(sets1, by = c("gene_class", "tool_ref" = "tool")) %>%
    rename(qc_ref = query) %>%
    left_join(sets1, by = c("gene_class", "tool_comp" = "tool")) %>%
    rename(qc_comp = query) %>%
    left_join(sets0, by = c( "tool_ref" = "tool")) %>%
    rename(q_ref = query) %>% 
    left_join(sets0, by = c( "tool_comp" = "tool")) %>%
    rename(q_comp = query)
  
  JI_class_other <- JI_class_other %>% 
    rowwise() %>%
    mutate(csc = length(intersect(qc_ref, qc_comp))/ length(qc_comp)) %>% 
    rowwise() %>% 
    mutate(ref_n_class = length(qc_ref), comp_n_class = length(qc_comp), ref_n_all = length(q_ref), comp_n_all = length(q_comp)) %>% 
    ungroup() %>% 
    filter(tool_ref != tool_comp) %>% 
    select(-c(qc_ref, qc_comp, q_ref, q_comp))
  return(JI_class_other)
}





plot_db_discrepancy <- function(unigenes,
                                  db_cluster,
                                  theme1,
                                  JI_db,
                                  JI_all_db,
                                  JI_all,
                                tool_a = "RGI-DIAMOND",
                                  tool_b = "ABRicate-CARD",
                                  pipeline_a_label = "RGI",
                                  pipeline_b_label = "ABRicate-CARD",
                                  pattern_a = "none",
                                  pattern_b = "stripe",
                                  both_label = "Both") {
  
  same_pair <- function(col1, col2, a, b) {
    (col1 == a & col2 == b) | (col1 == b & col2 == a)
  }
  
  # Jaccard index between tool_a/tool_b restricted to rows where rank_col == rank_value
  rank_jaccard <- function(data, rank_col, rank_value, a, b) {
    subset_data <- data %>%
      ungroup() %>%
      filter(.data[[rank_col]] == rank_value) %>%
      select(tool, query) %>%
      distinct()
    set_a <- subset_data$query[subset_data$tool == a]
    set_b <- subset_data$query[subset_data$tool == b]
    denom <- length(union(set_a, set_b))
    if (denom == 0) return(NA_real_)
    length(intersect(set_a, set_b)) / denom
  }
  
  db_cluster_query <- db_cluster %>%
    mutate(tools_db = factor(tools_db[tool], levels = tools_db_factor)) %>%
    filter(tool %in% c(.env$tool_a, .env$tool_b)) %>%
    group_by(cluster_99) %>%
    mutate(n_tool = n_distinct(tool)) %>%
    mutate(label_gene = factor(
      ifelse(n_tool == 2, .env$both_label,
             ifelse(.env$tool_a %in% tool, .env$pipeline_a_label, .env$pipeline_b_label)),
      levels = c(.env$pipeline_a_label, .env$pipeline_b_label, .env$both_label)
    )) %>%
    slice_head(n = 1)
  
  n_tools_by_cluster <- db_cluster %>%
    mutate(tools_db = factor(tools_db[tool], levels = tools_db_factor)) %>%
    filter(tool %in% c(.env$tool_a, .env$tool_b)) %>%
    group_by(cluster_99) %>% 
    mutate(n = n_distinct(tool))
  
  unigenes_ref <- unigenes %>%
    filter(tool %in% basic_tools, tool %in% c(.env$tool_a, .env$tool_b)) %>%
    mutate(label_gene = db_cluster_query$label_gene[match(cluster_99, db_cluster_query$cluster_99)]) %>%
    group_by(query, cluster_99) %>%
    mutate(present_in_others_db = n_tools_by_cluster$n[match(cluster_99, n_tools_by_cluster$cluster_99)]) %>% 
    #mutate(detected_by_both = ifelse(detected_by_both > 1, "Both pipelines", "Single pipeline")) %>% 
    mutate(detected_by_both = n_distinct(tool) > 1) %>% #   , 
                                     #"Both pipelines", "Single pipeline")) %>% 
    ungroup() %>%
    group_by(query) %>%
    mutate(detected_by_both_query = n_distinct(tool) > 1) %>%  
                                                  #"Reported under different ref gene", 
                                                  #"Single pipeline hit")) %>%
    ungroup() %>% 
    mutate(detected_by_both = ifelse(detected_by_both, "Reported by both pipelines \nsame reference gene",
                                     ifelse(detected_by_both_query, 
                                            "Reported by both pipelines \ndifferent reference gene", 
                                            ifelse(present_in_others_db > 1, "Reported by a single pipeline \nbut ref gene exists in both", "Reported by a \nsingle pipeline")))) %>% 
    mutate(detected_by_both = factor(detected_by_both, levels = c("Reported by a \nsingle pipeline",
                                                                  "Reported by a single pipeline \nbut ref gene exists in both",
                                                                  "Reported by both pipelines \ndifferent reference gene",
                                                                  "Reported by both pipelines \nsame reference gene")))
  
  pipeline_levels <- c(pipeline_a_label, pipeline_b_label)
  pattern_values  <- setNames(c(pattern_a, pattern_b), pipeline_levels)
  
  fill_values <- setNames(pal_10_complete[c(1, 2, 3)], c(pipeline_a_label, pipeline_b_label, both_label))
  
  unigenes_ref_plot <- unigenes_ref %>%
    ggplot(aes(x = detected_by_both,
               pattern = factor(ifelse(tool == .env$tool_a, .env$pipeline_a_label, .env$pipeline_b_label),
                                levels = pipeline_levels),
               fill = label_gene)) +
    geom_bar_pattern(position = position_dodge2(preserve = "single"),
                     color = "black",
                     width = 0.8, pattern_color = "black", pattern_fill = pattern_fill,
                     pattern_density = 0.001,
                     pattern_spacing = 0.01,
                     pattern_size = 0.3) +
    scale_pattern_manual(values = pattern_values) +
    xlab("Unigene found in:") +
    ylab("Number of unigenes reported as ARG") +
    scale_fill_manual(values = fill_values, guide = guide_legend(nrow = 1)) +
    labs(fill = "Hit reference gene is in", pattern = "Pipeline")
  
  id_ref <- unigenes_ref %>%
    ggplot(aes(x = detected_by_both, y = id,
               pattern = factor(ifelse(tool == .env$tool_a, .env$pipeline_a_label, .env$pipeline_b_label),
                                levels = pipeline_levels),
               fill = label_gene)) +
    geom_boxplot_pattern(position = position_dodge2(preserve = "single"), outlier.shape = NA,
                         color = "black",
                         width = 0.8, pattern_color = "black",
                         pattern_density = 0.001,
                         pattern_spacing = 0.01,
                         pattern_size = 0.3) +
    scale_pattern_manual(values = pattern_values) +
    xlab("") +
    ylab("Identity level alignment") +
    labs(fill = "Hit reference gene is in", pattern = "Pipeline") +
    scale_fill_manual(values = fill_values, guide = guide_legend(nrow = 1)) +
    ylim(c(0, 100))
  
  jaccard_data <- bind_rows(
    JI_db %>%
      filter(same_pair(tool1, tool2, .env$tool_a, .env$tool_b)) %>%
      mutate(lev = "Reference DB"),
    JI_all_db %>%
      filter(same_pair(tool_ref, tool_comp, .env$tool_a, .env$tool_b)) %>%
      rename(tool1 = tool_ref, tool2 = tool_comp) %>%
      select(-csc) %>%
      mutate(lev = "Reference gene \nwith hits in GMGC"),
    JI_all %>%
      filter(same_pair(tool_ref, tool_comp, .env$tool_a, .env$tool_b)) %>%
      rename(tool1 = tool_ref, tool2 = tool_comp) %>%
      select(-csc) %>%
      mutate(lev = "Reported unigenes")
  ) %>%
    mutate(point_type = "Observed")
  
  expected_value <- unigenes_ref %>%
    ungroup() %>%
    summarise(expected = mean(label_gene == .env$both_label)) %>%
    pull(expected)
  
  expected_row <- tibble(lev = "Expected from \nreference genes hits",
                         jaccard = expected_value,
                         point_type = "Expected")
  
  risk_levels <- c("I", "II", "III", "IV")
  
  blastp_rows <- tibble(
    lev = paste0("BLASTp risk ", risk_levels),
    jaccard = sapply(risk_levels, function(r) rank_jaccard(unigenes_ref, "rank_highest_bit_80", r, tool_a, tool_b)),
    point_type = "BLASTp risk"
  )
  
  aro_rows <- tibble(
    lev = paste0("ARO risk ", risk_levels),
    jaccard = sapply(risk_levels, function(r) rank_jaccard(unigenes_ref, "rank_aro", r, tool_a, tool_b)),
    point_type = "ARO risk"
  )
  
  jaccard_all <- bind_rows(jaccard_data, expected_row, blastp_rows, aro_rows) %>%
    mutate(lev = factor(lev, levels = c(
      "Reference DB", "Reference gene \nwith hits in GMGC", "Reported unigenes",
      "Expected from \nreference genes hits",
      paste0("BLASTp risk ", risk_levels),
      paste0("ARO risk ", risk_levels)
    )))
  
  jaccard_ref <- jaccard_all %>%
    ggplot(aes(x = lev, y = jaccard, shape = point_type)) +
    geom_point(size = 3) +
    scale_shape_manual(values = c("Observed" = 16, "Expected" = 17,
                                  "BLASTp risk" = 15, "ARO risk" = 8)) +
    ylim(c(0, 1)) +
    xlab("") +
    ylab("Jaccard index") +
    guides(fill = "none", alpha = "none") +
    theme1
  
  card_discrepancy <-
    (unigenes_ref_plot + theme1 | id_ref + theme1 | jaccard_ref) /
    (patchwork::wrap_elements(
      full = g_legend(unigenes_ref_plot + theme1 +
                        theme(legend.position = "bottom", legend.box = "vertical")))) +
    patchwork::plot_layout(heights = c(5, 1))
  
  return(card_discrepancy)
}




plot_rank_distribution <- function(unigenes,
                                   db_cluster,
                                   theme1,
                                   tool_a = "RGI-DIAMOND",
                                   tool_b = "ABRicate-CARD",
                                   pipeline_a_label = "RGI",
                                   pipeline_b_label = "ABRicate-CARD",
                                   pattern_a = "none",
                                   pattern_b = "stripe",
                                   both_label = "Both") {
  
  db_cluster_card <- db_cluster %>%
    mutate(tools_db = factor(tools_db[tool], levels = tools_db_factor)) %>%
    filter(tool %in% c(.env$tool_a, .env$tool_b)) %>%
    group_by(cluster_99) %>%
    mutate(n_tool = n_distinct(tool)) %>%
    mutate(label_gene = factor(
      ifelse(n_tool == 2, .env$both_label,
             ifelse(.env$tool_a %in% tool, .env$pipeline_a_label, .env$pipeline_b_label)),
      levels = c(.env$pipeline_a_label, .env$pipeline_b_label, .env$both_label)
    )) %>%
    slice_head(n = 1)
  
  unigenes_card <- unigenes %>%
    filter(tool %in% basic_tools, tool %in% c(.env$tool_a, .env$tool_b)) %>%
    mutate(label_gene = db_cluster_card$label_gene[match(cluster_99, db_cluster_card$cluster_99)]) %>%
    group_by(query, cluster_99) %>%
    mutate(detected_by_both = factor(ifelse(n_distinct(tool) > 1, "Both pipelines", "Single pipeline"),
                                     levels = c("Single pipeline", "Both pipelines"))) %>%
    ungroup()
  
  pipeline_levels <- c(pipeline_a_label, pipeline_b_label)
  pattern_values  <- setNames(c(pattern_a, pattern_b), pipeline_levels)
  fill_values     <- setNames(pal_10_complete[c(1, 2, 3)], c(pipeline_a_label, pipeline_b_label, both_label))
  
  rank_levels     <- c("I", "II", "III", "IV")
  detected_levels <- c("Single pipeline", "Both pipelines")
  x_levels <- unlist(lapply(paste0("Rank ", rank_levels), function(r) paste0(r, "\n", detected_levels)))
  
  long_data <- bind_rows(
    unigenes_card %>% mutate(rank_type = "rank_aro", rank_value = as.character(rank_aro)),
    unigenes_card %>% mutate(rank_type = "rank_highest_bit_80", rank_value = as.character(rank_highest_bit_80))
  ) %>%
    filter(rank_value %in% rank_levels) %>%
    mutate(
      x_lab = factor(paste0("Rank ", rank_value, "\n", detected_by_both), levels = x_levels),
      rank_type = factor(rank_type, levels = c("rank_aro", "rank_highest_bit_80"))
    )
  
  long_data %>%
    ggplot(aes(x = x_lab,
               pattern = factor(ifelse(tool == .env$tool_a, .env$pipeline_a_label, .env$pipeline_b_label),
                                levels = pipeline_levels),
               fill = label_gene)) +
    geom_bar_pattern(position = position_dodge2(preserve = "single"),
                     color = "black", pattern_color = "black", pattern_fill = pattern_fill,
                     pattern_density = 0.001, pattern_spacing = 0.01, pattern_size = 0.3) +
    scale_pattern_manual(values = pattern_values) +
    scale_fill_manual(values = fill_values, guide = guide_legend(nrow = 1)) +
    facet_grid(rank_type ~ ., scales = "free_y") +
    xlab("") +
    ylab("Number of unigenes") +
    labs(fill = "Hit reference gene is in", pattern = "Pipeline") +
    theme1
}




compute_expected_jaccard <- function(unigenes, db_cluster) {
  
  tools <- unique(as.character(unigenes$tool[unigenes$tool %in% basic_tools]))
  pairs <- combn(tools, 2, simplify = FALSE)
  
  purrr::map_dfr(pairs, function(p) {
    a <- p[1]; b <- p[2]
    
    # which cluster_99 (reference-DB level) are hit by BOTH tools of this pair
    both_clusters <- db_cluster %>%
      filter(tool %in% c(a, b)) %>%
      group_by(cluster_99) %>%
      summarise(n_tool = n_distinct(tool), .groups = "drop") %>%
      filter(n_tool == 2) %>%
      pull(cluster_99)
    
    # all unigenes rows for this pair (both tools' rows combined)
    unigenes_pair <- unigenes %>%
      filter(tool %in% basic_tools, tool %in% c(a, b))
    
    expected <- mean(unigenes_pair$cluster_99 %in% both_clusters)
    
    tibble(tool1 = a, tool2 = b, expected_jaccard = expected)
  })
}



get_unigene_classification <- function(unigenes,
                                       db_cluster,
                                       tool_a = "RGI-DIAMOND",
                                       tool_b = "ABRicate-CARD",
                                       pipeline_a_label = "RGI",
                                       pipeline_b_label = "ABRicate-CARD",
                                       both_label = "Both") {
  
  db_cluster_query <- db_cluster %>%
    filter(tool %in% c(.env$tool_a, .env$tool_b)) %>%
    group_by(cluster_99) %>%
    mutate(n_tool = n_distinct(tool)) %>%
    mutate(label_gene = factor(
      ifelse(n_tool == 2, .env$both_label,
             ifelse(.env$tool_a %in% tool, .env$pipeline_a_label, .env$pipeline_b_label)),
      levels = c(.env$pipeline_a_label, .env$pipeline_b_label, .env$both_label)
    )) %>%
    slice_head(n = 1) %>%
    ungroup()
  
  n_tools_by_cluster <- db_cluster %>%
    filter(tool %in% c(.env$tool_a, .env$tool_b)) %>%
    group_by(cluster_99) %>%
    summarise(n = n_distinct(tool), .groups = "drop")
  
  unigenes_ref <- unigenes %>%
    filter(tool %in% basic_tools, tool %in% c(.env$tool_a, .env$tool_b)) %>%
    mutate(label_gene = db_cluster_query$label_gene[match(cluster_99, db_cluster_query$cluster_99)]) %>%
    mutate(present_in_others_db = n_tools_by_cluster$n[match(cluster_99, n_tools_by_cluster$cluster_99)]) %>%
    group_by(query, cluster_99) %>%
    mutate(detected_by_both = n_distinct(tool) > 1) %>%
    ungroup() %>%
    group_by(query) %>%
    mutate(detected_by_both_query = n_distinct(tool) > 1) %>%
    ungroup() %>%
    mutate(detected_by_both = ifelse(
      detected_by_both, "Reported by both pipelines \nsame reference gene",
      ifelse(detected_by_both_query, "Reported by both pipelines \ndifferent reference gene",
             ifelse(present_in_others_db > 1, "Reported by a single pipeline \nbut ref gene exists in both",
                    "Reported by a \nsingle pipeline")))) %>%
    mutate(detected_by_both = factor(detected_by_both, levels = c(
      "Reported by a \nsingle pipeline",
      "Reported by a single pipeline \nbut ref gene exists in both",
      "Reported by both pipelines \ndifferent reference gene",
      "Reported by both pipelines \nsame reference gene"
    )))
  
  query_classification <- unigenes_ref %>%
    distinct(query, detected_by_both, label_gene) %>%
    mutate(tool_a = .env$tool_a, tool_b = .env$tool_b)
  
  n_dup <- query_classification %>% count(query) %>% filter(n > 1) %>% nrow()
  if (n_dup > 0) {
    warning(sprintf("%d queries have more than one (x, color) classification — likely multiple cluster_99 assignments for the same query. Inspect with count(query) %%>%% filter(n > 1).", n_dup))
  }
  
  query_classification
}




