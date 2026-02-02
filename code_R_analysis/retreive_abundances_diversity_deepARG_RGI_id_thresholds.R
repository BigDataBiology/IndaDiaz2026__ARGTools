library(dplyr)
library(tidyr)

################################################################################################################################################
# working directory

setwd("~/Documents/GitHub/arg_compare/")
lst <- readRDS("code_R_analysis/output_abundance_diversity_resistome/results_tools.rds")
lst2 <- list(deeparg.norm.prot = lst$deeparg.norm.prot, rgi.diamond.prot = lst$rgi.diamond.prot)
lst <- lst2
rm(lst2)
args_abundances <- read.delim("data/abundances/args_abundances.tsv")
metadata <- read.delim("data/metadata_GMGC10.sample.meta.tsv")
metadata <- metadata %>% mutate(sample = sample_id)
args_abundances <- args_abundances %>% left_join(metadata[,c("sample","insertsHQ", "insertsRaw")], by = "sample")
metadata <- metadata %>% filter(sample %in% args_abundances$sample)

unigenes60 <- unique(c(lst$deeparg.norm.prot$query[lst$deeparg.norm.prot$id >= 60],
                       lst$rgi.diamond.prot$query[lst$rgi.diamond.prot$id >= 60]))

unigenes70 <- unique(c(lst$deeparg.norm.prot$query[lst$deeparg.norm.prot$id >= 70],
                       lst$rgi.diamond.prot$query[lst$rgi.diamond.prot$id >= 70]))

unigenes80 <- unique(c(lst$deeparg.norm.prot$query[lst$deeparg.norm.prot$id >= 80],
                       lst$rgi.diamond.prot$query[lst$rgi.diamond.prot$id >= 80]))

args_abundances60 <- args_abundances %>% filter(X %in% unigenes60)
args_abundances70 <- args_abundances %>% filter(X %in% unigenes70)
args_abundances80 <- args_abundances %>% filter(X %in% unigenes80)

rarefaction <- function(X, raw, inserts, depth = 5e6, seed = 2025){
  set.seed(seed)
  reads_count <- ceiling(raw)
  arg_reads <- sum(reads_count)
  not_arg_reads <- inserts[1] - arg_reads
  expanded_reads <- rep(c(X,"not an ARG"), times = c(reads_count, not_arg_reads))
  
  if(depth > inserts[1]){
    depth <- inserts[1]
  }
  
  subsample <- sample(expanded_reads, depth, replace = FALSE)
  tab <- table(subsample)
  return(tab)
}

# do the rarefaction across samples
rarefied_counts60 <- args_abundances60 %>% filter(raw > 0) %>% ungroup() %>% group_by(sample) %>%
  summarise(
    rarefied = list(
      rarefaction(X, raw, insertsHQ, 5e6)
    ),
    .groups = "drop"
  ) %>%
  unnest_longer(rarefied, indices_include = TRUE) %>%
  rename(X = rarefied_id, count = rarefied) %>% 
  mutate(count = as.integer(count)) %>% 
  rename(rarified_count = count)

rarefied_counts70 <- args_abundances70 %>% filter(raw > 0) %>% ungroup() %>% group_by(sample) %>%
  summarise(
    rarefied = list(
      rarefaction(X, raw, insertsHQ, 5e6)
    ),
    .groups = "drop"
  ) %>%
  unnest_longer(rarefied, indices_include = TRUE) %>%
  rename(X = rarefied_id, count = rarefied) %>% 
  mutate(count = as.integer(count)) %>% 
  rename(rarified_count = count)

rarefied_counts80 <- args_abundances80 %>% filter(raw > 0) %>% ungroup() %>% group_by(sample) %>%
  summarise(
    rarefied = list(
      rarefaction(X, raw, insertsHQ, 5e6)
    ),
    .groups = "drop"
  ) %>%
  unnest_longer(rarefied, indices_include = TRUE) %>%
  rename(X = rarefied_id, count = rarefied) %>% 
  mutate(count = as.integer(count)) %>% 
  rename(rarified_count = count)


##


# add them to the abundance object   
args_abundances60 <- args_abundances60 %>% left_join(rarefied_counts60, by = c("sample", "X"))
args_abundances60 <- args_abundances60 %>% mutate(rarified_count = ifelse(is.na(rarified_count), 0, rarified_count))

args_abundances70 <- args_abundances70 %>% left_join(rarefied_counts70, by = c("sample", "X"))
args_abundances70 <- args_abundances70 %>% mutate(rarified_count = ifelse(is.na(rarified_count), 0, rarified_count))

args_abundances80 <- args_abundances80 %>% left_join(rarefied_counts80, by = c("sample", "X"))
args_abundances80 <- args_abundances80 %>% mutate(rarified_count = ifelse(is.na(rarified_count), 0, rarified_count))


abundance_parent <- function(abund_df, d){
  d <- d %>% filter(!is.na(parent)) 
  Y <- abund_df %>% mutate(aro = d$ARO[match(X, d$query)], 
                                  parent_description = d$parent_description[match(X, d$query)],
                                  new_level = d$new_level[match(X, d$query)])
  # we sum by ARO and gene class (not by unigene), the scale counts, 
  # raw counts, rounded raw counts, raw unique counts, normed10m, 
  # the distinct_unigenes_rarefied (diversity), 
  # distinct_unigenes_raw (diversity without rarefaction), 
  # distinct_unigenes_raw_unique (anothe type of diversity))
  
  Y_aro <- Y %>% filter(!is.na(new_level))  %>% group_by(sample, aro) %>% 
    summarise(scaled = sum(scaled), raw = sum(raw), raw_rounded = sum(ceiling(raw)), 
              raw_unique = sum(raw_unique), normed10m = sum(normed10m), 
              distinct_unigenes_rarefied = n_distinct(X[rarified_count > 0]), 
              distinct_unigenes_raw = n_distinct(X), 
              distinct_unigenes_raw_unique = n_distinct(X[raw_unique > 0])) %>% 
    mutate(tool = d$tool[1]) %>% ungroup() %>% 
    select(sample, aro, tool, scaled, raw, raw_unique, 
           normed10m, distinct_unigenes_rarefied, 
           distinct_unigenes_raw, distinct_unigenes_raw_unique) %>%
    rename(gene = aro) %>% mutate(aggregation = "ARO") %>% 
    select(sample, gene, aggregation, tool, scaled, raw, 
           raw_unique, normed10m, distinct_unigenes_rarefied, 
           distinct_unigenes_raw, distinct_unigenes_raw_unique)
  
  Y_new_level <- Y %>% filter(!is.na(new_level)) %>% group_by(sample, new_level) %>% 
    summarise(scaled = sum(scaled), raw = sum(raw), raw_rounded = sum(ceiling(raw)), 
              raw_unique = sum(raw_unique), normed10m = sum(normed10m), 
              distinct_unigenes_rarefied = n_distinct(X[rarified_count > 0]), 
              distinct_unigenes_raw = n_distinct(X), 
              distinct_unigenes_raw_unique = n_distinct(X[raw_unique > 0])) %>% 
    mutate(tool = d$tool[1]) %>% ungroup() %>% 
    select(sample, new_level, tool, scaled, raw, raw_unique, 
           normed10m, distinct_unigenes_rarefied, 
           distinct_unigenes_raw, distinct_unigenes_raw_unique) %>%
    rename(gene = new_level) %>% mutate(aggregation = "new_level") %>% 
    select(sample, gene, aggregation, tool, scaled, raw, 
           raw_unique, normed10m, distinct_unigenes_rarefied, 
           distinct_unigenes_raw, distinct_unigenes_raw_unique)
  
  Y_aro <- Y_aro %>% bind_rows(Y_new_level)
  
  return(Y_aro)
}


lst60 <- list(deeparg.norm.prot = lst$deeparg.norm.prot[lst$deeparg.norm.prot$id >= 60,], 
              rgi.diamond.prot = lst$rgi.diamond.prot[lst$rgi.diamond.prot$id >= 60,])
lst70 <- list(deeparg.norm.prot = lst$deeparg.norm.prot[lst$deeparg.norm.prot$id >= 70,], 
              rgi.diamond.prot = lst$rgi.diamond.prot[lst$rgi.diamond.prot$id >= 70,])
lst80 <- list(deeparg.norm.prot = lst$deeparg.norm.prot[lst$deeparg.norm.prot$id >= 80,], 
              rgi.diamond.prot = lst$rgi.diamond.prot[lst$rgi.diamond.prot$id >= 80,])


# calculate the abundance and diversity for all genes and habitats by tool.
lst_abundance_diversity60 <- do.call(rbind, lapply(lst60, function(d) {abundance_parent(args_abundances60, d) }))
lst_abundance_diversity70 <- do.call(rbind, lapply(lst70, function(d) {abundance_parent(args_abundances70, d) }))
lst_abundance_diversity80 <- do.call(rbind, lapply(lst80, function(d) {abundance_parent(args_abundances80, d) }))

# add habitat info
lst_abundance_diversity60 <- lst_abundance_diversity60 %>% 
  mutate(habitat = metadata$habitat[match(sample, metadata$sample_id)])

lst_abundance_diversity70 <- lst_abundance_diversity70 %>% 
  mutate(habitat = metadata$habitat[match(sample, metadata$sample_id)])

lst_abundance_diversity80 <- lst_abundance_diversity80 %>% 
  mutate(habitat = metadata$habitat[match(sample, metadata$sample_id)])

# Unigenes detected as ARG by tool
unigenes <- do.call(rbind, lapply(lst, function(x) x[,c("query", "tool", "ARO", "parent", "parent_description", "new_level", "id")])) 

# habitats
EN <- c("human gut", "human oral",  "human skin", "human nose", "human vagina", 
        "dog gut", "cat gut", "mouse gut", "pig gut", "wastewater", "marine", "freshwater",  
        "soil" , "amplicon", "isolate",  "built-environment" )

# higher habitat classifiction
SO <- c(rep("humans", 5), rep("mammals", 4),  "wastewater", "marine", "freshwater", "soil", rep("other", 3))
names(SO) <- EN

# higher habitat classifiction for abundance
lst_abundance_diversity60 <- lst_abundance_diversity60 %>% 
  mutate(habitat2 = SO[lst_abundance_diversity60$habitat])

lst_abundance_diversity70 <- lst_abundance_diversity70 %>% 
  mutate(habitat2 = SO[lst_abundance_diversity70$habitat])

lst_abundance_diversity80 <- lst_abundance_diversity80 %>% 
  mutate(habitat2 = SO[lst_abundance_diversity80$habitat])

# save abundance and diversity
saveRDS(lst_abundance_diversity60, file = "code_R_analysis/output_abundance_diversity_resistome/abundance_diversity_60.rds", compress = T)
saveRDS(lst_abundance_diversity70, file = "code_R_analysis/output_abundance_diversity_resistome/abundance_diversity_70.rds", compress = T)
saveRDS(lst_abundance_diversity80, file = "code_R_analysis/output_abundance_diversity_resistome/abundance_diversity_80.rds", compress = T)

# save the results per tool
#for(j in 1:length(lst)){
#  write.csv(lst[[j]], file = paste0("code_R_analysis/output_abundance_diversity_resistome/processed_tool.",names(lst)[j],".csv"), row.names = F)
#}


#### CORE AND PAN 
## load the unigenes clusterd at 90% with vsearch
clusters <- read.delim("cluster_vsearch/clusters.uc", header = F)
clusters <- clusters %>% filter(V1 != "C")
clusters <- clusters %>% mutate(centroid = ifelse(V10 == "*", V9, V10))

# backup abundances
#args_abundances0 <- args_abundances
#args_abundances <- args_abundances0

# add habitat to abundances
args_abundances <- args_abundances %>% 
  mutate(habitat = metadata$habitat[match(sample, metadata$sample_id)])

args_abundances60 <- args_abundances60 %>% 
  mutate(habitat = metadata$habitat[match(sample, metadata$sample_id)])

args_abundances70 <- args_abundances70 %>% 
  mutate(habitat = metadata$habitat[match(sample, metadata$sample_id)])

args_abundances80 <- args_abundances80 %>% 
  mutate(habitat = metadata$habitat[match(sample, metadata$sample_id)])

# add centroid to abundances


args_abundances <- args_abundances %>% mutate(centroid = clusters$centroid[match(X, clusters$V9)])
args_abundances <- args_abundances %>% rename(query = X)

args_abundances60 <- args_abundances60 %>% mutate(centroid = clusters$centroid[match(X, clusters$V9)])
args_abundances60 <- args_abundances60 %>% rename(query = X)

args_abundances70 <- args_abundances70 %>% mutate(centroid = clusters$centroid[match(X, clusters$V9)])
args_abundances70 <- args_abundances70 %>% rename(query = X)

args_abundances80 <- args_abundances80 %>% mutate(centroid = clusters$centroid[match(X, clusters$V9)])
args_abundances80 <- args_abundances80 %>% rename(query = X)

# add centroid to the result of each tool
lst60 <- lapply(lst60, function(x) x %>% mutate(centroid = clusters$centroid[match(query, clusters$V9)]))
lst70 <- lapply(lst70, function(x) x %>% mutate(centroid = clusters$centroid[match(query, clusters$V9)]))
lst80 <- lapply(lst80, function(x) x %>% mutate(centroid = clusters$centroid[match(query, clusters$V9)]))

# which centroids have more than 1 gene class in the cluster 

lapply(lst, function(x) x %>% 
         group_by(centroid) %>% mutate(n = n_distinct(new_level))  %>% 
         filter(n>1) %>% arrange(centroid, desc(n)) %>% 
         select(query, centroid, new_level))

# assign majority rule for the class of the centroid per tool
lst60 <- lapply(lst60, function(x) x %>% 
         group_by(centroid) %>% mutate(n = n_distinct(new_level))  %>% 
         mutate(new_level_centroid = new_level[match(centroid, query)]) %>%
         mutate(new_level_majority = names(which.max(table(new_level)))) %>%
         mutate(new_level_centroid = ifelse(is.na(new_level_centroid), new_level_majority, new_level_centroid)))


lst70 <- lapply(lst70, function(x) x %>% 
                group_by(centroid) %>% mutate(n = n_distinct(new_level))  %>% 
                mutate(new_level_centroid = new_level[match(centroid, query)]) %>%
                mutate(new_level_majority = names(which.max(table(new_level)))) %>%
                mutate(new_level_centroid = ifelse(is.na(new_level_centroid), new_level_majority, new_level_centroid)))

lst80 <- lapply(lst80, function(x) x %>% 
                group_by(centroid) %>% mutate(n = n_distinct(new_level))  %>% 
                mutate(new_level_centroid = new_level[match(centroid, query)]) %>%
                mutate(new_level_majority = names(which.max(table(new_level)))) %>%
                mutate(new_level_centroid = ifelse(is.na(new_level_centroid), new_level_majority, new_level_centroid)))



# save the result of all tools
lst60$deeparg.norm.prot$tool <- "DeepARG"
lst70$deeparg.norm.prot$tool <- "DeepARG"
lst80$deeparg.norm.prot$tool <- "DeepARG"

lst60$rgi.diamond.prot$tool <- "RGI-DIAMOND"
lst70$rgi.diamond.prot$tool <- "RGI-DIAMOND"
lst80$rgi.diamond.prot$tool <- "RGI-DIAMOND"

saveRDS(lst60,  file = "code_R_analysis/output_abundance_diversity_resistome/results_tools60.rds", compress = T)
saveRDS(lst70,  file = "code_R_analysis/output_abundance_diversity_resistome/results_tools70.rds", compress = T)
saveRDS(lst80,  file = "code_R_analysis/output_abundance_diversity_resistome/results_tools80.rds", compress = T)



# functions for core resistome

cut_size_core <- function(df, cut) {
  return(df %>% ungroup() %>% filter(p >= cut) %>% # filter by proportion of sammples with a specific gene
           mutate(cut = cut, cnt = 1) %>% 
           ungroup() %>% 
           select(centroid, new_level_centroid, tool, habitat, cut, cnt))
}

filter_samples_core <- function(args_abundances_core, d){
  DF <- args_abundances_core %>% 
    filter(query %in% d$query) %>% # filter unigenes with abundances that are in the results of the tool d
    mutate(tool = d$tool[1]) %>% # fetch tool name
    mutate(new_level_centroid = d$new_level_centroid[match(query, d$query)]) # fetch class of the centroid
  
  if (sum(is.na(DF$new_level_centroid)) > 0) {
    print(DF[is.na(DF$new_level_centroid),])
    stop("NA values detected in new_level for tool: ", d$tool[1])
  }
  
  DF <- DF %>%
    select(centroid, new_level_centroid, sample, tool, habitat) %>% 
    ungroup() %>% group_by(tool, habitat) %>% mutate(N = n_distinct(sample)) %>% # total number of samples per habitat
    ungroup() %>% group_by(centroid, habitat) %>% mutate(n = n_distinct(sample)) %>% # number of samples where a centroid appears per habitat
    ungroup() %>% mutate( p = n / N) %>% filter(p >= 0.1) %>% # proportion of samples with the centroid
    group_by(habitat, centroid) %>% 
    slice_head(n = 1) %>% ungroup() %>% select(-sample)
    return(DF)
}


core_resistome <- function(args_abundances, samples_to_collect, sed, lst, mx_sample_size, j, cuts, df) {
  set.seed(seed = sed) # random seed
  samples_to_collect2 <- samples_to_collect  %>%
    slice_sample(n = mx_sample_size,  replace = FALSE) %>%
    ungroup() %>% select(sample) %>% pull() # subsample the metagenomic samples of a habitat

  args_abundances_core <- args_abundances %>% filter(sample %in% samples_to_collect2) # subsample the metagenomic samples of a habitat
  args_abundances_core <- do.call(rbind, lapply(lst, function(d) filter_samples_core(args_abundances_core, d))) # find the proportion of samples with the centroid 
  args_abundances_core_cut <- do.call(rbind, lapply(seq_along(cuts), function(k) {cut_size_core(args_abundances_core, cuts[k])})) # filter by proportion

  df <- df %>% bind_rows(args_abundances_core_cut)
  df <- df %>% ungroup() %>% group_by(centroid, new_level_centroid, tool, habitat, cut) %>%
    summarise(cnt = sum(cnt))
  return(df)
}


seeds <- seq(2001, 2500, 1)
#seeds <- c(2001)
cuts <- c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
#cuts <- c(0.5)

samples_to_collect <- args_abundances %>%
  group_by(habitat) %>%
  distinct(sample) %>%  # Get distinct values of samples per habitat
  mutate(n_unique = n()) %>%  # Calculate the number of unique values in each group
  ungroup() %>%
  group_by(habitat)

df60 <- data.frame(centroid = NULL, new_level_centroid = NULL, tool = NULL, habitat = NULL, cut = NULL, cnt = NULL)
df70 <- data.frame(centroid = NULL, new_level_centroid = NULL, tool = NULL, habitat = NULL, cut = NULL, cnt = NULL)
df80 <- data.frame(centroid = NULL, new_level_centroid = NULL, tool = NULL, habitat = NULL, cut = NULL, cnt = NULL)

for(j in 1:length(seeds)){
  print(j)
  df60 <- core_resistome(args_abundances60 %>% filter(rarified_count > 0),  samples_to_collect, seeds[j], lst, 100, j, cuts, df60)
  df70 <- core_resistome(args_abundances70 %>% filter(rarified_count > 0),  samples_to_collect, seeds[j], lst, 100, j, cuts, df70)
  df80 <- core_resistome(args_abundances80 %>% filter(rarified_count > 0),  samples_to_collect, seeds[j], lst, 100, j, cuts, df80)
}

df60 <- df60 %>% mutate(tool = ifelse(tool == "DeepARG-aa", "DeepARG", "RGI-DIAMOND"))
df70 <- df70 %>% mutate(tool = ifelse(tool == "DeepARG-aa", "DeepARG", "RGI-DIAMOND"))
df80 <- df80 %>% mutate(tool = ifelse(tool == "DeepARG-aa", "DeepARG", "RGI-DIAMOND"))

saveRDS(df60, file = "code_R_analysis/output_abundance_diversity_resistome/core_resistome60.rds", compress = T)
saveRDS(df70, file = "code_R_analysis/output_abundance_diversity_resistome/core_resistome70.rds", compress = T)
saveRDS(df80, file = "code_R_analysis/output_abundance_diversity_resistome/core_resistome80.rds", compress = T)



filter_samples_pan <- function(args_abundances_pan, d, j){

  DF <- args_abundances_pan %>% 
    filter(query %in% d$query) %>% 
    mutate(tool = d$tool[1]) %>%
    mutate(new_level_centroid = d$new_level_centroid[match(query, d$query)])
  
  Y_new_level <- DF %>% group_by(tool, habitat, new_level_centroid) %>% 
    summarise(unigenes = n_distinct(centroid)) %>%
    mutate(aggregation = "new_level_centroid", epoch = j) %>% ungroup() %>%
    rename(gene_class = new_level_centroid)
  
  return(Y_new_level)
}

pan_resistome <- function(df, samples_to_collect, sed, lst, mx_sample_size, j) {
  set.seed(seed = sed)
  
  samples_to_collect2 <- samples_to_collect  %>%
    slice_sample(n = mx_sample_size,  replace = FALSE) %>%  
    ungroup() %>% pull(sample)
  
  args_abundances_pan <- df %>% filter(sample %in% samples_to_collect2)
  args_abundances_pan_filter <- bind_rows(lapply(lst, function(d) filter_samples_pan(args_abundances_pan, d, j)))
  
  return(args_abundances_pan_filter)
}

args_abundances_rarified60 <- args_abundances60 %>% filter(rarified_count > 0)
args_abundances_rarified70 <- args_abundances70 %>% filter(rarified_count > 0)
args_abundances_rarified80 <- args_abundances80 %>% filter(rarified_count > 0)

df.pan.list60 <- vector("list", length(seeds))
df.pan.list70 <- vector("list", length(seeds))
df.pan.list80 <- vector("list", length(seeds))

for(j in seq_along(seeds)){
  print(j)
  df.pan.list60[[j]] <- pan_resistome(args_abundances_rarified60, samples_to_collect, seeds[j], lst60, 100, j)
  df.pan.list70[[j]] <- pan_resistome(args_abundances_rarified70, samples_to_collect, seeds[j], lst70, 100, j)
  df.pan.list80[[j]] <- pan_resistome(args_abundances_rarified80, samples_to_collect, seeds[j], lst80, 100, j)
}

df.pan2.60 <- bind_rows(df.pan.list60)
df.pan2.70 <- bind_rows(df.pan.list70)
df.pan2.80 <- bind_rows(df.pan.list80)

df.pan2.60 <- df.pan2.60 %>% mutate(tool = ifelse(tool == "DeepARG-aa", "DeepARG", "RGI-DIAMOND"))
df.pan2.70 <- df.pan2.70 %>% mutate(tool = ifelse(tool == "DeepARG-aa", "DeepARG", "RGI-DIAMOND"))
df.pan2.80 <- df.pan2.80 %>% mutate(tool = ifelse(tool == "DeepARG-aa", "DeepARG", "RGI-DIAMOND"))


saveRDS(df.pan2.60, file = "code_R_analysis/output_abundance_diversity_resistome/pan_resistome60.rds", compress = T)
saveRDS(df.pan2.70, file = "code_R_analysis/output_abundance_diversity_resistome/pan_resistome70.rds", compress = T)
saveRDS(df.pan2.80, file = "code_R_analysis/output_abundance_diversity_resistome/pan_resistome80.rds", compress = T)


