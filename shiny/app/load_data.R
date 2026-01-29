
load_results_tools <- function(DATA_DIR = "../../code_R_analysis/output_abundance_diversity_resistome", 
                               file = "results_tools.rds") {
  
  tryCatch({
    
    lst <- readRDS(file.path(DATA_DIR, file))
    
    unigenes <- as_tibble(
      do.call(rbind,
              lapply(lst, function(x)
                x[, c("query","tool","ARO","parent",
                      "parent_description","new_level","id")]))
    ) %>%
      filter(tool %in% tools_levels) %>%
      mutate(tool = factor(tool, levels = tools_levels))
    
    levels_unigenes <- unigenes %>%
      group_by(query) %>%
      slice_head(n = 1) %>%
      ungroup() %>%
      group_by(new_level) %>%
      summarise(n = n(), .groups = "drop") %>%
      arrange(desc(n)) %>%
      pull(new_level)
    
    recall_fnr <- create_class_overlaps(unigenes)
    
    list(
      lst = lst,
      unigenes = unigenes,
      recall_fnr = recall_fnr,
      levels_unigenes = levels_unigenes
    )
    
  }, error = function(e) {
    
    message("Error in load_results_tools(): ", e$message)
    
    NULL
  })
}



load_metadata <- function(DATA_DIR = "../../data", 
                          file = "metadata_GMGC10.sample.meta.tsv") {
  
  tryCatch({
    
    metadata <- read.delim("../../data/metadata_GMGC10.sample.meta.tsv")
    
    list(
      metadata = metadata
    )
    
  }, error = function(e) {
    
    message("Error in load_metadata(): ", e$message)
    
    NULL
  })
}


load_abundances <- function(DATA_DIR = "../../code_R_analysis/output_abundance_diversity_resistome", 
                            file1 = "abundance_diversity_part1.rds", 
                            file2 = "abundance_diversity_part2.rds",
                            file3 = "abundance_diversity_part3.rds",
                            metadata){
  tryCatch({
    # Load abundance files 
    abundance <- bind_rows(
      readRDS(file.path(DATA_DIR, "abundance_diversity_part1.rds")),
      readRDS(file.path(DATA_DIR, "abundance_diversity_part2.rds")),
      readRDS(file.path(DATA_DIR, "abundance_diversity_part3.rds"))
    )
    
    # the column distinct_unigenes_rarefied <- alpha diversity (number of different genes after rarefaction) 
    # no need to complete information here
    
    abundance <- abundance %>% 
      mutate(unigenes = distinct_unigenes_rarefied) %>%  # alpha diversity (number of different genes after rarefaction) 
      mutate(habitat = factor(habitat, levels = EN), # convert to factors for ordering the plots
             habitat2 = factor(SO[habitat], levels = h2), # convert to factors for ordering the plots
             tool = factor(tool, levels = tools_levels)) %>%  # convert to factors for ordering the plots
      mutate(location = ifelse(habitat2 %in% c("humans","mammals","wastewater","built-environment"), "human-related","external")) %>% 
      mutate(location = factor(location, levels = c("human-related","external"))) %>% 
      filter(tool %in% tools_levels & !habitat %in% not_env) %>%
      filter(aggregation %in% "new_level") %>%  # take only gene class aggregation
      mutate(habitat2 = factor(as.character(habitat2), h2))
    
    # normed10m is abundance, unigenes is the diversity (number of different genes )
    # we need to complete information, we need to make sure all samples appear in all tools 
    
    abundance_tool_sample <- abundance %>%
      group_by(tool, sample, habitat, habitat2) %>%  
      summarise(normed10m = sum(normed10m), unigenes = sum(unigenes)) %>% # sum the abundance and diversity
      ungroup() %>% 
      complete(sample, tool) %>% # complete with NAs
      left_join(abundance %>% select(sample, habitat, habitat2) %>% 
                  distinct(), by = "sample") %>% # get habitat and habitat2 
      mutate(habitat  = coalesce(habitat.x, habitat.y), 
             habitat2 = coalesce(habitat2.x, habitat2.y)) %>%
      select(-habitat.x, -habitat.y, -habitat2.x, -habitat2.y) %>% 
      mutate(normed10m = replace_na(normed10m, 0)) %>%  # change NAs to 0
      mutate(unigenes = replace_na(unigenes, 0)) %>% # change NAs to 0
      arrange(tool, sample)
    
    # normed10m is abundance, unigenes is the diversity (number of different genes )
    # we need to complete information, we need to make sure all samples have all gene classes and appear in all tools 
    
    abundance_class <- abundance %>% 
      ungroup() %>% 
      filter(tool %in% tools_levels, 
             !habitat %in% not_env, aggregation %in% "new_level") %>% 
      mutate(gene = factor(gene), sample = factor(sample), 
             tool = factor(tool, levels = tools_levels)) %>%
      select(!c(raw, raw_unique, scaled, distinct_unigenes_rarefied, 
                distinct_unigenes_raw, distinct_unigenes_raw_unique)) %>% 
      complete( sample, gene, tool, 
                fill = list(normed10m = 0,  unigenes = 0)) %>% ## complete for all samples, tools, and gene classes
      mutate(aggregation = "new_level", 
             gene = as.character(gene),
             habitat = metadata$habitat[match(sample, metadata$sample_id)],
             habitat2 = metadata$habitat2[match(sample, metadata$sample_id)],
             location = metadata$location[match(sample, metadata$sample_id)])  # not so dplyr way to fetch habitat and environment
    
    list(
      abundance = abundance,
      abundance_tool_sample = abundance_tool_sample,
      abundance_class = abundance_class)
    
  }, error = function(e) {
    message("Error in load_abundances(): ", e$message)
    NULL
  })
}


# load_pan_core <- function(DATA_DIR = "../../code_R_analysis/output_abundance_diversity_resistome", 
#                           core_file = "core_resistome.rds", 
#                           pan_file = "pan_resistome.rds", 
#                           data_list = list()){
#   tryCatch({
#     core <- readRDS(file.path(DATA_DIR, "core_resistome.rds"))
#     pan <- readRDS(file.path(DATA_DIR, "pan_resistome.rds"))
#     
#     core <- core %>% 
#       rename(new_level = new_level_centroid, 
#              X = centroid) %>% 
#       filter(tool %in% tools_levels, 
#              !habitat %in% not_env) %>% 
#       mutate(habitat = factor(habitat, levels = EN2), 
#              tool = factor(tool, levels =  tools_levels)) %>% 
#       mutate(tool = factor(tool, levels = tools_levels))
#     
#     pan <- pan %>% 
#       filter(tool %in% tools_levels, 
#              !habitat %in% not_env, 
#              aggregation %in% "new_level_centroid") %>% 
#       mutate(habitat = factor(habitat, levels = EN2), 
#              tool = factor(tool, levels =  tools_levels)) %>% 
#       mutate(tool = factor(tool, levels = tools_levels))
#     
#     sumpan2 <- pan %>% ungroup() %>% 
#       group_by(tool, habitat, aggregation, epoch) %>% 
#       summarise(s = sum(unigenes)) %>%
#       ungroup() %>% 
#       group_by(tool, habitat, aggregation) %>% 
#       summarise(md = median(s), mn = mean(s), sd = sd(s))
#     
#     
#     data_list$core = core
#     data_list$pan = pan
#     data_list$sumpan2 = sumpan2
#     
#   }, error = function(e) {
#     print("Error in core/pan")
#     data_list = ifelse(length(data_list) == 0, NULL, data_list)
#   })
#   return(data_list)
# }
# 



