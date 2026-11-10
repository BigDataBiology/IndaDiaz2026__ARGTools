library(dplyr)
library(tidyverse)


#setwd("~/Documents/GitHub/arg_compare/") 
options(dplyr.summarise.inform = FALSE)
# source("helper.R")


# results per individual tool 
lst <- readRDS("code_R_analysis/output_abundance_diversity_resistome/results_tools.rds")

# metagenomes' metadata
metadata <- read.delim("data/metadata_GMGC10.sample.meta.tsv")
metadata <- metadata %>% select(sample_id, habitat, insertsHQ)

# abundance <- readRDS("output_abundance_diversity_resistome/abundance_diversity.rds")
abundance <- readRDS("code_R_analysis/output_abundance_diversity_resistome/abundance_diversity.rds")
abundance <- abundance %>% rename(geneclass = gene, pipeline = tool) %>% 
  select(-c(aggregation,new_level))

conversion_file <- readRDS("code_R_analysis/output_abundance_diversity_resistome/conversion_ARO_parent_new_level.rds")
conversion_file <- conversion_file %>% rename(geneclass = new_level)


lst2 <- lapply(lst, function(x) x <- x %>% rename(geneclass_argcompare = new_level,
                                                  pipeline = tool))
lst2 <- lst2[1:15]

lapply(names(lst2), function(nm) {
  write.csv(lst2[[nm]], paste0("code_R_analysis/data_to_Zenodo/pipelines_output/pipeline_",nm,".csv"), row.names = FALSE)
})


write.csv(abundance, file = gzfile("code_R_analysis/data_to_Zenodo/abundance_richness.csv.gz"), row.names = F)
write.csv(metadata, file = "code_R_analysis/data_to_Zenodo/metagenomes_metadata.csv", row.names = F)
write.csv(conversion_file, file = "code_R_analysis/data_to_Zenodo/conversion_aro_geneclass.csv", row.names = F)


reported_unigenes_as_ARG_per_habitat <- read.delim("code_R_analysis/output_abundance_diversity_resistome/reported_unigenes_as_ARG_per_habitat.csv", sep = ",")
reported_unigenes_as_ARG_per_habitat <- reported_unigenes_as_ARG_per_habitat %>% filter(habitat %in% c("human gut", "human oral",  "human skin", 
                                                               "human nose", "human vagina", 
                                                               "dog gut", "cat gut", "mouse gut", 
                                                               "pig gut", "wastewater", "marine", 
                                                               "freshwater", "soil" ))
reported_unigenes_as_ARG_per_habitat <- reported_unigenes_as_ARG_per_habitat %>% rename(unigene = X)

write.csv(reported_unigenes_as_ARG_per_habitat, file = "code_R_analysis/data_to_Zenodo/reported_unigenes_as_ARG_per_habitat.csv", row.names = F)


basic_tools <- c(
  "RGI-DIAMOND", "ABRicate-CARD","AMRFinderPlus", "ABRicate-NCBI",
  "ResFinder", "ABRicate-ResFinder", 
  "ABRicate-ARGANNOT", "ABRicate-MEGARes", "fARGene", "DeepARG" )
tool_map <- c(
  "abricate-argannot" = "ABRicate-ARGANNOT" ,
  "abricate-card"      = "ABRicate-CARD",
  "abricate-megares"   = "ABRicate-MEGARes",
  "abricate-ncbi"      =  "ABRicate-NCBI",
  "abricate-resfinder" = "ABRicate-ResFinder",
  "amrfinderplus"      = "AMRFinderPlus",
  "deeparg"            = "DeepARG",
  "resfinder"          = "ResFinder",
  "rgi-card"        = "RGI-DIAMOND",
  "fargene" = "fARGene"
)
tag_to_tool <- c(
  "resfinder"          = "ResFinder",
  "abricate-resfinder" = "ABRicate-ResFinder",
  "abricate-argannot"  = "ABRicate-ARGANNOT",
  "abricate-megares"   = "ABRicate-MEGARes",
  "abricate-card"      = "ABRicate-CARD",
  "abricate-ncbi"      = "ABRicate-NCBI",
  "deeparg"            = "DeepARG",
  "rgi-card"           = "RGI-DIAMOND",
  "amrfinderplus"      = "AMRFinderPlus",
  "fargene" = "fARGene"
)

# db_cluster <- read.delim("db_cluster/nested_with_coverage_out/nested_cluster_membership.tsv") 
# db_cluster <- db_cluster %>% mutate(tool = sapply(strsplit(protein_id, split = "@@@"), function(x) x[1]),
#                                     gene = sapply(strsplit(protein_id, split = "@@@"), function(x) x[2])) %>% 
#   mutate(tool = factor(tool_map[tool], levels = basic_tools)) 
# 
# 
# db_cluster <- db_cluster %>% mutate(tool = sapply(strsplit(protein_id, split = "@@@"), function(x) x[1]),
#                                     gene = sapply(strsplit(protein_id, split = "@@@"), function(x) x[2]))
# 
# norm_prime <- function(x) {
#   x <- gsub("###", "'", x, fixed = TRUE)
#   x <- gsub("\u2019", "'", x, fixed = TRUE)
#   x <- gsub("\u2032", "'", x, fixed = TRUE)
#   trimws(x)
# }
# 
# db_cluster <- db_cluster %>%
#   mutate(gene_revised = ifelse(tool %in% "rgi-card", 
#                                as.numeric(sub(".*\\|ID:([0-9]+)\\|.*", "\\1", gene)),
#                                NA))
# 
# db_cluster <- db_cluster %>%
#   mutate(gene_revised = ifelse(tool %in% c("abricate-argannot", "abricate-card","abricate-ncbi", "abricate-megares","abricate-resfinder"), 
#                                sapply(strsplit(gene, split = "~~~", fixed = TRUE), `[`, 2),gene_revised)) %>%
#   mutate(gene_revised = case_when(
#     tool == "abricate-argannot"  ~ sub("\\([^)]*\\)", "", gene_revised),
#     tool == "abricate-resfinder" ~ sub("_[0-9]+$", "", gene_revised),
#     TRUE ~ gene_revised
#   ))
# 
# extract_resfinder_accession <- function(x) {
#   refseq_style <- str_extract(x, "[A-Z]{1,4}_[0-9A-Za-z.]+$")
#   ifelse(!is.na(refseq_style), refseq_style, sub(".*_", "", x))
# }
# 
# db_cluster <- db_cluster %>%
#   mutate(gene = trimws(gene)) %>%
#   mutate(gene_revised = ifelse(tool %in% c("resfinder"), extract_resfinder_accession(gene), gene_revised)) %>%
#   mutate(gene_revised = ifelse(tool == "resfinder" & grepl("^fosA", gene), sub("_.*", "", gene), gene_revised)) %>% 
#   mutate(gene_revised = ifelse(tool %in% c("resfinder") & gene_revised == "Y15705", "1_Y15705", gene_revised)) %>%
#   mutate(gene_revised = ifelse(tool %in% c("resfinder") & gene_revised == "Y15704", "1_Y15704", gene_revised))
# 
# db_cluster <- db_cluster %>%
#   mutate(gene_revised = ifelse(tool %in% c("resfinder"), trimws(extract_resfinder_accession(gene)), gene_revised))
# 
# db_cluster$gene_revised[db_cluster$tool == "resfinder" & db_cluster$gene_revised == "Y15704"] <- "1_Y15704"
# db_cluster$gene_revised[db_cluster$tool == "resfinder" & db_cluster$gene_revised == "Y15705"] <- "1_Y15705"
# 
# db_cluster <- db_cluster %>%
#   mutate(gene_revised = ifelse(tool %in% "amrfinderplus", 
#                                sapply(strsplit(gene, split = "|", fixed = TRUE), function(x) x[2]),
#                                gene_revised))
# 
# db_cluster <- db_cluster %>%
#   mutate(gene_revised = ifelse(tool %in% "deeparg", gene, gene_revised))
# db_cluster <- db_cluster %>% mutate(gene_revised = ifelse(tool %in% "fargene", norm_prime(gene), gene_revised))
# db_cluster$tool <- as.vector(tag_to_tool[db_cluster$tool])
# 
# 
# write.csv(db_cluster, file = "code_R_analysis/data_to_Zenodo/reference_database_clustering.csv", row.names = F)

