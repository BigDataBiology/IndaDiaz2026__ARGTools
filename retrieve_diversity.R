library(dplyr)
library(tidyr)



lst <- readRDS( "~/GitHub/arg_compare/results_tools.rds")

################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
# diversity 
# The abundances for each unigene had already been filterd with the file genes_prot_dna.csv

args_abundances <- read.delim("~/GitHub/arg_compare/data/abundances/args_abundances.tsv")
df2 <- readRDS( "~/df2.rds")

sum(unique(args_abundances$X) %in% df2$Term_ID)
sum(!unique(args_abundances$X) %in% df2$Term_ID)

# rgi.diamond
d <- lst$rgi.diamond %>% filter(!is.na(parent)) 
args_abundances <- args_abundances %>% mutate(rgi.diamond.ARO = d$ARO[match(X, d$query)], rgi.diamond.parent = d$parent[match(X, d$query)])

# rgi.diamond.prot
d <- lst$rgi.diamond.prot %>% filter(!is.na(parent)) 
args_abundances <- args_abundances %>% mutate(rgi.diamond.prot.ARO = d$ARO[match(X, d$query)], rgi.diamond.prot.parent = d$parent[match(X, d$query)])

# rgi.blast
d <- lst$rgi.blast %>% filter(!is.na(parent)) 
args_abundances <- args_abundances %>% mutate(rgi.blast.ARO = d$ARO[match(X, d$query)], rgi.blast.parent = d$parent[match(X, d$query)])

# deeparg
d <- lst$deeparg.norm %>% filter(!is.na(parent)) 
args_abundances <- args_abundances %>% mutate(deeparg.norm.ARO = d$ARO[match(X, d$query)], deeparg.norm.parent = d$parent[match(X, d$query)])

d <- lst$deeparg.norm.prot %>% filter(!is.na(parent)) 
args_abundances <- args_abundances %>% mutate(deeparg.norm.prot.ARO = d$ARO[match(X, d$query)], deeparg.norm.prot.parent = d$parent[match(X, d$query)])

# fargene
d <- lst$fargene %>% filter(!is.na(parent)) 
args_abundances <- args_abundances %>% mutate(fargene.ARO = d$ARO[match(X, d$query)], fargene.parent = d$parent[match(X, d$query)])

d <- lst$fargene.prot %>% filter(!is.na(parent)) 
args_abundances <- args_abundances %>% mutate(fargene.prot.ARO = d$ARO[match(X, d$query)], fargene.prot.parent = d$parent[match(X, d$query)])


# amrfinder
d <- lst$amrfinder.norm %>% filter(!is.na(parent))
args_abundances <- args_abundances %>% mutate(amrfinder.norm.ARO = d$ARO[match(X, d$query)], amrfinder.norm.parent = d$parent[match(X, d$query)])

d <- lst$amrfinder.norm.prot %>% filter(!is.na(parent))
args_abundances <- args_abundances %>% mutate(amrfinder.norm.prot.ARO = d$ARO[match(X, d$query)], amrfinder.norm.prot.parent = d$parent[match(X, d$query)])

gc()
# abricate 
d <- lst$abricate.argannot.norm %>% filter(!is.na(parent)) 
args_abundances <- args_abundances %>% mutate(abricate.argannot.norm.ARO = d$ARO[match(X, d$query)], abricate.argannot.norm.parent = d$parent[match(X, d$query)])

d <- lst$abricate.card.norm %>% filter(!is.na(parent))
args_abundances <- args_abundances %>% mutate(abricate.card.norm.ARO = d$ARO[match(X, d$query)], abricate.card.norm.parent = d$parent[match(X, d$query)])

d <- lst$abricate.megares.norm %>% filter(!is.na(parent)) 
args_abundances <- args_abundances %>% mutate(abricate.megares.norm.ARO = d$ARO[match(X, d$query)], abricate.megares.norm.parent = d$parent[match(X, d$query)])

d <- lst$abricate.ncbi.norm %>% filter(!is.na(parent)) 
args_abundances <- args_abundances %>% mutate(abricate.ncbi.norm.ARO = d$ARO[match(X, d$query)], abricate.ncbi.norm.parent = d$parent[match(X, d$query)])

d <- lst$abricate.resfinder.norm %>% filter(!is.na(parent)) 
args_abundances <- args_abundances %>% mutate(abricate.resfinder.norm.ARO = d$ARO[match(X, d$query)], abricate.resfinder.norm.parent = d$parent[match(X, d$query)])

# resfinder
d <- lst$resfinder.norm %>% filter(!is.na(parent)) 
args_abundances <- args_abundances %>% mutate(resfinder.norm.ARO = d$ARO[match(X, d$query)], resfinder.norm.parent = d$parent[match(X, d$query)])

gc()

#saveRDS(args_abundances, file = "abundances_per_tool_and_unigene.rds", compress = T)

diversity_parent <- bind_rows(args_abundances %>% group_by(sample, rgi.diamond.parent) %>% 
                                summarise(unigenes = n_distinct(X)) %>% 
                                mutate(tool = "rgi.diamond", parent = rgi.diamond.parent) %>% ungroup() %>% select(sample, parent, tool, unigenes),
                              
                              args_abundances %>% group_by(sample, rgi.diamond.prot.parent) %>% 
                                summarise(unigenes = n_distinct(X)) %>% 
                                mutate(tool = "rgi.diamond.prot", parent = rgi.diamond.prot.parent) %>% ungroup() %>% select(sample, parent, tool, unigenes),
                              
                              args_abundances %>% group_by(sample, rgi.blast.parent) %>% 
                                summarise(unigenes = n_distinct(X)) %>% 
                                mutate(tool = "rgi.blast", parent = rgi.blast.parent) %>% ungroup() %>% select(sample, parent, tool, unigenes),
                              
                              args_abundances %>% group_by(sample, deeparg.norm.parent) %>% 
                                summarise(unigenes = n_distinct(X)) %>% 
                                mutate(tool = "deeparg", parent = deeparg.norm.parent) %>% ungroup() %>% select(sample, parent, tool, unigenes),
                              
                              args_abundances %>% group_by(sample, deeparg.norm.prot.parent) %>% 
                                summarise(unigenes = n_distinct(X)) %>% 
                                mutate(tool = "deeparg.prot", parent = deeparg.norm.prot.parent) %>% ungroup() %>% select(sample, parent, tool, unigenes),
                              
                              args_abundances %>% group_by(sample, fargene.parent) %>% 
                                summarise(unigenes = n_distinct(X)) %>% 
                                mutate(tool = "fargene", parent = fargene.parent) %>% ungroup() %>% select(sample, parent, tool, unigenes),
                              
                              args_abundances %>% group_by(sample, fargene.prot.parent) %>% 
                                summarise(unigenes = n_distinct(X)) %>% 
                                mutate(tool = "fargene.prot", parent = fargene.prot.parent) %>% ungroup() %>% select(sample, parent, tool, unigenes),
                              
                              args_abundances %>% group_by(sample, amrfinder.norm.parent) %>% 
                                summarise(unigenes = n_distinct(X)) %>% 
                                mutate(tool = "amrfinder", parent = amrfinder.norm.parent) %>% ungroup() %>% select(sample, parent, tool, unigenes),
                              
                              args_abundances %>% group_by(sample, amrfinder.norm.prot.parent) %>% 
                                summarise(unigenes = n_distinct(X)) %>% 
                                mutate(tool = "amrfinder.prot", parent = amrfinder.norm.prot.parent) %>% ungroup() %>% select(sample, parent, tool, unigenes),
                              
                              args_abundances %>% group_by(sample, abricate.argannot.norm.parent) %>% 
                                summarise(unigenes = n_distinct(X)) %>% 
                                mutate(tool = "abricate.argannot", parent = abricate.argannot.norm.parent) %>% ungroup() %>% select(sample, parent, tool, unigenes),
                              
                              args_abundances %>% group_by(sample, abricate.card.norm.parent) %>% 
                                summarise(unigenes = n_distinct(X)) %>% 
                                mutate(tool = "abricate.card", parent = abricate.card.norm.parent) %>% ungroup() %>% select(sample, parent, tool, unigenes),
                              
                              args_abundances %>% group_by(sample, abricate.megares.norm.parent) %>% 
                                summarise(unigenes = n_distinct(X)) %>% 
                                mutate(tool = "abricate.megares", parent = abricate.megares.norm.parent) %>% ungroup() %>% select(sample, parent, tool, unigenes),
                              
                              args_abundances %>% group_by(sample, abricate.ncbi.norm.parent) %>% 
                                summarise(unigenes = n_distinct(X)) %>% 
                                mutate(tool = "abricate.ncbi", parent = abricate.ncbi.norm.parent) %>% ungroup() %>% select(sample, parent, tool, unigenes),
                              
                              args_abundances %>% group_by(sample, abricate.resfinder.norm.parent) %>% 
                                summarise(unigenes = n_distinct(X)) %>% 
                                mutate(tool = "abricate.resfinder", parent = abricate.resfinder.norm.parent) %>% ungroup() %>% select(sample, parent, tool, unigenes),
                              
                              args_abundances %>% group_by(sample, resfinder.norm.parent) %>% 
                                summarise(unigenes = n_distinct(X)) %>% 
                                mutate(tool = "resfinder", parent = resfinder.norm.parent) %>% ungroup() %>% select(sample, parent, tool, unigenes))

gc()

metadata <- read.delim("metadata_GMGC10.sample.meta.tsv")
diversity_parent <- diversity_parent %>% filter(!is.na(parent)) %>% 
  mutate(habitat = metadata$habitat[match(sample, metadata$sample_id)])

saveRDS(diversity_parent, file = "~/GitHub/arg_compare/diversity_per_tool.rds", compress = T)





