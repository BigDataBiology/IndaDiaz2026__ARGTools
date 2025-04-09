library(dplyr)
library(tidyr)
library(stringr)
library(rdflib)
library(stringr)
library(jsonlite)


################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
# wanted ontologies from OLS
# manually curated to match several tools 

# setwd("GitHub/arg_compare/")

higher_ontology <- paste0("ARO:", c("3000557", "0000010", "3000519", "3000185", "3000381", "3000159", "3000012"))

lowest_ontology <- c(paste0("ARO:", c("3000341", "3000322", "3000345", "3004257", "3007419", "3004276", "3004275", "3000229", "3000225", "3000228",
                                      "3000128", "3000127", "3000126", "3000155", "3000151", "3000154", "3000153", "3004260", "3000078", "3000004", 
                                      "3000076", "3000075", "3007074", "3000122", "3000249", "3004467", "3004064", "3000342", "3003025", "3000221", 
                                      "3000320", "3000458", "3000333", "3007103", "3000576", "3000233", "3000869", "3000036", "3004261", "3000231", 
                                      "3000234", "3007428", "0000031", "3000560", "3004469", "3000419", "3000507", "3005086", "0000002", "3003425", 
                                      "3001208", "3000210", "3004238", "3003040", "0010002", "3003580", "3003768", "3001207", "3000492", "3007610", 
                                      "3000100", "3007429", "3000270", "3000451", "3002976", "3007425", "3004916")))

ontologies <- c(lowest_ontology, higher_ontology)
# description of each ARO
# antibiotic inactivation enzyme ARO:3000557 NOT IN 
#   # AAC(2') ARO:3000341 (alternatively only ARO:3000121)
#   # AAC(3) ARO:3000322 (alternatively only ARO:3000121)
#   # AAC(6') ARO:3000345 (alternatively only ARO:3000121)
#   # cpa ARO:3004257
#   # aminoglycoside bifunctional ARO:3007419
#   # ant2 ARO:3004276 (alternatively only ARO:3000218)
#    # ant3 ARO:3004275 (alternatively only ARO:3000218)
#   # ant4 ARO:3000229 (alternatively only ARO:3000218)
#   # ant6 ARO:3000225 (alternatively only ARO:3000218)
#   # ant9 ARO:3000228 (alternatively only ARO:3000218)
#   # aph2'' ARO:3000128 (alternatively only ARO:3000114)
#   # aph3'' ARO:3000127 (alternatively only ARO:3000114)
#   # aph3'  ARO:3000126 (alternatively only ARO:3000114)
#   # aph4   ARO:3000155 (alternatively only ARO:3000114)
#   # aph6   ARO:3000151 (alternatively only ARO:3000114)
#   # aph7'' ARO:3000154 (alternatively only ARO:3000114)
#   # aph9   ARO:3000153 (alternatively only ARO:3000114)
#   # Bah ARO:3004260
#   # beta-lactam A ARO:3000078    (alternatively only ARO:3000001)
#   # beta-lactam B ARO:3000004    (alternatively only ARO:3000001)
#   # beta-lactam C ARO:3000076    (alternatively only ARO:3000001)
#   # beta-lactam D ARO:3000075    (alternatively only ARO:3000001)
#   # capreomicyn ARO:3007074
#   # CAT cloranphenicol acetyltransferase ARO:3000122
#   # cloranphenicol phosphotransferase ARO:3000249
#   # ciprofloxacin phosphotransferase ARO:3004467
#   # Edeine acetyltransferase ARO:3004064
#   # fosfomycin inactivation enzyme ARO:3000342
#   # fusidic acid inactivation enzyme ARO:3003025
#   # lincosamide nucleotidyltransferase ARO:3000221
#   # macrolide esterase ARO:3000320
#   # macrolide glycosyltransferase ARO:3000458
#   # macrolide phosphotransferase (MPH) ARO:3000333
#   # nitroimidazole reductase ARO:3007103
#   # rifampin inactivation enzyme ARO:3000576
#   # streptogramin inactivation enzyme ARO:3000233
#   # streptothricin acetyltransferase (SAT) ARO:3000869
#   # tetracycline inactivation enzyme ARO:3000036
#   # viomycin phosphotransferase ARO:3004261
# antibiotic resistance gene cluster, cassette, or operon ARO:0000010 NOT IN 
#   # beta-lactam resistance operon ARO:3000231
#   # glycopeptide resistance gene cluster ARO:3000234
#    # polymyxin resistance operon ARO:3007428
# antibiotic resistant gene variant or mutant ARO:0000031
# antibiotic target modifying enzyme ARO:3000519 NOT IN 
#    # Erm 23S ribosomal RNA methyltransferase ARO:3000560
# antibiotic target protection protein ARO:3000185 NOT IN
#   # ABC-F ATP-binding cassette ribosomal protection protein ARO:3004469
#   # quinolone resistance protein (qnr)  ARO:3000419
#   # rifampin-resistant RNA polymerase-binding protein ARO:3000507
#   # Target protecting FusB-type protein conferring resistance to Fusidic acid ARO:3005086
#   # tetracycline-resistant ribosomal protection protein ARO:0000002
# antibiotic target replacement protein ARO:3000381 NOT IN 
#   # antibiotic resistant dihydrofolate reductase                    ARO:3003425
#   # methicillin resistant PBP2                                      ARO:3001208
#   # rifamycin-resistant beta-subunit of RNA polymerase (rpoB)        ARO:3000210
#   # sulfonamide resistant sul                                        ARO:3004238
# beta-lactam resistant penicillin-binding proteins                ARO:3003040
# efflux pump complex or subunit conferring antibiotic resistance ARO:3000159 NOT IN
    # TET EFFLUX PUMP GOES HERE! major facilitator superfamily (MFS) antibiotic efflux pump ARO:0010002
# gene altering cell wall charge                                   ARO:3003580
# gene conferring resistance via absence                           ARO:3003768
# gene involved in antibiotic sequestration                        ARO:3001207
# gene involved in self-resistance to antibiotic                   ARO:3000492
# gene modulating azole resistance                                 ARO:3007610
# gene modulating beta-lactam resistance                           ARO:3000100
# gene(s) or protein(s) associated with polymyxin resistance operon ARO:3007429
# protein modulating permeability to antibiotic                    ARO:3000270
# protein(s) and two-component regulatory system modulating antibiotic efflux ARO:3000451
# protein(s) conferring antibiotic resistance via molecular bypass ARO:3000012 NOT IN 
    # gene(s) or protein(s) associated with a glycopeptide resistance cluster ARO:3002976
# protein(s) conferring resistance via host-dependent nutrient acquisition ARO:3007425
# subunits of secretion system conferring antibiotic resistance  ARO:3004916



################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
# colors for plots

pal <- c("#8c510a","#bf812d","#dfc27d","#f6e8c3","#c7eae5","#80cdc1","#35978f","#01665e")

################################################################################################################################################
################################################################################################################################################
################################################################################################################################################

## FARGENE CONVERSIONS


fg_class <- c("aac2p","aac6p","aac6p","aac3", "aph3p","aph6p","class_a","aac3", "mph",
             "class_b1_b2", "class_b3","class_d","aac6p","erm","class_d","class_c",
             "aph2b","erm","tet_enzyme", "tet_rpg", "tet_efflux", "qnr")

#  Since there are different HMM's per class, we reduce them to just the class 
names(fg_class) <- c("aac2p","aac6p_1","aac6p_2","aac3_2", "aph3p","aph6p","class_a","aac3_1", "mph",
                    "class_b1_b2", "class_b3","class_d1","aac6p_3","erm_2","class_d2","class_c",
                    "aph2b","erm_1","tet_enzyme", "tet_rpg", "tet_efflux", "qnr")

hmm_models <- c("aac2p", "aac3_1", "aac3_2", "aac6p_1", "aac6p_2", "aac6p_3", "aph2b", "aph3p", "aph6p", "class_b1_b2", "class_b3",
               "class_c", "class_a", "tet_efflux", "tet_enzyme", 
               "mph", "erm_1", "erm_2", "class_d1","class_d2", "tet_rpg", "qnr")

names(hmm_models) <- c("aac2p-aligned", "aac3_class1-aligned", "aac3_class2-aligned", "aac6p_class1-aligned", "aac6p_class2-aligned",
"aac6p_class3-aligned", "aph2b-aligned", "aph3p-aligned", "aph6-aligned", "b1_b2_70_centroids-aligned", "b3_70_centroids-aligned", 
"class_C_70_centroids-aligned", "classA_70_centroids-aligned", "efflux_model_group_1-aligned", "enzyme_reduced_tetX1_X3-aligned",
"macrolide_phosphotransferases-aligned", "methyltransferase_grp1-aligned", "methyltransferase_grp2-aligned", 
"oxa_g1_70_centroids-aligned", "oxa_g2_70_centroids-aligned", "rpg_reference_sequences-aligned", "pmqnr_20120719.pfa")

# The fargene classes corresponding to CARD aros 
fargene2ARO <- c("ARO:3000341", "ARO:3000322", "ARO:3000345", "ARO:3000128", "ARO:3000126", "ARO:3000151",
                 "ARO:3000078", "ARO:3000004", "ARO:3000004", "ARO:3000076", "ARO:3000075",
                 "ARO:3000560", "ARO:3000333", "ARO:3000036", "ARO:0000002", "ARO:0010002", "ARO:3000419")

names(fargene2ARO) <- c("aac2p", "aac3", "aac6p", "aph2b", "aph3p","aph6",
                        "class_a", "class_b1_b2", "class_b3", "class_c", "class_d", 
                        "erm", "mph", "tet_enzyme", "tet_rpg", "tet_efflux", "qnr")

################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
## FETCH THE ONTOLOGY
## We use this to simplify the ontology and summarise the abundance per class and not per unigene
## For each ARO identifided in any of the tools, we will assign it's lowest_ontology, 
# and if the ARO is higher than the lowest_ontology, we will assign higher_ontology 

aro_url <- "https://raw.githubusercontent.com/arpcard/aro/master/aro.owl"
g <- rdf_parse(aro_url, format = "rdfxml")

# Define necessary namespaces
rdfs_ns <- "http://www.w3.org/2000/01/rdf-schema#"

# **STEP 1: Preload Parent-Child Relationships**
query_parents <- paste0("PREFIX rdfs: <", rdfs_ns, "> ",
                        "SELECT ?child ?parent WHERE { ?child rdfs:subClassOf ?parent }")
parent_map <- rdf_query(g, query_parents)  # DataFrame with columns "child" and "parent"

# Convert URIs to ARO:XXXX format
parent_map <- parent_map %>%
  mutate(Child_ID = str_extract(child, "ARO_[0-9]+") %>% str_replace_all("_", ":"),
         Parent_ID = str_extract(parent, "ARO_[0-9]+") %>% str_replace_all("_", ":")) %>%
  select(Child_ID, Parent_ID) %>%
  filter(!is.na(Child_ID) & !is.na(Parent_ID))  # Remove NA values

# **STEP 2: Preload Labels**
query_labels <- paste0("PREFIX rdfs: <", rdfs_ns, "> ",
                       "SELECT ?term ?label WHERE { ?term rdfs:label ?label }")
label_map <- rdf_query(g, query_labels)
# Convert URIs to ARO format
label_map <- label_map %>%
  mutate(Term_ID = str_extract(term, "ARO_[0-9]+") %>% str_replace_all("_", ":"),
         Term_Label = label) %>%
  select(Term_ID, Term_Label) %>%
  filter(!is.na(Term_ID))

get_parents_fast <- function(term_id, parent_map) {
  parents <- parent_map %>% filter(Child_ID == term_id) %>% pull(Parent_ID)
  all_parents <- parents
  while (length(parents) > 0) {
    parents <- parent_map %>% filter(Child_ID %in% parents) %>% pull(Parent_ID)
    all_parents <- unique(c(all_parents, parents))
  }
  return(all_parents)
}

process_terms_fast <- function(term_list, parent_map, label_map) {
  results <- data.frame(Term_ID = character(), Term_Label = character(),
                        Parent_ID = character(), Parent_Label = character(),
                        stringsAsFactors = FALSE)
  for (term_id in term_list) {
    term_label <- label_map %>% filter(Term_ID == term_id) %>% pull(Term_Label)
    term_label <- ifelse(length(term_label) > 0, term_label, "Unknown")
    parents <- get_parents_fast(term_id, parent_map)
    for (parent_id in parents) {
      parent_label <- label_map %>% filter(Term_ID == parent_id) %>% pull(Term_Label)
      parent_label <- ifelse(length(parent_label) > 0, parent_label, "Unknown")
      results <- rbind(results, data.frame(Term_ID = term_id, Term_Label = term_label,
                                           Parent_ID = parent_id, Parent_Label = parent_label,
                                           stringsAsFactors = FALSE))
    }
  }
  return(results)
}


################################################################################################################################################
################################################################################################################################################
################################################################################################################################################

# read the results per tool
# amrfinder

amrfinder.norm <- read.delim("dna/amrfinder.norm.tsv") %>% 
  rename(query = Contig.id) %>% 
  rename(ARG.class = Class) %>% 
  rename(ARG.subclass = Subclass)

amrfinder.norm.prot <- read.delim("protein/amrfinder.norm.tsv") %>% 
  rename(query = Protein.id) %>% 
  rename(ARG.class = Class) %>% 
  rename(ARG.subclass = Subclass)

# deeparg
deeparg.norm <- read.delim("dna/deeparg.norm.tsv") %>% 
  rename(query = read_id) %>% 
  rename(ARG.class = predicted_ARG.class)

deeparg.norm.prot <- read.delim("protein/deeparg.norm.tsv") %>% 
  rename(query = read_id) %>% 
  rename(ARG.class = predicted_ARG.class)

# rgi
rgi.diamond <- read.delim("dna/rgi_diamond.tsv")
rgi.diamond <- rgi.diamond[,-c(1,18,19,20)]
rgi.diamond$query <- gsub('.{2}$', '', rgi.diamond$Contig)
rgi.diamond <- rgi.diamond %>% rename(ARG.class = AMR.Gene.Family) %>%
  mutate(ARO = paste0("ARO:", ARO))

rgi.blast <- read.delim("dna/rgi_blast.tsv")
rgi.blast <- rgi.blast[,-c(1,18,19,20)]
rgi.blast$query <- gsub('.{2}$', '', rgi.blast$Contig)  
rgi.blast <- rgi.blast %>% rename(ARG.class = AMR.Gene.Family) %>%
  mutate(ARO = paste0("ARO:", ARO))

rgi.diamond.prot <- read.delim("protein/rgi_diamond.tsv")
rgi.diamond.prot <- rgi.diamond.prot[,-c(18,19,20)]
rgi.diamond.prot$query <- rgi.diamond.prot$ORF_ID
rgi.diamond.prot <- rgi.diamond.prot %>% rename(ARG.class = AMR.Gene.Family) %>%
  mutate(ARO = paste0("ARO:", ARO))

# fargene
fargene <- read.delim("dna/fargene_results.tsv", header = F)
fargene$query <- gsub('.{2}$', '', fargene$V1)
fargene$new_class = fg_class[fargene$V5]

## remove duplicated queries with different classes 
d1 <- duplicated(paste(fargene$query, fargene$new_class))
d2 <- duplicated(fargene$query)
number_remove_fargene <- sum(!(!fargene$query %in% fargene$query[!fargene$query %in% fargene$query[d1] & fargene$query %in% fargene$query[d2]]))
fargene <- fargene[!fargene$query %in% fargene$query[!fargene$query %in% fargene$query[d1] & fargene$query %in% fargene$query[d2]],]

## remove duplicated queries with same class
fargene <- fargene[!duplicated(fargene$query),]
number_remove_fargene_same_class <- sum(duplicated(fargene$query))
rm(d1, d2)

# load the hmm scores
hmm <- read.table("dna/fargene_hmm.txt", quote="\"", comment.char="")
hmm$query <- gsub('.{2}$', '', hmm$V1)
hmm$q1 <- sapply(strsplit(hmm$V1, split = "GMGC10"), function(x) paste0("GMGC10",x[length(x)]))
hmm$new_class <- hmm_models[hmm$V4]
hmm$q1 <- gsub('.{2}$', '', hmm$q1)
hmm$q1 <- sapply(strsplit(hmm$q1, split = "_seq"), function(x) x[1])
hmm <- hmm[order(hmm$V14, decreasing = T),]
fargene$hmm <- hmm$V14[match(paste(fargene$query, fargene$V5), paste(hmm$q1, hmm$new_class))]
rm(hmm)
#

fargene.prot <- read.delim("protein/fargene_results.tsv", header = F)
fargene.prot$query <- fargene.prot$V1 
fargene.prot$new_class = fg_class[fargene.prot$V5]
# load the HMM scores 
hmm.prot <- read.table("protein/fargene_hmm.txt", quote="\"", comment.char="")
hmm.prot$new_class = hmm_models[hmm.prot$V4]
hmm.prot$q1 = hmm.prot$V1
hmm.prot <- hmm.prot[order(hmm.prot$V14, decreasing = T),]
fargene.prot$hmm <- hmm.prot$V14[match(paste(fargene.prot$V1, fargene.prot$V5), paste(hmm.prot$q1, hmm.prot$new_class))]

number_remove_fargene.prot_same_class <- sum(duplicated(fargene.prot$V1))
fargene.prot <- fargene.prot[!duplicated(fargene.prot$V1),]
rm(hmm.prot)

fargene$new_class[fargene$new_class %in% "aph6p"] <- "aph6"
fargene.prot$new_class[fargene.prot$new_class %in% "aph6p"] <- "aph6"

# abricate
abricate.argannot.norm <- read.delim("dna/abricate-argannot.norm.tsv", header=FALSE, comment.char="#")
abricate.argannot.norm <- abricate.argannot.norm[,-1]
abricate.argannot.norm <- abricate.argannot.norm %>% rename(query = V2, ARO = V16)
abricate.argannot.norm$ARG.class <- str_match(abricate.argannot.norm$V14, "\\(([^)]+)\\)")[,2]
abricate.argannot.norm$gene <- sub("\\([^)]*\\)", "", abricate.argannot.norm$V14)

abricate.card.norm <- read.delim("dna/abricate-card.tsv", header=FALSE, comment.char="#")
abricate.card.norm <- abricate.card.norm[,-1]
abricate.card.norm <- abricate.card.norm %>% rename(query = V2, ARG.class = V15)

abricate.megares.norm <- read.delim("dna/abricate-megares.norm.tsv", header=FALSE, comment.char="#")
abricate.megares.norm <- abricate.megares.norm[,-1] %>% rename(query = V2, ARO = V16)

abricate.ncbi.norm <- read.delim("dna/abricate-ncbi.norm.tsv", header=FALSE, comment.char="#")
abricate.ncbi.norm <- abricate.ncbi.norm[,-1] %>% rename(query = V2) %>% rename(ARG.class = V15, gene = V14, ARO = V16)

abricate.resfinder.norm <- read.delim("dna/abricate-resfinder.norm.tsv", header=FALSE, comment.char="#")
abricate.resfinder.norm <- abricate.resfinder.norm[,-1] %>% rename(query = V2, drug = V15, gene = V14, ARO = V16)

# resfinder
resfinder.norm <- read.delim("dna/resfinder_table.norm.tsv")
resfinder.norm <- resfinder.norm %>% rename(query = Contig) %>% rename(gene = Resistance.gene) 


################################################################################################################################################
# complement resfinder with its phenotype (from /work/microbiome/users/juan/resfinder_databases/resfinder_db/phenotypes.txt)

phenotypes_resfiner <- read.delim("check_missing_annot/phenotypes.txt")
phenotypes_resfiner <- phenotypes_resfiner[!duplicated(phenotypes_resfiner$Gene_accession.no.),]
v <- vector()
for(j in 1:nrow(resfinder.norm)){
  sect <- intersect(grep(resfinder.norm$gene[j], phenotypes_resfiner$Gene_accession.no., fixed = TRUE),
                    grep(resfinder.norm$Accession.no.[j], phenotypes_resfiner$Gene_accession.no., fixed = TRUE))
  if(length(sect)>1) {
    print(j)
    print(sect)}
  if(3122 %in% sect){
    v <- c(v,3122)
  } else {
    v <- c(v,sect[1])
  }
}

resfinder.norm$ARG.class <- phenotypes_resfiner$Class[v]
rm(phenotypes_resfiner, v)


################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
# ARO for abricate card
# Innformation from CARD models (ex. aro identifiers, bit socores per ARO, etc.)
# It could be useful to obtain the AROs for abricate.CARD 
# CARD JSON


card_data <- fromJSON("check_missing_annot/card.json")


card_abricate_accession <- unique(sapply(strsplit(abricate.card.norm$V13, split = ":"), function(x) x[1]))
card_abricate_gene <- unique(abricate.card.norm$V6)

# Retrieve the aro name and aro accession for all entries in the json file
card_gene_aro <- data.frame(gene=NULL, aro=NULL)
for(j in 1:(length(card_data)-3)){
  if(card_data[[j]]$ARO_name %in%  card_abricate_gene){
    card_gene_aro <- rbind(card_gene_aro, data.frame(gene=card_data[[j]]$ARO_name, aro=card_data[[j]]$ARO_accession))  
  } else {
    if(card_data[[j]]$CARD_short_name %in%  card_abricate_gene){
      card_gene_aro <- rbind(card_gene_aro, data.frame(gene=card_data[[j]]$CARD_short_name, aro=card_data[[j]]$ARO_accession))  
    } else {
      if(card_data[[j]]$model_name %in%  card_abricate_gene){
        card_gene_aro <- rbind(card_gene_aro, data.frame(gene=card_data[[j]]$model_name, aro=card_data[[j]]$ARO_accession))  
      }
    }
  }
}

rm(card_data)

# load the manually curated AROS
manually_curated <- read.table("check_missing_annot/abricate-card-manual-curation.txt", stringsAsFactors = F, text = T)
manually_curated <- data.frame(gene = sapply(strsplit(manually_curated$V1, split = ","), function(x) x[1]), aro = sapply(strsplit(manually_curated$V1, split = ","), function(x) x[2]))

# add them to the data retrieved from json file
card_gene_aro <- rbind(card_gene_aro, manually_curated)

# add the ARO to abricate.card.norm
abricate.card.norm$ARO <- paste0("ARO:",card_gene_aro$aro[match(abricate.card.norm$V6, card_gene_aro$gene)])


################################################################################################################################################
# ARO for abricate.argannot.norm
# manually assign the following 
# abricate.argannot.norm[abricate.argannot.norm$ARO == "",]
abricate.argannot.norm$ARO[abricate.argannot.norm$V6 == "(Ntmdz)nimj_Nitroimidazole_Gene"] = "ARO:3007112"
abricate.argannot.norm$ARO[abricate.argannot.norm$V6 == "(Bla)blaACT-16"] = "ARO:3001827"

################################################################################################################################################
# ARO for abricate.resfinder.norm
resfinder_ARO_mapping <- read.delim("check_missing_annot/resfinder_ARO_mapping.tsv")

# correct the capital letters in these two entries
abricate.resfinder.norm$V6[abricate.resfinder.norm$V6 == "oqxA_1"] <- "OqxA_1"
abricate.resfinder.norm$V6[abricate.resfinder.norm$V6 == "oqxB_1"] <- "OqxB_1"

j <- which(abricate.resfinder.norm$ARO == "")

resfinder_missing_aro <- tibble(gene=abricate.resfinder.norm$V6[j], accession=abricate.resfinder.norm$V13[j]) %>% distinct() %>%
  filter(!is.na(match(paste(gene, accession, sep = "_"), resfinder_ARO_mapping$Original.ID)))

k <- match(paste(abricate.resfinder.norm$V6[j], abricate.resfinder.norm$V13[j], sep = "_"), resfinder_ARO_mapping$Original.ID)

abricate.resfinder.norm$ARO[j]  <- paste0("ARO:", resfinder_ARO_mapping$ARO[k])
abricate.resfinder.norm$ARO[abricate.resfinder.norm$ARO == "ARO:NA"] <- ""

rm(j, k)


################################################################################################################################################
# ARO for resfinder.norm

d <- unique(resfinder.norm$Accession.no.[resfinder.norm$ARO == ""])
d2 <- cbind(d, NA)
for(j in 1:length(d)){
  d2[j,2] <- paste0("ARO:", resfinder_ARO_mapping$ARO[grep(d[j], resfinder_ARO_mapping$Original.ID)])
}
d2 <- as.data.frame(d2)
d2$V2[d2$V2=="ARO:"] <- ""
j <- which(resfinder.norm$ARO == "")
resfinder.norm$ARO[j] <- d2$V2[match(resfinder.norm$Accession.no.[j], d2$d)]

# to debug arg_norm
resfinder_missing_aro_resfinder <- d2[d2$V2 != "",]
resfinder_missing_aro_resfinder$gene <- resfinder.norm$gene[match(resfinder_missing_aro_resfinder$d, resfinder.norm$Accession.no.)]
resfinder_missing_aro_resfinder <- resfinder_missing_aro_resfinder[,c(3,1)]
names(resfinder_missing_aro_resfinder) <- c("gene","accession")
rm(resfinder_ARO_mapping, j, d, d2)


################################################################################################################################################
# ARO for abricate megares

abricate_annotation <- read.delim("check_missing_annot/abricate_annotation_megares.tsv")
abricate_annotation <- abricate_annotation[,-c(19,20,21)]
abricate_annotation$MEG <- sapply(strsplit(abricate_annotation$ORF, split = "~~~"), function(x) x[3])
abricate_annotation$gene <- sapply(strsplit(abricate_annotation$ORF, split = "~~~"), function(x) x[2])
j <- which(abricate.megares.norm$ARO == "")
abricate.megares.norm$ARO[j] <- paste0("ARO:", abricate_annotation$ARO[match(abricate.megares.norm$V13[j], abricate_annotation$MEG)])
abricate.megares.norm$ARO[abricate.megares.norm$ARO == "ARO:NA"] <- ""

# to debug arg_norm
abricate_megares_missing_aro <- tibble(gene = abricate.megares.norm$V6[j], accession = abricate.megares.norm$V13[j], aro = abricate.megares.norm$ARO[j])
abricate_megares_missing_aro <- abricate_megares_missing_aro %>% filter(aro != "")
abricate_megares_missing_aro <- abricate_megares_missing_aro[,c(1,2)]
abricate_megares_missing_aro <- abricate_megares_missing_aro %>% distinct()

rm(j, abricate_annotation)

################################################################################################################################################
# ARO for abricate NCBI

abricate_annotation <- read.delim("check_missing_annot/abricate_annotation.tsv")
abricate_annotation <- abricate_annotation[,-c(19,20,21)]
x2 <- sapply(strsplit(abricate_annotation$ORF, split = "~~~"), function(x) x[2])
x3 <- sapply(strsplit(abricate_annotation$ORF, split = "~~~"), function(x) x[3])

j <- abricate.ncbi.norm$ARO ==""
d <- paste0("ARO:", abricate_annotation$ARO[match(abricate.ncbi.norm$V6[j], x2)])
d[d=="ARO:NA"] <- ""
abricate.ncbi.norm$ARO[j] <- d

# to debug arg_norm
abricate_ncbi_missing_aro <- tibble(gene = abricate.ncbi.norm$V6[j], accession = abricate.ncbi.norm$V13[j], aro = abricate.ncbi.norm$ARO[j])
abricate_ncbi_missing_aro <- abricate_ncbi_missing_aro %>% filter(aro != "")
abricate_ncbi_missing_aro <- abricate_ncbi_missing_aro[,c(1,2)]
abricate_ncbi_missing_aro <- abricate_ncbi_missing_aro %>% distinct()

rm(j, d, abricate_annotation, x2, x3)
# I could not find anything extra for amrfinder 

# to debug arg_norm
# write.csv(file = "retrieved_manual_rgi.csv", bind_rows(resfinder_missing_aro %>% mutate(tool = "abricate_resfinder"), 
# resfinder_missing_aro_resfinder %>% mutate(tool = "resfinder"),
# abricate_megares_missing_aro %>% mutate(tool = "abricate_megares"),
# abricate_ncbi_missing_aro %>% mutate(tool = "abricate_ncbi")), row.names = F)

rm(abricate_ncbi_missing_aro, abricate_megares_missing_aro, resfinder_missing_aro, resfinder_missing_aro_resfinder)



################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
### AROs

# Create a list of unique AROs in all tools 
# I need to include those from abricate.CARD at some point 
# I need to fix by hand those AROs missing in each tool

fargene$ARO <- as.vector(fargene2ARO[fargene$new_class])
fargene.prot$ARO <- as.vector(fargene2ARO[fargene.prot$new_class])

aros <- tibble(data.frame(rbind(cbind("RGI - diamond NT", rgi.diamond$Best_Hit_ARO, rgi.diamond$ARO),
                                cbind("RGI - diamond AA", rgi.diamond.prot$Best_Hit_ARO, rgi.diamond.prot$ARO),
                                cbind("RGI - blast NT", rgi.blast$Best_Hit_ARO, rgi.blast$ARO),
                                cbind("deepArg - NT", deeparg.norm$X.ARG, deeparg.norm$ARO),
                                cbind("deepArg - AA", deeparg.norm.prot$X.ARG, deeparg.norm.prot$ARO),
                                cbind("amrfinder NT", amrfinder.norm$Element.symbol, amrfinder.norm$ARO),
                                cbind("amrfinder - AA", amrfinder.norm.prot$Element.symbol, amrfinder.norm.prot$ARO),
                                cbind("abricate - ARGANNOT NT", abricate.argannot.norm$V6, abricate.argannot.norm$ARO),
                                cbind("abricate - MEGARES NT", abricate.megares.norm$V6, abricate.megares.norm$ARO),
                                cbind("abricate - NCBI NT", abricate.ncbi.norm$V6, abricate.ncbi.norm$ARO),
                                cbind("abricate - ResFinder NT", abricate.resfinder.norm$V6, abricate.resfinder.norm$ARO),
                                cbind("abricate - CARD NT", abricate.card.norm$V6, abricate.card.norm$ARO),
                                cbind("ResFinder - NT", resfinder.norm$gene, resfinder.norm$ARO),
                                cbind("fargene - NT", fargene$new_class, fargene$ARO),
                                cbind("fargene - AA", fargene.prot$new_class, fargene.prot$ARO))))

aros <- aros %>% distinct()
aros %>% group_by(X1) %>% summarise(n = n(), e = sum(X3 == "")) %>% mutate(p=e/n) %>% arrange(desc(n))

# Information shared with the group
# no_aros_summary <- aros %>% group_by(X1) %>% summarise(total_genes = n(), genes_with_no_ARO = sum(X3 == "")) %>% mutate(percentage_no_aro = round(genes_with_no_ARO/total_genes*100,1))
# df_no_ARO <- data.frame(aros %>% filter(X3=="")) %>% select(X1, X2)
# write.csv(df_no_ARO, file = "No_ARO_per_tool.csv", row.names = F, quote=F)
# rm(no_aros_summary, df_no_ARO)

# here we use the ontology information g, parent_map, label_map, rdfs_ns, aro_url, query_parents, parent_map
df <- process_terms_fast(unique(aros$X3[aros$X3!=""]), parent_map, label_map)
rm(g, rdfs_ns, aro_url, query_parents, parent_map, query_labels, label_map)


# THIS MIGHT NEED TO CHENGE IF MORE GENES APPEAR
# SOME GENES SEEM TO APPEAR IN TWO OR MORE DIFFERENT BRANCHES OF ONTOLOGY
# ARBITRARLY ASSIGNED ONE OF THE TWO OR MORE BRANCHES

# methicillin resistant PBP2  ARO:3001208 --->  pbp2 ARO:3003040
# Penicillin-binding protein mutations conferring resistance to beta-lactam antibiotics ARO:3003938 --->  pbp2 ARO:3003040
# gene(s) or protein(s) associated with polymyxin resistance operon ARO:3007429 ---> ARO:3003580
# ARO:0010004 resistance-nodulation-cell division (RND) antibiotic efflux pump removed ---> replaced by ARO:3000451
# ARO:3000219 mutant efflux regulatory protein conferring antibiotic resistance ---> replaced by ARO:3000451
# ARO:3000270 protein modulating permeability to antibiotic ---> ARO:3000451
# ARO:3000492  gene involved in self-resistance to antibiotic ---> ARO:3004064
# child ARO:3002484 gets ARO:3000076   || remove ARO:3000075 (classs D) for ARO:3000076 class C
# antibiotic resistant gene variant or mutant ARO:0000031 ---> ARO:3000451 and ARO:3003040

overrul <- c("ARO:3003040", "ARO:3003040", "ARO:3003580", "ARO:3000451", "ARO:3004064", "ARO:3000076", "ARO:3000451", "ARO:3003040", "ARO:3000004", "ARO:3000210")

# FIRST WE DEAL WITH THOSE AROS THAT ARE ALREADY IN THE LOWEST DESIRED ONTOLOGY
df_lowest_term <- df %>% filter(Term_ID %in% lowest_ontology) %>% mutate(Parent_ID = Term_ID, Parent_Label = Term_Label) %>% distinct()
length(unique(df_lowest_term$Term_ID))

# Now we deal with the AROS which parent is in the lowest ontology
df_lowest_parent <- df %>% filter(Parent_ID %in% lowest_ontology & !Term_ID %in% df_lowest_term$Term_ID) 
length(unique(df_lowest_parent$Term_ID))

# Now we deal with the AROS which parent is in the higher ontology
df_not_lowest <- df %>% filter(!Term_ID %in% c(df_lowest_parent$Term_ID, df_lowest_term$Term_ID))
length(unique(df_not_lowest$Term_ID))

# Fort lowest parent 
# this are the aros appearing in several branches
repeated_wanted <- df_lowest_parent %>% filter(Parent_ID %in% lowest_ontology) %>% group_by(Term_ID) %>% summarise(n = n()) %>% filter(n >1)
# those repeated after selecting overrul aro, if this had entries, we need to modify overrul
df_lowest_parent %>% filter(Parent_ID %in% lowest_ontology) %>%  group_by(Term_ID) %>% filter(n()>1)%>% group_by(Term_ID) %>% summarise(n=sum(Parent_ID %in% overrul)) %>% filter(n<1)
# unique lowest for df lowest parent 
df_lowest_parent_2 <- df_lowest_parent %>% filter(Parent_ID %in% lowest_ontology) %>% group_by(Term_ID) %>% group_by(Term_ID) %>% filter(n()<2) %>% ungroup()
# add those with double ontology
df_lowest_parent_2 <- df_lowest_parent_2 %>% bind_rows(df_lowest_parent %>% filter(Parent_ID %in% lowest_ontology) %>% group_by(Term_ID) %>% group_by(Term_ID) %>% filter(n()>1) %>% ungroup() %>% filter(Parent_ID %in% overrul))


# For not lowest
# this are the aros appearing in several branches
repeated_wanted <- df_not_lowest %>% filter(Parent_ID %in% ontologies) %>% group_by(Term_ID) %>% summarise(n = n()) %>% filter(n >1)
# those repeated after selecting overrul aro, if this had entries, we need to modify overrul
df_not_lowest %>% filter(Parent_ID %in% ontologies) %>%  group_by(Term_ID) %>% filter(n()>1)%>% group_by(Term_ID) %>% summarise(n=sum(Parent_ID %in% overrul)) %>% filter(n<1)
# unique lowest for df lowest parent 
df_not_lowest_2 <- df_not_lowest %>% filter(Parent_ID %in% ontologies) %>% group_by(Term_ID) %>% filter(n()<2) %>% ungroup()
# add those with double ontology
df_not_lowest_2 <- df_not_lowest_2 %>% bind_rows(df_not_lowest %>% filter(Parent_ID %in% ontologies) %>% group_by(Term_ID) %>% group_by(Term_ID) %>% filter(n()>1) %>% ungroup() %>% filter(Parent_ID %in% overrul))
dim(df_not_lowest_2)

df2 <- df_lowest_term %>% bind_rows(df_lowest_parent_2) %>% bind_rows(df_not_lowest_2)
rm(df_lowest_term, df_lowest_parent, df_not_lowest, df_lowest_parent_2, df_not_lowest_2, aros)

# Did we catch all of the aros?
sum(!df$Term_ID %in% df2$Term_ID)


### Complement tools 
deeparg.norm <- deeparg.norm %>% mutate(parent = df2$Parent_ID[match(ARO, df2$Term_ID)], parent_description = df2$Parent_Label[match(ARO, df2$Term_ID)])
deeparg.norm.prot  <- deeparg.norm.prot %>% mutate(parent = df2$Parent_ID[match(ARO, df2$Term_ID)], parent_description = df2$Parent_Label[match(ARO, df2$Term_ID)])
rgi.blast  <- rgi.blast %>% mutate(parent = df2$Parent_ID[match(ARO, df2$Term_ID)], parent_description = df2$Parent_Label[match(ARO, df2$Term_ID)])
rgi.diamond  <- rgi.diamond %>% mutate(parent = df2$Parent_ID[match(ARO, df2$Term_ID)], parent_description = df2$Parent_Label[match(ARO, df2$Term_ID)])
rgi.diamond.prot  <- rgi.diamond.prot %>% mutate(parent = df2$Parent_ID[match(ARO, df2$Term_ID)], parent_description = df2$Parent_Label[match(ARO, df2$Term_ID)])
fargene  <- fargene %>% mutate(parent = df2$Parent_ID[match(ARO, df2$Term_ID)], parent_description = df2$Parent_Label[match(ARO, df2$Term_ID)])
fargene.prot  <- fargene.prot %>% mutate(parent = df2$Parent_ID[match(ARO, df2$Term_ID)], parent_description = df2$Parent_Label[match(ARO, df2$Term_ID)])
amrfinder.norm  <- amrfinder.norm %>% mutate(parent = df2$Parent_ID[match(ARO, df2$Term_ID)], parent_description = df2$Parent_Label[match(ARO, df2$Term_ID)])
amrfinder.norm.prot  <- amrfinder.norm.prot %>% mutate(parent = df2$Parent_ID[match(ARO, df2$Term_ID)], parent_description = df2$Parent_Label[match(ARO, df2$Term_ID)])
abricate.argannot.norm  <- abricate.argannot.norm %>% mutate(parent = df2$Parent_ID[match(ARO, df2$Term_ID)], parent_description = df2$Parent_Label[match(ARO, df2$Term_ID)])
abricate.card.norm  <- abricate.card.norm %>% mutate(parent = df2$Parent_ID[match(ARO, df2$Term_ID)], parent_description = df2$Parent_Label[match(ARO, df2$Term_ID)])
abricate.megares.norm  <- abricate.megares.norm %>% mutate(parent = df2$Parent_ID[match(ARO, df2$Term_ID)], parent_description = df2$Parent_Label[match(ARO, df2$Term_ID)])
abricate.ncbi.norm  <- abricate.ncbi.norm %>% mutate(parent = df2$Parent_ID[match(ARO, df2$Term_ID)], parent_description = df2$Parent_Label[match(ARO, df2$Term_ID)])
abricate.resfinder.norm  <- abricate.resfinder.norm %>% mutate(parent = df2$Parent_ID[match(ARO, df2$Term_ID)], parent_description = df2$Parent_Label[match(ARO, df2$Term_ID)])
resfinder.norm  <- resfinder.norm %>% mutate(parent = df2$Parent_ID[match(ARO, df2$Term_ID)], parent_description = df2$Parent_Label[match(ARO, df2$Term_ID)])


saveRDS(list(deeparg.norm = deeparg.norm, deeparg.norm.prot = deeparg.norm.prot, 
     rgi.blast = rgi.blast, rgi.diamond = rgi.diamond, rgi.diamond.prot = rgi.diamond.prot,
     fargene = fargene, fargene.prot = fargene.prot, amrfinder.norm = amrfinder.norm, amrfinder.norm.prot = amrfinder.norm.prot,
     abricate.argannot.norm = abricate.argannot.norm, abricate.card.norm = abricate.card.norm, abricate.megares.norm= abricate.megares.norm,
     abricate.ncbi.norm = abricate.ncbi.norm, abricate.resfinder.norm = abricate.resfinder.norm, resfinder.norm = resfinder.norm), 
     file = "results_tools.rds", compress = T)

################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
# ABUNDANCES 
# The abundances for each unigene had already been filterd with the file genes_prot_dna.csv

args_abundances <- read.delim("data/abundances/args_abundances.tsv")
sum(unique(args_abundances$X) %in% df2$Term_ID)
sum(!unique(args_abundances$X) %in% df2$Term_ID)

# rgi.diamond
d <- rgi.diamond %>% filter(!is.na(parent)) 
args_abundances <- args_abundances %>% mutate(rgi.diamond.ARO = d$ARO[match(X, d$query)], rgi.diamond.parent = d$parent[match(X, d$query)])

# rgi.diamond.prot
d <- rgi.diamond.prot %>% filter(!is.na(parent)) 
args_abundances <- args_abundances %>% mutate(rgi.diamond.prot.ARO = d$ARO[match(X, d$query)], rgi.diamond.prot.parent = d$parent[match(X, d$query)])

# rgi.blast
d <- rgi.blast %>% filter(!is.na(parent)) 
args_abundances <- args_abundances %>% mutate(rgi.blast.ARO = d$ARO[match(X, d$query)], rgi.blast.parent = d$parent[match(X, d$query)])

# deeparg
d <- deeparg.norm %>% filter(!is.na(parent)) 
args_abundances <- args_abundances %>% mutate(deeparg.norm.ARO = d$ARO[match(X, d$query)], deeparg.norm.parent = d$parent[match(X, d$query)])

d <- deeparg.norm.prot %>% filter(!is.na(parent)) 
args_abundances <- args_abundances %>% mutate(deeparg.norm.prot.ARO = d$ARO[match(X, d$query)], deeparg.norm.prot.parent = d$parent[match(X, d$query)])

# fargene
d <- fargene %>% filter(!is.na(parent)) 
args_abundances <- args_abundances %>% mutate(fargene.ARO = d$ARO[match(X, d$query)], fargene.parent = d$parent[match(X, d$query)])

d <- fargene.prot %>% filter(!is.na(parent)) 
args_abundances <- args_abundances %>% mutate(fargene.prot.ARO = d$ARO[match(X, d$query)], fargene.prot.parent = d$parent[match(X, d$query)])


# amrfinder
d <- amrfinder.norm %>% filter(!is.na(parent))
args_abundances <- args_abundances %>% mutate(amrfinder.norm.ARO = d$ARO[match(X, d$query)], amrfinder.norm.parent = d$parent[match(X, d$query)])

d <- amrfinder.norm.prot %>% filter(!is.na(parent))
args_abundances <- args_abundances %>% mutate(amrfinder.norm.prot.ARO = d$ARO[match(X, d$query)], amrfinder.norm.prot.parent = d$parent[match(X, d$query)])

# abricate 
d <- abricate.argannot.norm %>% filter(!is.na(parent)) 
args_abundances <- args_abundances %>% mutate(abricate.argannot.norm.ARO = d$ARO[match(X, d$query)], abricate.argannot.norm.parent = d$parent[match(X, d$query)])

d <- abricate.card.norm %>% filter(!is.na(parent))
args_abundances <- args_abundances %>% mutate(abricate.card.norm.ARO = d$ARO[match(X, d$query)], abricate.card.norm.parent = d$parent[match(X, d$query)])

d <- abricate.megares.norm %>% filter(!is.na(parent)) 
args_abundances <- args_abundances %>% mutate(abricate.megares.norm.ARO = d$ARO[match(X, d$query)], abricate.megares.norm.parent = d$parent[match(X, d$query)])

d <- abricate.ncbi.norm %>% filter(!is.na(parent)) 
args_abundances <- args_abundances %>% mutate(abricate.ncbi.norm.ARO = d$ARO[match(X, d$query)], abricate.ncbi.norm.parent = d$parent[match(X, d$query)])

d <- abricate.resfinder.norm %>% filter(!is.na(parent)) 
args_abundances <- args_abundances %>% mutate(abricate.resfinder.norm.ARO = d$ARO[match(X, d$query)], abricate.resfinder.norm.parent = d$parent[match(X, d$query)])

# resfinder
d <- resfinder.norm %>% filter(!is.na(parent)) 
args_abundances <- args_abundances %>% mutate(resfinder.norm.ARO = d$ARO[match(X, d$query)], resfinder.norm.parent = d$parent[match(X, d$query)])



#saveRDS(args_abundances, file = "abundances_per_tool_and_unigene.rds", compress = T)

abundance_parent <- bind_rows(args_abundances %>% group_by(sample, rgi.diamond.parent) %>% 
  summarise(scaled = sum(scaled), raw = sum(raw), raw_unique = sum(raw_unique), normed10m = sum(normed10m)) %>% 
  mutate(tool = "rgi.diamond", parent = rgi.diamond.parent) %>% ungroup() %>% select(sample, parent, tool, scaled, raw, raw_unique, normed10m),
  
  args_abundances %>% group_by(sample, rgi.diamond.prot.parent) %>% 
    summarise(scaled = sum(scaled), raw = sum(raw), raw_unique = sum(raw_unique), normed10m = sum(normed10m)) %>% 
    mutate(tool = "rgi.diamond.prot", parent = rgi.diamond.prot.parent) %>% ungroup() %>% select(sample, parent, tool, scaled, raw, raw_unique, normed10m),
  
  args_abundances %>% group_by(sample, rgi.blast.parent) %>% 
    summarise(scaled = sum(scaled), raw = sum(raw), raw_unique = sum(raw_unique), normed10m = sum(normed10m)) %>% 
    mutate(tool = "rgi.blast", parent = rgi.blast.parent) %>% ungroup() %>% select(sample, parent, tool, scaled, raw, raw_unique, normed10m),
  
  args_abundances %>% group_by(sample, deeparg.norm.parent) %>% 
    summarise(scaled = sum(scaled), raw = sum(raw), raw_unique = sum(raw_unique), normed10m = sum(normed10m)) %>% 
    mutate(tool = "deeparg", parent = deeparg.norm.parent) %>% ungroup() %>% select(sample, parent, tool, scaled, raw, raw_unique, normed10m),
  
  args_abundances %>% group_by(sample, deeparg.norm.prot.parent) %>% 
    summarise(scaled = sum(scaled), raw = sum(raw), raw_unique = sum(raw_unique), normed10m = sum(normed10m)) %>% 
    mutate(tool = "deeparg.prot", parent = deeparg.norm.prot.parent) %>% ungroup() %>% select(sample, parent, tool, scaled, raw, raw_unique, normed10m),
  
  args_abundances %>% group_by(sample, fargene.parent) %>% 
    summarise(scaled = sum(scaled), raw = sum(raw), raw_unique = sum(raw_unique), normed10m = sum(normed10m)) %>% 
    mutate(tool = "fargene", parent = fargene.parent) %>% ungroup() %>% select(sample, parent, tool, scaled, raw, raw_unique, normed10m),
  
  args_abundances %>% group_by(sample, fargene.prot.parent) %>% 
    summarise(scaled = sum(scaled), raw = sum(raw), raw_unique = sum(raw_unique), normed10m = sum(normed10m)) %>% 
    mutate(tool = "fargene.prot", parent = fargene.prot.parent) %>% ungroup() %>% select(sample, parent, tool, scaled, raw, raw_unique, normed10m),
  
  args_abundances %>% group_by(sample, amrfinder.norm.parent) %>% 
    summarise(scaled = sum(scaled), raw = sum(raw), raw_unique = sum(raw_unique), normed10m = sum(normed10m)) %>% 
    mutate(tool = "amrfinder", parent = amrfinder.norm.parent) %>% ungroup() %>% select(sample, parent, tool, scaled, raw, raw_unique, normed10m),
  
  args_abundances %>% group_by(sample, amrfinder.norm.prot.parent) %>% 
    summarise(scaled = sum(scaled), raw = sum(raw), raw_unique = sum(raw_unique), normed10m = sum(normed10m)) %>% 
    mutate(tool = "amrfinder.prot", parent = amrfinder.norm.prot.parent) %>% ungroup() %>% select(sample, parent, tool, scaled, raw, raw_unique, normed10m),
  
  args_abundances %>% group_by(sample, abricate.argannot.norm.parent) %>% 
    summarise(scaled = sum(scaled), raw = sum(raw), raw_unique = sum(raw_unique), normed10m = sum(normed10m)) %>% 
    mutate(tool = "abricate.argannot", parent = abricate.argannot.norm.parent) %>% ungroup() %>% select(sample, parent, tool, scaled, raw, raw_unique, normed10m),
  
  args_abundances %>% group_by(sample, abricate.card.norm.parent) %>% 
    summarise(scaled = sum(scaled), raw = sum(raw), raw_unique = sum(raw_unique), normed10m = sum(normed10m)) %>% 
    mutate(tool = "abricate.card", parent = abricate.card.norm.parent) %>% ungroup() %>% select(sample, parent, tool, scaled, raw, raw_unique, normed10m),
  
  args_abundances %>% group_by(sample, abricate.megares.norm.parent) %>% 
    summarise(scaled = sum(scaled), raw = sum(raw), raw_unique = sum(raw_unique), normed10m = sum(normed10m)) %>% 
    mutate(tool = "abricate.megares", parent = abricate.megares.norm.parent) %>% ungroup() %>% select(sample, parent, tool, scaled, raw, raw_unique, normed10m),
  
  args_abundances %>% group_by(sample, abricate.ncbi.norm.parent) %>% 
    summarise(scaled = sum(scaled), raw = sum(raw), raw_unique = sum(raw_unique), normed10m = sum(normed10m)) %>% 
    mutate(tool = "abricate.ncbi", parent = abricate.ncbi.norm.parent) %>% ungroup() %>% select(sample, parent, tool, scaled, raw, raw_unique, normed10m),
  
  args_abundances %>% group_by(sample, abricate.resfinder.norm.parent) %>% 
    summarise(scaled = sum(scaled), raw = sum(raw), raw_unique = sum(raw_unique), normed10m = sum(normed10m)) %>% 
    mutate(tool = "abricate.resfinder", parent = abricate.resfinder.norm.parent) %>% ungroup() %>% select(sample, parent, tool, scaled, raw, raw_unique, normed10m),
  
  args_abundances %>% group_by(sample, resfinder.norm.parent) %>% 
    summarise(scaled = sum(scaled), raw = sum(raw), raw_unique = sum(raw_unique), normed10m = sum(normed10m)) %>% 
    mutate(tool = "resfinder", parent = resfinder.norm.parent) %>% ungroup() %>% select(sample, parent, tool, scaled, raw, raw_unique, normed10m))

metadata <- read.delim("metadata_GMGC10.sample.meta.tsv")
abundance_parent <- abundance_parent %>% filter(!is.na(parent)) %>% 
  mutate(habitat = metadata$habitat[match(sample, metadata$sample_id)])

saveRDS(abundance_parent, file = "abundance_per_tool.rds", compress = T)
write.csv(abundance_parent, file = "abundance_per_tool.csv", row.names = F)




