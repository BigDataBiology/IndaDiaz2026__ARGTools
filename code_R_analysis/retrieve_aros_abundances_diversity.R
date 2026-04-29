library(dplyr)
library(tidyr)
library(stringr)
library(rdflib)
library(stringr)
library(jsonlite)

################################################################################################################################################
# working directory
setwd("~/Documents/GitHub/arg_compare/")


# OLS
# lower ontology is the preferred level, for example 
# anything within (with lower ontology) antibiotic target modifying enzyme 
# ARO:3000519 will be assigned ARO:3000519, except for 
# Erm 23S ribosomal RNA methyltransferase ARO:3000560

higher_ontology <- paste0("ARO:", c("3000557", "0000010", "3000519", "3000185", "3000381", "3000159", "3000012"))
lowest_ontology <- c(paste0("ARO:", c("3000341", "3000322", "3000345", "3004257", "3007419", "3004276", "3004275", "3000229", "3000225", "3000228",
                                      "3000128", "3000127", "3000126", "3000155", "3000151", "3000154", "3000153", "3004260", "3000078", "3000004", 
                                      "3000076", "3000075", "3007074", "3000122", "3000249", "3004467", "3004064", "3000342", "3003025", "3000221", 
                                      "3000320", "3000458", "3000333", "3007103", "3000576", "3000233", "3000869", "3000036", "3004261", "3000231", 
                                      "3000234", "3007428", "0000031", "3000560", "3004469", "3000419", "3000507", "3005086", "0000002", "3003425", 
                                      "3001208", "3000210", "3004238", "3003040", "0010002", "3003580", "3003768", "3001207", "3000492", "3007610", 
                                      "3000100", "3007429", "3000270", "3000451", "3002976", "3007425", "3004916")))

ontologies <- c(lowest_ontology, higher_ontology)
# This are the original levels 
# description of each ARO
# antibiotic inactivation enzyme ARO:3000557 NOT IN 
#   # AAC(2') ARO:3000341 (alternatively only ARO:3000121)
#   # AAC(3) ARO:3000322 (alternatively only ARO:3000121)
#   # AAC(6') ARO:3000345 (alternatively only ARO:3000121)
#   # cpa ARO:3004257
#   # aminoglycoside bifunctional ARO:3007419
#   # ant2 ARO:3004276 (alternatively only ARO:3000218)
#   # ant3 ARO:3004275 (alternatively only ARO:3000218)
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


## FARGENE CONVERSIONS
## simplified names for fargene models 
fg_class <- c("aac2p","aac6p","aac6p","aac3", "aph3p","aph6p","class_a","aac3", "mph",
             "class_b1_b2", "class_b3","class_d","aac6p","erm","class_d","class_c",
             "aph2b","erm","tet_enzyme", "tet_rpg", "tet_efflux", "qnr")

# Output names of the models, since there are different HMM's per class, we reduce them to just the class 

names(fg_class) <- c("aac2p","aac6p_1","aac6p_2","aac3_2", "aph3p","aph6p","class_a","aac3_1", "mph",
                    "class_b1_b2", "class_b3","class_d1","aac6p_3","erm_2","class_d2","class_c",
                    "aph2b","erm_1","tet_enzyme", "tet_rpg", "tet_efflux", "qnr")

# This applies for both fargene_results and fargene_hmm files. To fetch the hmm score we also need the model name from the hmm file below

hmm_models <- c("aac2p", "aac3_1", "aac3_2", "aac6p_1", "aac6p_2", "aac6p_3", "aph2b", "aph3p", "aph6p", "class_b1_b2", "class_b3",
               "class_c", "class_a", "tet_efflux", "tet_enzyme", 
               "mph", "erm_1", "erm_2", "class_d1","class_d2", "tet_rpg", "qnr")

names(hmm_models) <- c("aac2p-aligned", "aac3_class1-aligned", "aac3_class2-aligned", "aac6p_class1-aligned", "aac6p_class2-aligned",
"aac6p_class3-aligned", "aph2b-aligned", "aph3p-aligned", "aph6-aligned", "b1_b2_70_centroids-aligned", "b3_70_centroids-aligned", 
"class_C_70_centroids-aligned", "classA_70_centroids-aligned", "efflux_model_group_1-aligned", "enzyme_reduced_tetX1_X3-aligned",
"macrolide_phosphotransferases-aligned", "methyltransferase_grp1-aligned", "methyltransferase_grp2-aligned", 
"oxa_g1_70_centroids-aligned", "oxa_g2_70_centroids-aligned", "rpg_reference_sequences-aligned", "pmqnr_20120719.pfa")

# Manually assigned aros to fargene classes
fargene2ARO <- c("ARO:3000341", "ARO:3000322", "ARO:3000345", "ARO:3000128", "ARO:3000126", "ARO:3000151", "ARO:3000151",
                 "ARO:3000078", "ARO:3000004", "ARO:3000004", "ARO:3000076", "ARO:3000075",
                 "ARO:3000560", "ARO:3000333", "ARO:3000036", "ARO:0000002", "ARO:0010002", "ARO:3000419")

names(fargene2ARO) <- c("aac2p", "aac3", "aac6p", "aph2b", "aph3p","aph6", "aph6p",
                        "class_a", "class_b1_b2", "class_b3", "class_c", "class_d", 
                        "erm", "mph", "tet_enzyme", "tet_rpg", "tet_efflux", "qnr")


################################################################################################################################################
################################################################################################################################################
## FETCH THE ONTOLOGY from CARD

## We use this to simplify the ontology and summarise the abundance per class and not per unigene
## For each identified ARO by any of the tools, we will assign it's lowest_ontology, 
# and if the ARO is higher than the lowest_ontology, we will assign higher_ontology 

# aro_url <- "https://raw.githubusercontent.com/arpcard/aro/master/aro.owl"
# g <- rdf_parse(aro_url, format = "rdfxml")
g <- rdf_parse("data/aro.owl", format = "rdfxml")
# necessary namespaces
rdfs_ns <- "http://www.w3.org/2000/01/rdf-schema#"

# Parent-Child Relationships**
query_parents <- paste0("PREFIX rdfs: <", rdfs_ns, "> ",
                        "SELECT ?child ?parent WHERE { ?child rdfs:subClassOf ?parent }")
parent_map <- rdf_query(g, query_parents)  # DataFrame with columns "child" and "parent"

# Convert URIs to ARO:XXXX format
parent_map <- parent_map %>%
  mutate(Child_ID = str_extract(child, "ARO_[0-9]+") %>% str_replace_all("_", ":"),
         Parent_ID = str_extract(parent, "ARO_[0-9]+") %>% str_replace_all("_", ":")) %>%
  select(Child_ID, Parent_ID) %>%
  filter(!is.na(Child_ID) & !is.na(Parent_ID))  # Remove NA values

# Preload Labels**
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


# LOAD THE RESULTS FROM EACH TOOL

# We standardize the unigene column to query
# if the tool reports gene class and subclass, we call rename them
# we assigned tool name as well

# amrfinder 

amrfinder.norm <- read.delim("dna/amrfinder.norm.tsv",skip = 1) %>% 
  rename(query = Contig.id) %>% 
  rename(ARG.class = Class) %>% 
  rename(ARG.subclass = Subclass) %>%
  mutate(tool = "AMRFinderPlus-nt", id = X..Identity.to.reference) %>% 
  group_by(query) %>% 
  arrange(desc(X..Identity.to.reference + X..Coverage.of.reference), 
          desc(Alignment.length)) %>%
  slice_head(n=1) %>% ungroup() # leave one match in repeated queries

amrfinder.norm.prot <- read.delim("protein/amrfinder.norm.tsv", skip = 1) %>% 
  rename(query = Protein.id) %>% 
  rename(ARG.class = Class) %>% 
  rename(ARG.subclass = Subclass) %>%
  mutate(tool = "AMRFinderPlus", id = X..Identity.to.reference) %>% 
  group_by(query) %>% 
  arrange(desc(X..Identity.to.reference + X..Coverage.of.reference), 
          desc(Alignment.length)) %>%
  slice_head(n=1) %>% ungroup() # leave one match in repeated queries

# deeparg
deeparg.norm <- read.delim("dna/deeparg.norm.tsv", skip = 1) %>% 
  rename(query = read_id) %>% 
  rename(ARG.class = predicted_ARG.class) %>%
  mutate(tool = "DeepARG", id = identity)

# no repeated unigenes 

deeparg.norm.prot <- read.delim("protein/deeparg.norm.tsv", skip = 1) %>% 
  rename(query = read_id) %>% 
  rename(ARG.class = predicted_ARG.class) %>%
  mutate(tool = "DeepARG-aa", id = identity)

# no repeated unigenes 

# rgi
rgi.diamond <- read.delim("dna/rgi_diamond.tsv") %>% 
  select(-c(ORF_ID, Predicted_DNA, Predicted_Protein, CARD_Protein_Sequence)) %>% 
  mutate(query = gsub('.{2}$', '', Contig)) %>% 
  rename(ARG.class = AMR.Gene.Family) %>%
  mutate(ARO = paste0("ARO:", ARO)) %>%
  mutate(tool = "RGI-DIAMOND", id = Best_Identities) %>% 
  group_by(query) %>% 
  arrange(desc(Best_Identities), 
          desc(Best_Hit_Bitscore), 
          desc(Percentage.Length.of.Reference.Sequence)) %>%
  slice_head(n=1) %>% ungroup()

# 4 repeated unigenes, 6 copies deleted 

rgi.blast <- read.delim("dna/rgi_blast.tsv") %>% 
  select(-c(ORF_ID, Predicted_DNA, Predicted_Protein, CARD_Protein_Sequence)) %>% 
  mutate(query = gsub('.{2}$', '', Contig)) %>% 
  rename(ARG.class = AMR.Gene.Family) %>%
  mutate(ARO = paste0("ARO:", ARO)) %>%
  mutate(tool = "RGI-BLAST", id = Best_Identities) %>% 
  group_by(query) %>% 
  arrange(desc(Best_Identities), 
          desc(Best_Hit_Bitscore), 
          desc(Percentage.Length.of.Reference.Sequence)) %>%
  slice_head(n=1) %>% ungroup()


# 4 repeated unigenes, 7 copies deleted 

rgi.diamond.prot <- read.delim("protein/rgi_diamond.tsv") %>% 
  select(-c(Predicted_DNA, Predicted_Protein, CARD_Protein_Sequence)) %>% 
  mutate(query = ORF_ID) %>% 
  rename(ARG.class = AMR.Gene.Family) %>% 
  mutate(ARO = paste0("ARO:", ARO)) %>%
  mutate(tool = "RGI-DIAMOND-aa", id = Best_Identities)

# no repeated unigenes


# fargene
fargene <- read.delim("dna/fargene_results.tsv", header = F) %>%
  mutate(query = gsub('.{2}$', '', V1)) %>% 
  mutate(new_class =  fg_class[V5]) %>% 
  group_by(query, new_class) %>% 
  slice_head(n=1) %>% ungroup()

# 48 queries with repeated class removed

## remove duplicated queries with different classes 
# load the hmm scores
hmm <- read.table("dna/fargene_hmm.txt", quote="\"", comment.char="") %>%
  mutate(query = gsub('.{2}$', '', V1),
         q1 = sapply(strsplit(V1, split = "GMGC10"), function(x) paste0("GMGC10",x[length(x)])),
         new_class = hmm_models[V4]) %>% 
  mutate(q1 = gsub('.{2}$', '', q1)) %>% 
  mutate(q1 = sapply(strsplit(q1, split = "_seq"), function(x) x[1])) %>% 
  arrange(query, desc(V14)) 

fargene <- fargene %>% 
  mutate(hmm = hmm$V14[match(paste(query, V5), paste(hmm$q1, hmm$new_class))],
         tool = "fARGene") 

fargene <- fargene %>% ungroup() %>% 
  group_by(query) %>% 
  arrange(desc(hmm)) %>% 
  slice_head(n=1) %>% 
  ungroup()
# 19 observations with different class were deleted, 1 observation left per unigene 

rm(hmm)

## identity levels with RGI (loose, strict, perfect) for those genes identified by fargene

fargene_with_rgi <- read.delim("check_missing_annot/fargene_predicted_fna.txt") %>% 
  select(-c(ORF_ID, Predicted_DNA, Predicted_Protein, CARD_Protein_Sequence)) %>% 
  mutate(query = gsub(' ', '', Contig)) %>% 
  mutate(query = gsub('.{2}$', '', query)) %>% 
  mutate(query = gsub('tmp_', '', query)) %>% 
  mutate(ARO = paste0("ARO:", ARO))


fargene <- fargene %>% 
  mutate(id = fargene_with_rgi$Best_Identities[match(query, fargene_with_rgi$query)],
         aro.rgi = fargene_with_rgi$ARO[match(query, fargene_with_rgi$query)],
         coverage.rgi = fargene_with_rgi$Percentage.Length.of.Reference.Sequence[match(query, fargene_with_rgi$query)],)

rm(fargene_with_rgi)


fargene <- fargene %>% mutate(manual.ARO = as.vector(fargene2ARO[new_class]),
                              ARO = ifelse(is.na(aro.rgi), "", aro.rgi))

## this part is after assigning new_levels and giving the manual new_level to the observations missing it (it is in retrospective)
# for fargene we took ARO:3002484 and made them all class D
# "ARO:3002600"  AAC(3)-Ib/AAC(6')-Ib3, this is aac(6) no matter what, the coverage is below 45%, id 2 below 70% 
# "ARO:3002599" this is aac(6) no matter what, the coverage is below 40%, id 2 are 100, 6 more below 70% 
# "ARO:3002597"AAC(6')-Ie-APH(2'')-Ia, most coverage are below 40%, 9 below 60%, we give the fargene class instead
# The  MPH that did not agree with fargene have id levels below 30%, they are aph, so very similar, we give the fargene class instead
# The beta-lactam modulation all have id levels below 50%, we give the fargene class instead
# The V/M all have id levels below 40%, we give the fargene class instead


fargene.prot <- read.delim("protein/fargene_results.tsv", header = F) %>% 
  mutate(query = V1, 
         new_class = fg_class[V5]) %>%
  group_by(query, new_class) %>% 
  slice_head(n = 1) %>% ungroup()

# 41 unigenes with repeated class, 41 observations deleted, 1 left per unigene 

# load the HMM scgroup_by()# load the HMM scores 
hmm.prot <- read.table("protein/fargene_hmm.txt", quote="\"", comment.char="") %>%
  mutate(new_class = hmm_models[V4], q1 = V1) %>% 
  arrange(q1, desc(V14))

fargene.prot <- fargene.prot %>% 
  mutate(hmm = hmm.prot$V14[match(paste(V1, V5), paste(hmm.prot$q1, hmm.prot$new_class))])

fargene.prot <- fargene.prot %>% ungroup() %>%
  group_by(query) %>% arrange(desc(hmm)) %>% slice_head(n=1) %>% ungroup()

rm(hmm.prot)

# 19 unigenes with different class, 19 observations removed, 1 left per unigene   
# fargene.prot %>% group_by(query) %>% mutate(n = n()) %>% filter(n >1) %>% select(query) %>% distinct()

# rename aph6p in fargene and fargene.prot

fargene.prot <- fargene.prot %>% 
  mutate(new_class = ifelse(new_class %in% "aph6p", "aph6", new_class))
fargene <- fargene %>% 
  mutate(new_class = ifelse(new_class %in% "aph6p", "aph6", new_class))

fargene.prot <- fargene.prot %>%
  mutate(tool = "fARGene-aa")

# identity levels with RGI 

fargene_with_rgi <- read.delim("check_missing_annot/fargene_prot_predicted_faa.txt") %>% 
  mutate(query = gsub(' ', '', ORF_ID)) %>% 
  mutate(query = gsub('tmp_', '', query),
         ARO = paste0("ARO:", ARO))

fargene.prot <- fargene.prot %>% 
  mutate(id = fargene_with_rgi$Best_Identities[match(query, fargene_with_rgi$query)],
         aro.rgi = fargene_with_rgi$ARO[match(query, fargene_with_rgi$query)],
         coverage.rgi = fargene_with_rgi$Percentage.Length.of.Reference.Sequence[match(query, fargene_with_rgi$query)],)

rm(fargene_with_rgi)


fargene.prot <- fargene.prot %>% mutate(manual.ARO = as.vector(fargene2ARO[new_class]),
                                        ARO = ifelse(is.na(aro.rgi), "", aro.rgi))

# use AROs from fargene.prot in fargene and vice versa
fargene <- fargene %>% mutate(ARO = ifelse(ARO == "", fargene.prot$ARO[match(query, fargene.prot$query)], ARO)) %>%
  mutate(ARO = ifelse(is.na(ARO), "", ARO))
fargene.prot <- fargene.prot %>% mutate(ARO = ifelse(ARO == "", fargene$ARO[match(query, fargene$query)], ARO)) %>%
  mutate(ARO = ifelse(is.na(ARO), "", ARO))

# after tblastx in the CARD website 
# "GMGC10.029_357_382.STRB_1", "GMGC10.008_012_191.SPH_1", GMGC10.170_456_610.STRB_1 are given  "ARO:3002658"

fargene <- fargene %>% 
  mutate(ARO = ifelse(query %in% c("GMGC10.029_357_382.STRB", "GMGC10.008_012_191.SPH", "GMGC10.170_456_610.STRB"), "ARO:3002658", ARO))
fargene.prot <- fargene.prot %>% 
  mutate(ARO = ifelse(query %in% c("GMGC10.029_357_382.STRB", "GMGC10.008_012_191.SPH", "GMGC10.170_456_610.STRB"), "ARO:3002658", ARO))


# manual blast search in CARD website
manual_blast_fg <- cbind(c("GMGC10.002_063_752.UNKNOWN", "GMGC10.014_038_935.UNKNOWN", "GMGC10.018_009_687.UNKNOWN", "GMGC10.039_734_695.UNKNOWN", "GMGC10.155_234_566.UNKNOWN",
  "GMGC10.155_509_478.UNKNOWN","GMGC10.165_346_918.UNKNOWN", "GMGC10.172_669_160.UNKNOWN", "GMGC10.182_555_377.UNKNOWN", "GMGC10.207_886_721.UNKNOWN",
  "GMGC10.210_614_745.UNKNOWN", "GMGC10.213_353_700.UNKNOWN", "GMGC10.218_878_822.UNKNOWN", "GMGC10.237_491_692.UNKNOWN", "GMGC10.279_981_496.UNKNOWN",
  "GMGC10.282_839_022.UNKNOWN", "GMGC10.287_300_498.UNKNOWN"),
  c("ARO:3002528", "ARO:3009063", "ARO:3003676", "ARO:3003199", "ARO:3002589",
  "ARO:3002644","ARO:3002645","ARO:3002524", "ARO:3002999", "ARO:3002572",
  "ARO:3004359", "ARO:3002571", "ARO:3003720", "ARO:3006942", "ARO:3004621",
  "ARO:3002549", "ARO:3002644"))

fargene <- fargene %>% 
  mutate(ARO = ifelse(query %in% manual_blast_fg[,1], manual_blast_fg[match(query, manual_blast_fg[ , 1]), 2], ARO))

fargene.prot <- fargene.prot %>% 
  mutate(ARO = ifelse(query %in% manual_blast_fg[,1], manual_blast_fg[match(query, manual_blast_fg[ , 1]), 2], ARO))

###
# abricate
abricate.argannot.norm <- read.delim("dna/abricate-argannot.norm.tsv", header=FALSE, comment.char="#") %>%
  select(-V1) %>% 
  rename(query = V2, ARO = V16) %>% 
  mutate(ARG.class = str_match(V14, "\\(([^)]+)\\)")[,2],
         gene = sub("\\([^)]*\\)", "", V14),
         tool = "ABRicate-ARGANNOT", 
         id = V11) %>% 
  group_by(query) %>% 
  arrange(desc(id) + desc(V10)) %>% 
  slice_head(n = 1) %>% 
  ungroup()

# 6 unigenes with repeated observations, 1 with same gene class, 6 observations removed 
# abricate.argannot.norm %>% group_by(query) %>% mutate(n = n()) %>% filter(n>1) %>% select(query) %>% distinct()


abricate.card.norm <- read.delim("dna/abricate-card.tsv", header=FALSE, comment.char="#") %>%
  select(-V1) %>% 
  rename(query = V2, ARG.class = V15) %>%
  mutate(tool = "ABRicate-CARD", id = V11) %>% 
  group_by(query) %>% 
  arrange(desc(id) + desc(V10)) %>% 
  slice_head(n = 1) %>% 
  ungroup()

# 18 unigenes with repeated observations, 2 with same gene class, 6 observations removed , 1 left per unigene
# abricate.card.norm  %>% group_by(query) %>% mutate(n = n()) %>% filter(n>1) %>% select(query) %>% distinct()
# abricate.card.norm  %>% group_by(query, ARG.class) %>% mutate(n = n()) %>% filter(n>1) 
# ARO COMES LATER 


abricate.megares.norm <- read.delim("dna/abricate-megares.norm.tsv", header=FALSE, comment.char="#") %>%
  select(-V1) %>% 
  rename(query = V2, ARO = V16) %>%
  mutate(tool = "ABRicate-MEGARes", id = V11) %>% 
  group_by(query) %>% 
  arrange(desc(id) + desc(V10)) %>% 
  slice_head(n = 1) %>% 
  ungroup()

# 20 unigenes with repeated observations, 20 observations removed , 1 left per unigene
# abricate.megares.norm  %>% group_by(query) %>% mutate(n = n()) %>% filter(n>1) %>% select(query) %>% distinct()

abricate.ncbi.norm <- read.delim("dna/abricate-ncbi.norm.tsv", header=FALSE, comment.char="#") %>%
  select(-V1) %>% 
  rename(query = V2) %>% rename(ARG.class = V15, gene = V14, ARO = V16) %>%
  mutate(tool = "ABRicate-NCBI", id = V11) %>% 
  group_by(query) %>% 
  arrange(desc(id) + desc(V10)) %>% 
  slice_head(n = 1) %>% 
  ungroup()

# 10 unigenes with repeated observations (1 same gene), 10 observations removed , 1 left per unigene
# abricate.ncbi.norm  %>% group_by(query) %>% mutate(n = n()) %>% filter(n>1) 
# abricate.ncbi.norm  %>% group_by(query) %>% mutate(n = n()) %>% filter(n>1) %>% select(query) %>% distinct()

abricate.resfinder.norm <- read.delim("dna/abricate-resfinder.norm.tsv", header=FALSE, comment.char="#") %>%
  select(-V1) %>% 
  rename(query = V2, drug = V15, gene = V14, ARO = V16) %>%
  mutate(tool = "ABRicate-ResFinder", id = V11) %>% 
  group_by(query) %>% 
  arrange(desc(id) + desc(V10)) %>% 
  slice_head(n = 1) %>% 
  ungroup()

# 7 unigenes with repeated observations, 7 observations removed , 1 left per unigene
# abricate.resfinder.norm  %>% group_by(query) %>% mutate(n = n()) %>% filter(n>1) 
# abricate.resfinder.norm  %>% group_by(query) %>% mutate(n = n()) %>% filter(n>1) %>% select(query) %>% distinct()


# resfinder
resfinder.norm <- read.delim("dna/resfinder_table.norm.tsv", skip = 1) %>% 
  rename(query = Contig) %>% rename(gene = Resistance.gene) %>%
  mutate(tool = "ResFinder", id = Identity)  


# resfinder.norm  %>% group_by(query) %>% mutate(n = n()) %>% filter(n>1) 
# resfinder.norm  %>% group_by(query) %>% mutate(n = n()) %>% filter(n>1) %>% select(query) %>% distinct()
# resfinder.norm  %>% group_by(query, gene) %>% mutate(n = n()) %>% filter(n>1)
# resfinder.norm  %>% group_by(query, gene) %>% mutate(n = n()) %>% filter(n>1) %>% select(query) %>% distinct()

# 220 repeated unigenes
# 109 have the same gene 
## 113 have different gene 

resfinder.norm <- resfinder.norm %>% 
  group_by(query, gene) %>% 
  arrange(desc(id) + desc(Coverage)) %>% 
  slice_head(n = 1) %>% 
  ungroup() %>% 
  group_by(query) %>% 
  arrange(desc(id) + desc(Coverage)) %>% 
  slice_head(n = 1) %>% 
  ungroup()


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

resfinder.norm <- resfinder.norm %>% 
  mutate(ARG.class = phenotypes_resfiner$Class[v])

rm(phenotypes_resfiner, v)


################################################################################################################################################
################################################################################################################################################

# Checking missing AROs per dataset 

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

card_gene_aro <- card_gene_aro %>% mutate(aro = paste0("ARO:",aro))

# load the manually curated AROS
manually_curated <- read.table("check_missing_annot/abricate-card-manual-curation.txt", stringsAsFactors = F, text = T)
manually_curated <- data.frame(gene = sapply(strsplit(manually_curated$V1, split = ","), function(x) x[1]), aro = sapply(strsplit(manually_curated$V1, split = ","), function(x) x[2]))
manually_curated <- manually_curated %>% mutate(aro = paste0("ARO:",aro))

# add them to the data retrieved from json file

card_gene_aro <- card_gene_aro %>% bind_rows(manually_curated)

# add the ARO to abricate.card.norm
abricate.card.norm <- abricate.card.norm %>% mutate(ARO = card_gene_aro$aro[match(V6, card_gene_aro$gene)])
# sum(abricate.card.norm$ARO == "")

################################################################################################################################################

# ARO for abricate.argannot.norm
# manually assign the following 
# abricate.argannot.norm[abricate.argannot.norm$ARO == "",]

abricate.argannot.norm <- abricate.argannot.norm %>% 
  mutate(ARO = ifelse(V6 %in% "(Ntmdz)nimj_Nitroimidazole_Gene", "ARO:3007112", 
                      ifelse(V6 %in% "(Bla)blaACT-16", "ARO:3001827", ARO)))

# sum(abricate.argannot.norm$ARO == "")

################################################################################################################################################
# ARO for abricate.resfinder.norm
resfinder_ARO_mapping <- read.delim("check_missing_annot/resfinder_ARO_mapping.tsv")
#resfinder_ARO_mapping$Gene.Name.in.CARD[1000:1100]

# correct the capital letters in these two entries

abricate.resfinder.norm <- abricate.resfinder.norm %>% 
  mutate(V6 = ifelse(V6 %in% "oqxA_1", "OqxA_1", 
                     ifelse(V6 %in% "oqxB_1", "OqxB_1", V6))) 

j <- which(abricate.resfinder.norm$ARO == "")

resfinder_missing_aro <- tibble(gene=abricate.resfinder.norm$V6[j], accession=abricate.resfinder.norm$V13[j]) %>% distinct() %>%
  filter(!is.na(match(paste(gene, accession, sep = "_"), resfinder_ARO_mapping$Original.ID)))

k <- match(paste(abricate.resfinder.norm$V6[j], abricate.resfinder.norm$V13[j], sep = "_"), resfinder_ARO_mapping$Original.ID)

abricate.resfinder.norm$ARO[j]  <- paste0("ARO:", resfinder_ARO_mapping$ARO[k])

abricate.resfinder.norm <- abricate.resfinder.norm %>% 
  mutate(ARO = ifelse(ARO == "ARO:NA", "", ARO))


others <- cbind(c("mdf(A)_1", "mcr-9_1", "blaCTX-M-63_1", "blaCARB-4_1", "dfrA19_1"), 
                  c("ARO:3001328", "ARO:3004684", "ARO:3001924", "ARO:3002243", "ARO:3003015"))

abricate.resfinder.norm <- abricate.resfinder.norm %>% mutate(ARO = ifelse(V6 %in% others[,1], others[match(V6, others[,1]),2], ARO))
# 30 unigenes left - 5 genes "mdf(A)"      "mcr-9"       "blaCTX-M-63" "blaCARB-4"   "dfrA19"
#sum(abricate.resfinder.norm$ARO == "")


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

colistin <- cbind(c("mcr-3.38", "mcr-3.36", "mcr-8.2", "mcr-3.33", "mcr-10.2"), 
      c("ARO:3007251", "ARO:3007257", "ARO:3007229", "ARO:3007248", "ARO:3007277"))

resfinder.norm <- resfinder.norm %>% mutate(ARO = ifelse(gene %in% colistin[,1], colistin[match(gene, colistin[,1]),2], ARO))
#sum(resfinder.norm$ARO == "")


################################################################################################################################################
# ARO for abricate megares

megares_annotation <- read.delim("check_missing_annot/megares_ARO_mapping.tsv", header = T)
megares_annotation <- megares_annotation %>% mutate(ref = sapply(strsplit(megares_annotation$Original.ID, split = "\\|"), function(x) x[1]))

j <- which(abricate.megares.norm$ARO == "")
d <- paste0("ARO:", megares_annotation$ARO[match(abricate.megares.norm$V13[j], megares_annotation$ref)])
d[d=="ARO:NA"] <- ""
abricate.megares.norm$ARO[j] <- d

j <- which(abricate.megares.norm$ARO == "")
#unique(abricate.megares.norm[j,"V6"])
#unique(abricate.megares.norm[j,"V13"])

abricate_annotation <- read.delim("check_missing_annot/abricate_annotation_megares.tsv")
abricate_annotation <- abricate_annotation[,-c(19,20,21)]
abricate_annotation$MEG <- sapply(strsplit(abricate_annotation$ORF, split = "~~~"), function(x) x[3])
abricate_annotation$gene <- sapply(strsplit(abricate_annotation$ORF, split = "~~~"), function(x) x[2])
j <- which(abricate.megares.norm$ARO == "")
abricate.megares.norm$ARO[j] <- paste0("ARO:", abricate_annotation$ARO[match(abricate.megares.norm$V13[j], abricate_annotation$MEG)])
abricate.megares.norm$ARO[abricate.megares.norm$ARO == "ARO:NA"] <- ""

j <- which(abricate.megares.norm$ARO == "")
#unique(abricate.megares.norm[j,"V6"])
#unique(abricate.megares.norm[j,"V13"])
#unique(abricate.megares.norm[j,"V14"])$V14

abricate.megares.norm <- abricate.megares.norm %>% mutate(code_gene = sapply(strsplit(V14, split = ":"), function(x) x[length(x)]))
missing <- unique(abricate.megares.norm$V14[j])
new_code <- rep(NA, length(missing))

# RND EFFLUX PUMPS
new_code[grepl("RND", missing)] <- "ARO:0010004"
# MATE EFFLUX PUMPS
new_code[grepl("MATE", missing)] <- "ARO:3000112"
# MFS EFFLUX PUMPS
new_code[grepl("MFS", missing)] <- "ARO:0010002"
# SMR EFFLUX PUMPS
new_code[grepl("SMR", missing)] <- "ARO:0010003"

abricate.megares.norm$ARO[j] <- new_code[match(abricate.megares.norm$V14[j], missing)]


# to debug arg_norm

rm(j, abricate_annotation, megares_annotation, new_code, missing)

################################################################################################################################################
# ARO for abricate NCBI

# from argnorm

ncbi_annotation <- read.delim("check_missing_annot/ncbi_ARO_mapping.tsv")
ncbi_annotation <- ncbi_annotation %>% mutate(ref = sapply(strsplit(ncbi_annotation$Original.ID, split = "\\|"), function(x) x[2]))
ncbi_annotation <- ncbi_annotation %>% mutate(genename = sapply(strsplit(ncbi_annotation$Original.ID, split = "\\|"), function(x) x[5]))
ncbi_annotation <- ncbi_annotation %>% mutate(genename2 = sapply(strsplit(ncbi_annotation$Original.ID, split = "\\|"), function(x) x[6]))

# using gene name 
j <- abricate.ncbi.norm$ARO ==""
d <- paste0("ARO:", ncbi_annotation$ARO[match(abricate.ncbi.norm$V6[j], ncbi_annotation$genename)])
d[d=="ARO:NA"] <- ""
abricate.ncbi.norm$ARO[j] <- d

#  
j <- abricate.ncbi.norm$ARO ==""
d <- paste0("ARO:", ncbi_annotation$ARO[match(abricate.ncbi.norm$V6[j], ncbi_annotation$genename2)])
d[d=="ARO:NA"] <- ""
abricate.ncbi.norm$ARO[j] <- d

# using reference number
j <- abricate.ncbi.norm$ARO ==""
d <- paste0("ARO:", ncbi_annotation$ARO[match(abricate.ncbi.norm$V13[j], ncbi_annotation$ref)])
d[d=="ARO:NA"] <- ""
abricate.ncbi.norm$ARO[j] <- d


ncbi_annotation <- read.delim("check_missing_annot/ncbi_curation.tsv")
ncbi_annotation <- ncbi_annotation %>% mutate(ref = sapply(strsplit(ncbi_annotation$Original.ID, split = "\\|"), function(x) x[2]))
ncbi_annotation <- ncbi_annotation %>% mutate(genename = sapply(strsplit(ncbi_annotation$Original.ID, split = "\\|"), function(x) x[5]))
ncbi_annotation <- ncbi_annotation %>% mutate(genename2 = sapply(strsplit(ncbi_annotation$Original.ID, split = "\\|"), function(x) x[6]))

j <- abricate.ncbi.norm$ARO ==""
d <- paste0("ARO:", ncbi_annotation$ARO[match(abricate.ncbi.norm$V6[j], ncbi_annotation$genename)])
d[d=="ARO:NA"] <- ""
abricate.ncbi.norm$ARO[j] <- d

j <- abricate.ncbi.norm$ARO ==""
d <- paste0("ARO:", ncbi_annotation$ARO[match(abricate.ncbi.norm$V6[j], ncbi_annotation$genename2)])
d[d=="ARO:NA"] <- ""
abricate.ncbi.norm$ARO[j] <- d

j <- abricate.ncbi.norm$ARO ==""
d <- paste0("ARO:", ncbi_annotation$ARO[match(abricate.ncbi.norm$V13[j], ncbi_annotation$ref)])
d[d=="ARO:NA"] <- ""
abricate.ncbi.norm$ARO[j] <- d


# MANUALLY ASSIGNED THESE ONTOLOGIES
abricate_annotation <- read.delim("check_missing_annot/abricate_annotation.tsv")
abricate_annotation <- abricate_annotation[,-c(19,20,21)]
x2 <- sapply(strsplit(abricate_annotation$ORF, split = "~~~"), function(x) x[2])
x3 <- sapply(strsplit(abricate_annotation$ORF, split = "~~~"), function(x) x[3])

j <- abricate.ncbi.norm$ARO ==""
d <- paste0("ARO:", abricate_annotation$ARO[match(abricate.ncbi.norm$V6[j], x2)])
d[d=="ARO:NA"] <- ""
abricate.ncbi.norm$ARO[j] <- d

j <- abricate.ncbi.norm$ARO ==""
unique(abricate.ncbi.norm[j,]$V6)

# MANUALLY ASSIGNED THESE ONTOLOGIES
#"penI_Bp" from fargene 
# ARO:3002956 from amrfinder 
# the rest from card 

others.ncbi <- cbind(c("penI_Bp", "aadE-Cc", "vanY-32487235", "dfr7"), 
                c("ARO:3009015", "ARO:3002626", "ARO:3002956", "ARO:3002862"))

abricate.ncbi.norm <- abricate.ncbi.norm %>% mutate(ARO = ifelse(V6 %in% others.ncbi[,1], others.ncbi[match(V6, others.ncbi[,1]),2], ARO))
j <- abricate.ncbi.norm$ARO ==""
unique(abricate.ncbi.norm[j,]$V6)
rm(others.ncbi, others, abricate_annotation, x2, x3, j)

# to debug arg_norm
rm(d, ncbi_annotation)

### AMRFINDER
# 
# 

ncbi_annotation <- read.delim("check_missing_annot/ncbi_ARO_mapping.tsv")
ncbi_annotation <- ncbi_annotation %>% mutate(ref = sapply(strsplit(ncbi_annotation$Original.ID, split = "\\|"), function(x) x[2]))

j <- amrfinder.norm.prot$ARO ==""
d <- paste0("ARO:", ncbi_annotation$ARO[match(amrfinder.norm.prot$Closest.reference.accession[j], ncbi_annotation$ref)])
d[d=="ARO:NA"] <- ""
amrfinder.norm.prot$ARO[j] <- d

# to debug arg_norm

rm(j, d, ncbi_annotation)


################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
### AROs

# Create a list of unique AROs in all tools 

aros <- tibble(data.frame(rbind(cbind(rgi.diamond$tool, rgi.diamond$Best_Hit_ARO, rgi.diamond$ARO),
                                cbind(rgi.diamond.prot$tool, rgi.diamond.prot$Best_Hit_ARO, rgi.diamond.prot$ARO),
                                cbind(rgi.blast$tool, rgi.blast$Best_Hit_ARO, rgi.blast$ARO),
                                cbind(deeparg.norm$tool, deeparg.norm$X.ARG, deeparg.norm$ARO),
                                cbind(deeparg.norm.prot$tool, deeparg.norm.prot$X.ARG, deeparg.norm.prot$ARO),
                                cbind(amrfinder.norm$tool, amrfinder.norm$Element.symbol, amrfinder.norm$ARO),
                                cbind(amrfinder.norm.prot$tool, amrfinder.norm.prot$Element.symbol, amrfinder.norm.prot$ARO),
                                cbind(abricate.argannot.norm$tool, abricate.argannot.norm$V6, abricate.argannot.norm$ARO),
                                cbind(abricate.megares.norm$tool, abricate.megares.norm$V6, abricate.megares.norm$ARO),
                                cbind(abricate.ncbi.norm$tool, abricate.ncbi.norm$V6, abricate.ncbi.norm$ARO),
                                cbind(abricate.resfinder.norm$tool, abricate.resfinder.norm$V6, abricate.resfinder.norm$ARO),
                                cbind(abricate.card.norm$tool, abricate.card.norm$V6, abricate.card.norm$ARO),
                                cbind(resfinder.norm$tool, resfinder.norm$gene, resfinder.norm$ARO),
                                cbind(fargene$tool, fargene$new_class, fargene$ARO),
                                cbind(fargene.prot$tool, fargene.prot$new_class, fargene.prot$ARO))))

aros <- aros %>% distinct()
#aros %>% group_by(X1) %>% summarise(n = n(), e = sum(X3 == "")) %>% mutate(p=e/n) %>% arrange(desc(n))

# here we use the ontology information g, parent_map, label_map, rdfs_ns, aro_url, query_parents, parent_map

df <- process_terms_fast(unique(aros$X3[aros$X3!=""]), parent_map, label_map)
rm(g, rdfs_ns, query_parents, parent_map, query_labels, label_map)

# ARO:3002870 tet(34) will be forcefully receive ARO:3000012 protein(s) conferring antibiotic resistance via molecular bypass and not tet enzyme

df <- df %>% filter(!(Term_ID %in% "ARO:3002870" & Parent_ID %in% "ARO:3000036"))
df <- df %>% filter(!(Term_ID %in% "ARO:3002870" & Parent_ID %in% "ARO:3000557"))

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
# child ARO:3002484 gets ARO:3000076   || remove ARO:3000076 (classs d) for ARO:3000075 class C
# antibiotic resistant gene variant or mutant ARO:0000031 ---> ARO:3000451 and ARO:3003040
# gene conferring resistance via absence 3003768 (16) protein modulating permeability to antibiotic (41) 3000270

overrul <- c("ARO:3003040", "ARO:3003040", "ARO:3003580", "ARO:3000451", "ARO:3004064", "ARO:3000075", "ARO:3000451", "ARO:3003040", "ARO:3000004", "ARO:3000210", "ARO:3000270")

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

# those repeated after selecting overrul aro, if this had entries, we need to modify overrul (it should be empty)

df_lowest_parent %>% filter(Parent_ID %in% lowest_ontology) %>%  group_by(Term_ID) %>% filter(n()>1)%>% group_by(Term_ID) %>% summarise(n=sum(Parent_ID %in% overrul)) %>% filter(n<1)

# unique lowest for df lowest parent 

df_lowest_parent_2 <- df_lowest_parent %>% filter(Parent_ID %in% lowest_ontology) %>% group_by(Term_ID) %>% group_by(Term_ID) %>% filter(n()<2) %>% ungroup()

# add those with double ontology

df_lowest_parent_2 <- df_lowest_parent_2 %>% bind_rows(df_lowest_parent %>% filter(Parent_ID %in% lowest_ontology) %>% group_by(Term_ID) %>% group_by(Term_ID) %>% filter(n()>1) %>% ungroup() %>% filter(Parent_ID %in% overrul))

# For not lowest
# this are the aros appearing in several branches

repeated_wanted <- df_not_lowest %>% filter(Parent_ID %in% ontologies) %>% group_by(Term_ID) %>% summarise(n = n()) %>% filter(n >1)

# those repeated after selecting overrul aro, if this had entries, we need to modify overrul (it should be empty)

df_not_lowest %>% filter(Parent_ID %in% ontologies) %>%  group_by(Term_ID) %>% filter(n()>1) %>% group_by(Term_ID) %>% summarise(n=sum(Parent_ID %in% overrul)) %>% filter(n<1)

# unique lowest for df lowest parent 

df_not_lowest_2 <- df_not_lowest %>% filter(Parent_ID %in% ontologies) %>% group_by(Term_ID) %>% filter(n()<2) %>% ungroup()

# add those with double ontology

df_not_lowest_2 <- df_not_lowest_2 %>% bind_rows(df_not_lowest %>% filter(Parent_ID %in% ontologies) %>% group_by(Term_ID) %>% group_by(Term_ID) %>% filter(n()>1) %>% ungroup() %>% filter(Parent_ID %in% overrul))
dim(df_not_lowest_2)

df2 <- df_lowest_term %>% bind_rows(df_lowest_parent_2, df_not_lowest_2)
rm(df_lowest_term, df_lowest_parent, df_not_lowest, df_lowest_parent_2, df_not_lowest_2, aros)

# Did we catch all of the aros? (it should be zero)
sum(!df$Term_ID %in% df2$Term_ID)

##############
############## add a new arbitrary level to avoid repetitive onthologies

new_level <- c("chloramphenicol phosphotransferase", "cpt",
               "ciprofloxacin phosphotransferase", "crpP",
               "viomycin phosphotransferase", "vph",
               "major facilitator superfamily (MFS) antibiotic efflux pump", "MFS efflux pump",
               "class C beta-lactamase", "class C beta-lactamase",
               "aminoglycoside bifunctional resistance protein", "bifunctional aminoglycoside",
               "streptothricin acetyltransferase (SAT)", "sat",
               "AAC(2')", "aac",
               "AAC(6')", "aac",
               "AAC(3)", "aac",
               "APH(3')", "aph",
               "APH(6)" , "aph",
               "class A beta-lactamase", "class A beta-lactamase",
               "macrolide phosphotransferase (MPH)", "mph",
               "class B (metallo-) beta-lactamase", "class B beta-lactamase",
               "class D beta-lactamase", "class D beta-lactamase",
               "quinolone resistance protein (qnr)", "qnr",
               "Erm 23S ribosomal RNA methyltransferase", "erm",
               "APH(2'')", "aph",
               "tetracycline inactivation enzyme", "tet enzyme",
               "tetracycline-resistant ribosomal protection protein", "tet RPG",
               "gene(s) or protein(s) associated with a glycopeptide resistance cluster", "van",
               "protein(s) and two-component regulatory system modulating antibiotic efflux", "efflux pump",
               "fosfomycin inactivation enzyme", "fos",
               "antibiotic resistant dihydrofolate reductase", "dfr",
               "rifampin-resistant RNA polymerase-binding protein", "rifampin Rbp",
               "antibiotic resistant gene variant or mutant", "variant or mutant",            
               "gene involved in antibiotic sequestration", "antibiotic sequestration",
               "gene modulating beta-lactam resistance", "beta-lactam modulation resistance",
               "rifampin inactivation enzyme", "rifampin inactivation enzyme",
               "nitroimidazole reductase", "nim",
               "cpa acetyltransferase", "cpa",
               "protein(s) conferring resistance via host-dependent nutrient acquisition", "host-dependent nutrient acquisition",
               "protein modulating permeability to antibiotic", "permeability modulation",
               "chloramphenicol acetyltransferase (CAT)", "cat",
               "streptogramin inactivation enzyme", "vat",
               "gene altering cell wall charge", "cell wall charge",
               "ANT(2'')", "ant",
               "ANT(9)", "ant",
               "ANT(4')", "ant",
               "APH(3'')", "aph",
               "lincosamide nucleotidyltransferase (LNU)", "lnu",
               "ANT(3'')", "ant",
               "gene involved in self-resistance to antibiotic", "self-resistance",
               "ANT(6)", "ant",
               "APH(9)", "aph",
               "APH(7'')", "aph",
               "sulfonamide resistant sul", "sul",
               "ABC-F ATP-binding cassette ribosomal protection protein", "abcF",
               "macrolide glycosyltransferase", "mgt",
               "macrolide esterase", "mel",
               "fusidic acid inactivation enzyme", "fai",
               "Bah amidohydrolase", "bah",
               "Target protecting FusB-type protein conferring resistance to Fusidic acid", "fusB-type",
               "APH(4)", "aph",
               "gene conferring resistance via absence", "resistance by absence",
               "glycopeptide resistance gene cluster", "van",
               "beta-lactam resistant penicillin-binding proteins", "PBP",
               "Edeine acetyltransferase", "edeQ",
               "rifamycin-resistant beta-subunit of RNA polymerase (rpoB)", "rpoB",
               "efflux pump complex or subunit conferring antibiotic resistance", "efflux pump",
               "protein(s) conferring antibiotic resistance via molecular bypass", "molecular bypass",
               "antibiotic target modifying enzyme", "target-modifying enzyme",
               "antibiotic inactivation enzyme", "antibiotic inactivation enzyme",
               "capreomycin phosphotransferase", "cph")

odd_vals  <- new_level[seq(1, length(new_level), by = 2)]
even_vals <- new_level[seq(2, length(new_level), by = 2)]
new_level_df <- data.frame(old = odd_vals, new = even_vals)
rm(new_level, even_vals, odd_vals)

df2 <- df2 %>% mutate(new_level = new_level_df$new[match(Parent_Label, new_level_df$old)])

saveRDS(df2, file = "code_R_analysis/output_abundance_diversity_resistome/conversion_ARO_parent_new_level.rds")
write.csv(df2, file = "code_R_analysis/output_abundance_diversity_resistome/conversion_ARO_parent_new_level.csv",  row.names = F)

fargene <- fargene %>% mutate(manual.parent = df2$Parent_ID[match(manual.ARO, df2$Parent_ID)],
                                      manual.parent_description = df2$Parent_Label[match(manual.ARO, df2$Parent_ID)],
                                      manual.new_level = df2$new_level[match(manual.ARO, df2$Parent_ID)])

fargene.prot <- fargene.prot %>% mutate(manual.parent = df2$Parent_ID[match(manual.ARO, df2$Parent_ID)],
                                                manual.parent_description = df2$Parent_Label[match(manual.ARO, df2$Parent_ID)],
                                                manual.new_level = df2$new_level[match(manual.ARO, df2$Parent_ID)])


################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
### Complement tools 

lst <- list(deeparg.norm = deeparg.norm, deeparg.norm.prot = deeparg.norm.prot, 
            rgi.blast = rgi.blast, rgi.diamond = rgi.diamond, rgi.diamond.prot = rgi.diamond.prot,
            fargene = fargene, fargene.prot = fargene.prot, amrfinder.norm = amrfinder.norm, amrfinder.norm.prot = amrfinder.norm.prot,
            abricate.argannot.norm = abricate.argannot.norm, abricate.card.norm = abricate.card.norm, abricate.megares.norm= abricate.megares.norm,
            abricate.ncbi.norm = abricate.ncbi.norm, abricate.resfinder.norm = abricate.resfinder.norm, resfinder.norm = resfinder.norm)

lst <- lapply(lst, function(x) x %>% mutate(parent = df2$Parent_ID[match(ARO, df2$Term_ID)], 
                                     parent_description = df2$Parent_Label[match(ARO, df2$Term_ID)],
                                     new_level = df2$new_level[match(ARO, df2$Term_ID)]))


# assign new level from fargene when rgi did not report one 
lst$fargene <- lst$fargene %>% 
  mutate(new_level = ifelse(is.na(new_level), manual.new_level, new_level),
         parent = ifelse(is.na(parent), manual.parent, parent),
         parent_description = ifelse(is.na(parent_description), manual.parent_description, parent_description))

lst$fargene.prot <- lst$fargene.prot %>% 
  mutate(new_level = ifelse(is.na(new_level), manual.new_level, new_level),
         parent = ifelse(is.na(parent), manual.parent, parent),
         parent_description = ifelse(is.na(parent_description), manual.parent_description, parent_description))

# change new level for those that are "wrong"
## this part is after assigning new_levels and giving the manual new_level to the observations missing it (it is in retrospective)
# for fargene we took ARO:3002484 and made them all class D
# "ARO:3002600"  AAC(3)-Ib/AAC(6')-Ib3, this is aac(6) no matter what, the coverage is below 45%, id 2 below 70% 
# "ARO:3002599" this is aac(6) no matter what, the coverage is below 40%, id 2 are 100, 6 more below 70% 
# "ARO:3002597"AAC(6')-Ie-APH(2'')-Ia, most coverage are below 40%, 9 below 60%, we give the fargene class instead
# The  MPH that did not agree with fargene have id levels below 30%, they are aph, so very similar, we give the fargene class instead
# The beta-lactam modulation all have id levels below 50%, we give the fargene class instead
# The V/M all have id levels below 40%, we give the fargene class instead
# cph ARO:3007075 is actually a aph2b

lst$fargene <- lst$fargene %>% 
  mutate(ARO = ifelse(coalesce(as.character(new_level), "NA") != coalesce(as.character(manual.new_level), "NA"), manual.ARO, ARO),
         parent = ifelse(coalesce(as.character(new_level), "NA") != coalesce(as.character(manual.new_level), "NA"), manual.parent, parent),
         parent_description = ifelse(coalesce(as.character(new_level), "NA") != coalesce(as.character(manual.new_level), "NA"), manual.parent_description, parent_description)) %>% 
  mutate(new_level = ifelse(coalesce(as.character(new_level), "NA") != coalesce(as.character(manual.new_level), "NA"), manual.new_level, new_level))

lst$fargene.prot <- lst$fargene.prot %>% 
  mutate(ARO = ifelse(coalesce(as.character(new_level), "NA") != coalesce(as.character(manual.new_level), "NA"), manual.ARO, ARO),
         parent = ifelse(coalesce(as.character(new_level), "NA") != coalesce(as.character(manual.new_level), "NA"), manual.parent, parent),
         parent_description = ifelse(coalesce(as.character(new_level), "NA") != coalesce(as.character(manual.new_level), "NA"), manual.parent_description, parent_description)) %>% 
  mutate(new_level = ifelse(coalesce(as.character(new_level), "NA") != coalesce(as.character(manual.new_level), "NA"), manual.new_level, new_level))


rm(list = setdiff(ls(), "lst"))


lst$deeparg.norm.id70 <- lst$deeparg.norm[lst$deeparg.norm$id>=70,]
lst$deeparg.norm.id70$tool <- "DeepARG70"
lst$deeparg.norm.id80 <- lst$deeparg.norm[lst$deeparg.norm$id>=80,]
lst$deeparg.norm.id80$tool <- "DeepARG80"
lst$deeparg.norm.id90 <- lst$deeparg.norm[lst$deeparg.norm$id>=90,]
lst$deeparg.norm.id90$tool <- "DeepARG90"

lst$rgi.diamond.id70 <- lst$rgi.diamond[lst$rgi.diamond$id>=70,]
lst$rgi.diamond.id70$tool <- "RGI-DIAMOND70"
lst$rgi.diamond.id80 <- lst$rgi.diamond[lst$rgi.diamond$id>=80,]
lst$rgi.diamond.id80$tool <- "RGI-DIAMOND80"
lst$rgi.diamond.id90 <- lst$rgi.diamond[lst$rgi.diamond$id>=90,]
lst$rgi.diamond.id90$tool <- "RGI-DIAMOND90"


# save the result of all tools
saveRDS(lst,  file = "code_R_analysis/output_abundance_diversity_resistome/results_tools_all_GMGC.rds", compress = T)

lst$deeparg.norm.id70 <- NULL
lst$deeparg.norm.id80 <- NULL
lst$deeparg.norm.id90 <- NULL
lst$rgi.diamond.id70 <- NULL
lst$rgi.diamond.id80 <- NULL
lst$rgi.diamond.id90 <- NULL


################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
# ABUNDANCES 
# The abundances for each unigene had already been filterd with the file genes_prot_dna.csv

args_abundances <- read.delim("data/abundances/args_abundances.tsv")
metadata <- read.delim("data/metadata_GMGC10.sample.meta.tsv")
metadata <- metadata %>% mutate(sample = sample_id)
args_abundances <- args_abundances %>% left_join(metadata[,c("sample","insertsHQ", "insertsRaw")], by = "sample")
metadata0 <- metadata

# metadata <- metadata %>% filter(!habitat %in% c("amplicon", "isolate", "built-environment"))
# metadata <- metadata %>% filter(sample %in% args_abundances$sample)

abund_habitat <- args_abundances %>% 
  select(c(X,sample)) %>% 
  mutate(habitat = metadata$habitat[match(sample, metadata$sample)]) 

genes_other_habitat <- abund_habitat %>% 
  filter(habitat %in% c("amplicon", "isolate",  "built-environment")) %>% 
  ungroup() %>%
  distinct() %>%
  pull(X)

genes_right_habitat <- abund_habitat %>% 
  ungroup() %>% 
  filter(!habitat %in% c("amplicon", "isolate",  "built-environment")) %>% 
  ungroup() %>%
  distinct() %>%
  pull(X)

genes_right_habitat <- unique(genes_right_habitat)

detected_unigenes_per_habitat <- abund_habitat %>% filter(!is.na(habitat)) %>% group_by(X, habitat) %>% slice_head(n = 1) %>% select(-sample)
write.csv(detected_unigenes_per_habitat, file = "code_R_analysis/output_abundance_diversity_resistome/reported_unigenes_as_ARG_per_habitat.csv", row.names = F)
rm(detected_unigenes_per_habitat)

lst0 <- lst
lst <- lapply(lst, function(x) x %>% filter(query %in% genes_right_habitat))

args_abundances <- args_abundances %>% 
  mutate(habitat = metadata$habitat[match(sample, metadata$sample)]) %>% 
  filter(!habitat %in% c("amplicon", "isolate",  "built-environment")) %>%
  filter(sample  %in% metadata$sample_id)

metadata <- metadata %>% filter(!habitat %in% c("amplicon", "isolate", "built-environment"))


set.seed(2025)

rarefaction <- function(X, raw, inserts, depth = 5e6, seed = 2025){
  #set.seed(seed)
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

rarefaction_fast <- function(X, raw, inserts, depth = 5e6, seed = 2025){
  # set.seed(seed)  # optional
  reads_count <- ceiling(raw)
  arg_reads <- sum(reads_count)
  not_arg_reads <- inserts[1] - arg_reads
  # adjust depth
  depth <- min(depth, inserts[1])
  # probabilities
  probs <- c(reads_count, not_arg_reads) / inserts[1]
  names(probs) <- c(X, "not an ARG")
  # multinomial sampling
  sampled <- rmultinom(1, size = depth, prob = probs)
  return(sampled[,1])
}

# do the rarefaction across samples

# rarefied_counts <- args_abundances %>% filter(raw > 0) %>% ungroup() %>% group_by(sample) %>%
#   summarise(
#     rarefied = list(
#       rarefaction(X, raw, insertsHQ, 5e6)
#     ),
#     .groups = "drop"
#   ) %>%
#   unnest_longer(rarefied, indices_include = TRUE) %>%
#   rename(X = rarefied_id, count = rarefied) %>% 
#   mutate(count = as.integer(count)) %>% 
#   rename(rarified_count = count)

rarefied_counts <- args_abundances %>%
  filter(raw > 0) %>%
  ungroup() %>%
  group_by(sample) %>%
  summarise(
    rarefied = list(
      rarefaction_fast(X, raw, insertsHQ, 5e6)
    ),
    .groups = "drop"
  ) %>%
  unnest_longer(rarefied, indices_include = TRUE) %>%
  rename(X = rarefied_id, rarified_count = rarefied)


# add them to the abundance object   
args_abundances <- args_abundances %>% left_join(rarefied_counts, by = c("sample", "X"))
args_abundances <- args_abundances %>% mutate(rarified_count = ifelse(is.na(rarified_count), 0, rarified_count))

lst$deeparg.norm.id70 <- lst$deeparg.norm[lst$deeparg.norm$id>=70,]
lst$deeparg.norm.id70$tool <- "DeepARG70"
lst$deeparg.norm.id80 <- lst$deeparg.norm[lst$deeparg.norm$id>=80,]
lst$deeparg.norm.id80$tool <- "DeepARG80"
lst$deeparg.norm.id90 <- lst$deeparg.norm[lst$deeparg.norm$id>=90,]
lst$deeparg.norm.id90$tool <- "DeepARG90"

lst$rgi.diamond.id70 <- lst$rgi.diamond[lst$rgi.diamond$id>=70,]
lst$rgi.diamond.id70$tool <- "RGI-DIAMOND70"
lst$rgi.diamond.id80 <- lst$rgi.diamond[lst$rgi.diamond$id>=80,]
lst$rgi.diamond.id80$tool <- "RGI-DIAMOND80"
lst$rgi.diamond.id90 <- lst$rgi.diamond[lst$rgi.diamond$id>=90,]
lst$rgi.diamond.id90$tool <- "RGI-DIAMOND90"

abundance_parent <- function(abund_df, d){
  d <- d %>% filter(!is.na(parent)) 
  Y <- abund_df %>% filter(X %in% d$query)
  Y <- Y %>% mutate(aro = d$ARO[match(X, d$query)], 
                                  parent_description = d$parent_description[match(X, d$query)],
                                  new_level = d$new_level[match(X, d$query)])

  Y_new_level <- Y %>% filter(!is.na(new_level)) %>% group_by(sample, new_level) %>% 
    summarise(
              normed10m = sum(normed10m), 
              distinct_unigenes_rarefied = n_distinct(X[rarified_count > 0]), 
              distinct_unigenes_raw = n_distinct(X)) %>% 
    mutate(tool = d$tool[1]) %>% ungroup() %>% 
    select(sample, new_level, tool, normed10m, 
           distinct_unigenes_rarefied, 
           distinct_unigenes_raw) %>%
    rename(gene = new_level) %>% mutate(aggregation = "new_level") %>% 
    select(sample, gene, aggregation, tool, normed10m, 
           distinct_unigenes_rarefied, 
           distinct_unigenes_raw)
  
#  Y_aro <- Y_aro %>% bind_rows(Y_new_level)
  return(Y_new_level)
}


# calculate the abundance and diversity for all genes and habitats by tool.
lst_abundance_diversity <- do.call(rbind, lapply(lst, function(d) {abundance_parent(args_abundances, d) }))
lst_abundance_diversity0 <- lst_abundance_diversity
lst_abundance_diversity <- lst_abundance_diversity %>% filter(aggregation %in% "new_level") %>% select(-c(aggregation))

# add habitat info
lst_abundance_diversity <- lst_abundance_diversity %>% 
  mutate(habitat = metadata$habitat[match(sample, metadata$sample_id)])

# Unigenes detected as ARG by tool

unigenes <- do.call(rbind, lapply(lst, function(x) x[,c("query", "tool", "ARO", "parent", "parent_description", "new_level", "id")])) 
unigenes <- unigenes %>% filter(query %in% genes_right_habitat)
rownames(unigenes) <- NULL
unigenes <- unigenes %>% select(c("query", "tool","new_level", "id"))


# higher habitat classification for abundance
# lst_abundance_diversity <- lst_abundance_diversity %>% 
#  mutate(habitat2 = SO[lst_abundance_diversity$habitat])

# save abundance and diversity
lst_abundance_diversity <- lst_abundance_diversity %>% mutate(abundance = normed10m/10) %>% 
  rename(richness = distinct_unigenes_rarefied, richness_no_rarified = distinct_unigenes_raw) %>% 
  select(c(sample, gene, tool, abundance, richness, richness_no_rarified))

saveRDS(lst_abundance_diversity, file = "code_R_analysis/output_abundance_diversity_resistome/abundance_diversity.rds", compress = T)
write.csv(lst_abundance_diversity, file = gzfile("code_R_analysis/output_abundance_diversity_resistome/abundance_diversity.csv.gz"), row.names = F)
saveRDS(unigenes, file = "code_R_analysis/output_abundance_diversity_resistome/unigenes_per_tool.rds", compress = T)
write.csv(unigenes, file = gzfile("code_R_analysis/output_abundance_diversity_resistome/unigenes_per_tool.csv.gz"), row.names = F)

# save the results per tool
#for(j in 1:length(lst)){
#  write.csv(lst[[j]], file = paste0("code_R_analysis/output_abundance_diversity_resistome/processed_tool.",names(lst)[j],".csv"), row.names = F)
#}


#### CORE AND PAN 
## load the unigenes clusterd at 90% with vsearch
clusters <- read.delim("cluster_vsearch/clusters.uc", header = F)
clusters <- clusters %>% filter(V1 != "C")
clusters <- clusters %>% mutate(centroid = ifelse(V10 == "*", V9, V10))


# add habitat to abundances
args_abundances <- args_abundances %>% 
  mutate(habitat = metadata$habitat[match(sample, metadata$sample_id)])

# add centroid to abundances
args_abundances <- args_abundances %>% mutate(centroid = clusters$centroid[match(X, clusters$V9)])
args_abundances <- args_abundances %>% rename(query = X)

# add centroid to the result of each tool
lst <- lapply(lst, function(x) x %>% mutate(centroid = clusters$centroid[match(query, clusters$V9)]))

# which centroids have more than 1 gene class in the cluster 

lapply(lst, function(x) x %>% 
         group_by(centroid) %>% mutate(n = n_distinct(new_level))  %>% 
         filter(n>1) %>% arrange(centroid, desc(n)) %>% 
         select(query, centroid, new_level))

# assign majority rule for the class of the centroid per tool
lst <- lapply(lst, function(x) x %>% 
         group_by(centroid) %>% mutate(n = n_distinct(new_level))  %>% 
         mutate(new_level_centroid = new_level[match(centroid, query)]) %>%
         mutate(new_level_majority = names(which.max(table(new_level)))) %>%
         mutate(new_level_centroid = ifelse(is.na(new_level_centroid), new_level_majority, new_level_centroid)))

lst2 <- lapply(lst, function(x) x %>% ungroup() %>% select(-c(centroid, n, new_level_centroid, new_level_majority)))

# save the result of all tools
saveRDS(lst2,  file = "code_R_analysis/output_abundance_diversity_resistome/results_tools.rds", compress = T)

# double_level <- lapply(lst, function(x) x %>%          
#      filter(n>1) %>% arrange(centroid, desc(n)) %>% 
#      ungroup %>%
#      group_by(centroid, new_level) %>%
#      slice_head(n = 1) %>%
#      select(tool, query, centroid, new_level, new_level_centroid, new_level_majority))


# functions for core resistome

cut_size_core <- function(df, cut) {
  return(df %>% ungroup() %>% filter(p >= cut) %>% # filter by proportion of samples with a specific gene
           mutate(cut = cut, cnt = 1) %>% 
           ungroup() %>% 
           select(centroid, new_level_centroid, tool, habitat, cut, cnt))
}

filter_samples_core <- function(args_abundances_core, d, samples_to_collect_number){
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
    #ungroup() %>% group_by(tool, habitat) %>% mutate(N = n_distinct(sample)) %>% # total number of samples per habitat
    ungroup() %>% 
    group_by(tool, habitat) %>% 
    mutate(N = samples_to_collect_number$n[match(habitat, samples_to_collect_number$habitat)]) %>% # total number of samples per habitat
    
    ungroup() %>% group_by(centroid, habitat) %>% mutate(n = n_distinct(sample)) %>% # number of samples where a centroid appears per habitat
    ungroup() %>% mutate( p = n / N) %>% filter(p >= 0.1) %>% # proportion of samples with the centroid
    group_by(habitat, centroid) %>% 
    slice_head(n = 1) %>% ungroup() %>% select(-sample)
    return(DF)
}


core_resistome <- function(args_abundances, samples_to_collect, sed, lst, mx_sample_size, j, cuts, df) {
  set.seed(seed = sed) # random seed
  
  samples_to_collect_number <-  samples_to_collect  %>% ungroup() %>% group_by(habitat) %>%
    slice_sample(n = mx_sample_size,  replace = FALSE) %>% summarise(n = n())
  
  samples_to_collect2 <- samples_to_collect  %>% ungroup() %>% group_by(habitat) %>%
    slice_sample(n = mx_sample_size,  replace = FALSE) %>%
    ungroup() %>% select(sample) %>% pull() 

  args_abundances_core <- args_abundances %>% filter(sample %in% samples_to_collect2) # subsample the metagenomic samples of a habitat
  args_abundances_core <- do.call(rbind, lapply(lst, function(d) filter_samples_core(args_abundances_core, d, samples_to_collect_number))) # find the proportion of samples with the centroid 
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

samples_to_collect <- metadata %>%
  group_by(habitat) %>%
  distinct(sample_id) %>%  # Get distinct values of samples per habitat
  rename(sample = sample_id) %>% 
  mutate(n_unique = n()) %>%  # Calculate the number of unique values in each group
  ungroup() %>%
  group_by(habitat)

df <- data.frame(centroid = NULL, new_level_centroid = NULL, tool = NULL, habitat = NULL, cut = NULL, cnt = NULL)

for(j in 1:length(seeds)){
  print(j)
  df <- core_resistome(args_abundances %>% filter(rarified_count > 0),  samples_to_collect, seeds[j], lst, 100, j, cuts, df)
}

saveRDS(df, file = "code_R_analysis/output_abundance_diversity_resistome/core_resistome.rds", compress = T)
write.csv(df, file = gzfile("code_R_analysis/output_abundance_diversity_resistome/core_resistome.csv.gz"), row.names = F)


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

args_abundances_rarified <- args_abundances %>% filter(rarified_count > 0)

df.pan.list <- vector("list", length(seeds))

for(j in seq_along(seeds)){
  print(j)
  df.pan.list[[j]] <- pan_resistome(args_abundances_rarified, samples_to_collect, seeds[j], lst, 100, j)
}

df.pan2 <- bind_rows(df.pan.list)

saveRDS(df.pan2, file = "code_R_analysis/output_abundance_diversity_resistome/pan_resistome.rds", compress = T)
write.csv(df.pan2, file = gzfile("code_R_analysis/output_abundance_diversity_resistome/pan_resistome.csv.gz"), row.names = F)

