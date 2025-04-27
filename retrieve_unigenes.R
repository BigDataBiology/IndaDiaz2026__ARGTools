library(dplyr)
library(tidyr)
library(stringr)
library(rdflib)
library(stringr)
library(jsonlite)


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
# produce a list with unique unigenes per tool 
# used to fetch abundances, lengths of genes, etc.

gene.tool <- tibble(data.frame(rbind(cbind("RGI - diamond NT", unique(rgi.diamond$query)),
                                     cbind("RGI - diamond AA", unique(rgi.diamond.prot$query)),
                                     cbind("RGI - blast NT", unique(rgi.blast$query)),
                                     cbind("deepArg - NT", unique(deeparg.norm$query)),
                                     cbind("deepArg - AA", unique(deeparg.norm.prot$query)),
                                     cbind("fargene - NT", unique(fargene$query)),
                                     cbind("fargene - AA", unique(fargene.prot$query)),
                                     cbind("amrfinder NT", unique(amrfinder.norm$query)),
                                     cbind("amrfinder - AA", unique(amrfinder.norm.prot$query)),
                                     cbind("abricate - ARGANNOT NT", unique(abricate.argannot.norm$query)),
                                     cbind("abricate - CARD NT", unique(abricate.card.norm$query)),
                                     cbind("abricate - MEGARES NT", unique(abricate.megares.norm$query)),
                                     cbind("abricate - NCBI NT", unique(abricate.ncbi.norm$query)),
                                     cbind("abricate - ResFinder NT", unique(abricate.resfinder.norm$query)),
                                     cbind("ResFinder - NT", unique(resfinder.norm$query)))))

gene.tool %>% group_by(X1) %>% summarise(n=n()) 
gene.tool <- gene.tool %>% distinct()
gene.tool %>% group_by(X1) %>% summarise(n=n()) 


# unique number of identified genes

length(unique(gene.tool[,2]$X2))

# save genes
write.csv(sort(unique(gene.tool[,2]$X2)), file ="genes_prot_dna.csv", quote=FALSE, row.names=FALSE)

