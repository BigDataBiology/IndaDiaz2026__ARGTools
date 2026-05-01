# ARG tools 

## Split the initial data

Create smaller `.faa` and `.fna` files to run Snakemake with. 

```bash
cd /work/microbiome/users/juan/arg_compare
python split_data.py /work/microbiome/global_data_spire/GMGC10.data/GMGC10.95nr.faa /work/microbiome/global_data_spire/GMGC10.data/GMGC10.95nr.fna 1000000 &
```

## Running analysis in `.fna` files

### Tools and run

- **fargene** for all the models available
- **deepARG**
- **RGI** with DIAMOND and BLAST
- **Abricate** with ARGANNOT, CARD, MEGARES, NCBI, and RESFINDER
- **Amrfinder**

Resfinder was run in a different folder due to an extended time to debug, and the instructions can be found further down in the document. 

- Running all the tools in the previous list  

```bash
cd /work/microbiome/users/juan/arg_compare/dna
# In a screen run
conda activate snakemake-env
snakemake -s Snakefile_dna -p -j 80 --use-conda
```

The results are found in each of the folders for the tools and are created for each of the "smaller" files.

#### Gather the results and run argAnnot

- RGI

RGI does not need argannot

```bash
cd /work/microbiome/users/juan/arg_compare/dna
# Copy the header
head -n1 rgi_diamond/GMGC10.95nr.fna_block_0000.fna.gz_rgi.tsv >> rgi_diamond.tsv
head -n1 rgi_blast/GMGC10.95nr.fna_block_0000.fna.gz_rgi.tsv >> rgi_blast.tsv
# Merge all the individual results in one file 
cat rgi_diamond/*rgi.tsv | grep -v "^ORF" >> rgi_diamond.tsv
cat rgi_blast/*rgi.tsv | grep -v "^ORF" >> rgi_blast.tsv
```

- amrfinder

```bash
cd /work/microbiome/users/juan/arg_compare/dna
# Copy the header
head -n1 amrfinder/GMGC10.95nr.fna_block_0000.fna.gz.tsv >> amrfinder.tsv
# Merge all the individual results in one file 
cat amrfinder/*.tsv | grep -v "^Protein id" >> amrfinder.tsv
# Manually  replace "Closest reference accession" for  "Accession of closest sequence"
# Obs: this replacement is not needed in the new version of argAnnot
#
# Run argAnnot
conda run -n argnorm argnorm amrfinderplus -i amrfinder.tsv -o amrfinder.norm.tsv
```

- deepArg

```bash
# Copy the header
head -n1 deeparg/GMGC10.95nr.fna_block_0000.fna.gz.ARG >> deeparg.tsv
# Merge all the individual results in one file 
cat deeparg/*gz.ARG | grep -v "^#ARG" >> deeparg.tsv
# Run argAnnot
conda run -n argnorm argnorm deeparg -i deeparg.tsv -o deeparg.norm.tsv
```

- fARGene
Obs: argAnnot does not support fARGene

``` bash
# Merge all the individual results in one file
# HMM has all the scores
# results has the list of queries identified as ARGs 
# predicted have the sequences
cat fargene_hmm/*.txt >> fargene_hmm.txt
cat fargene_results/*.tsv >> fargene_results.tsv
cat fargene_predicted/*.fasta >> fargene_predicted.fasta
```

- abricate 
Obs: argAnnot does not support abricate with CARD database.

```bash
# Copy the header
head -n1 abricate/GMGC10.95nr.fna_block_0000.fna.gz-resfinder.tsv > abricate-resfinder.tsv
head -n1 abricate/GMGC10.95nr.fna_block_0000.fna.gz-megares.tsv > abricate-megares.tsv
head -n1 abricate/GMGC10.95nr.fna_block_0000.fna.gz-argannot.tsv > abricate-argannot.tsv
head -n1 abricate/GMGC10.95nr.fna_block_0000.fna.gz-card.tsv > abricate-card.tsv
head -n1 abricate/GMGC10.95nr.fna_block_0000.fna.gz-ncbi.tsv > abricate-ncbi.tsv
# Merge all the individual results in one file 
cat abricate/*-resfinder.tsv | grep -v "^#FILE" >> abricate-resfinder.tsv
cat abricate/*-megares.tsv | grep -v "^#FILE" >> abricate-megares.tsv
cat abricate/*-argannot.tsv | grep -v "^#FILE" >> abricate-argannot.tsv
cat abricate/*-card.tsv | grep -v "^#FILE" >> abricate-card.tsv
cat abricate/*-ncbi.tsv | grep -v "^#FILE" >> abricate-ncbi.tsv
# Run argAnnot
conda run -n argnorm argnorm abricate -i abricate-resfinder.tsv --db resfinder -o abricate-resfinder.norm.tsv
conda run -n argnorm argnorm abricate -i abricate-megares.tsv --db megares -o abricate-megares.norm.tsv
conda run -n argnorm argnorm abricate -i abricate-argannot.tsv --db argannot -o abricate-argannot.norm.tsv
conda run -n argnorm argnorm abricate -i abricate-ncbi.tsv --db ncbi -o abricate-ncbi.norm.tsv
```

- The RGI, fARGene, and normalized results are exported to be analyzed locally.
- From the local analysis we need to retrieve the name of unique unigenes identified to get their length in aa and nt, and their abundance.

## Running analysis in `.faa` files

### Tools and run

- **fargene** for all the models available
- **deepARG**
- **RGI** with DIAMOND (BLAST had previoulsy been run) 
- **Amrfinder**

- Running all the tools in the previous list  

```bash
cd /work/microbiome/users/juan/arg_compare/protein
# In a screen run
conda activate snakemake-env
snakemake -s Snakefile_protein -p -j 80 --use-conda
```

The results are found in each of the folders for the tools and are created for each of the "smaller" files.

#### Gather the results and run argAnnot

- RGI

RGI does not need argannot

```bash
# Copy the header
head -n1 rgi_diamond/GMGC10.95nr.faa_block_0000.faa.gz_rgi.tsv >> rgi_diamond.tsv
# Merge all the individual results in one file 
cat rgi_diamond/*rgi.tsv | grep -v "^ORF" >> rgi_diamond.tsv
```

- amrfinder 

```bash
# Copy the header
head -n1 amrfinder/GMGC10.95nr.faa_block_0000.faa.gz.tsv >> amrfinder.tsv
# Merge all the individual results in one file 
cat amrfinder/*.tsv | grep -v "^Protein id" >> amrfinder.tsv
# Manually  replace "Closest reference accession" for  "Accession of closest sequence"
# Obs: this replacement is not needed in the new version of argAnnot
#
# Run argAnnot
argnorm amrfinderplus -i amrfinder.tsv -o amrfinder.norm.tsv
```

- deepARG

```bash
# Copy the header
head -n1 deeparg/GMGC10.95nr.faa_block_0000.faa.gz.ARG >> deeparg.tsv
# Merge all the individual results in one file 
cat deeparg/*gz.ARG | grep -v "^#ARG" >> deeparg.tsv
# Run argAnnot
argnorm deeparg -i deeparg.tsv -o deeparg.norm.tsv
```

- fARGene
Obs: argAnnot does not support fARGene

``` bash
# Merge all the individual results in one file
# HMM has all the scores
# results has the list of queries identified as ARGs 
# predicted have the sequences
cat fargene_hmm/*.txt >> fargene_hmm.txt
cat fargene_results/*.tsv >> fargene_results.tsv
cat fargene_predicted/*.fasta >> fargene_predicted.fasta
```


### Running resfinder

Resfinder for `fna` files has its own folder to run snakemake but the results go to `/work/microbiome/users/juan/arg_compare/dna`

```bash
cd /work/microbiome/users/juan/arg_compare/resfinder_dna
# In a screen run
conda activate snakemake-env
snakemake -s Snakefile_resfinder_dna -p -j 80 --use-conda
```

Copy the phenotypes to our working directory to complement the report
```bash
cd /work/microbiome/users/juan/arg_compare/check_missing_annot  
cp /work/microbiome/users/juan/resfinder_databases/resfinder_db/phenotypes.txt .
```

#### Gather the results and run argAnnot


- We cath the files `ResFinder_results_table.txt` in file `resfinder_table.tsv` and run argannot

```bash
# Copy the header
head -n1 resfinder/GMGC10.95nr.fna_block_0031.fna.gz_table.tsv > resfinder_table_tmp.tsv
# Merge all the individual results in one file 
cat resfinder/*table.tsv | tail -n +2 >> resfinder_table_tmp.tsv
# Filter unwanted rows
(head -n 1 resfinder_table_tmp.tsv && tail -n +2 resfinder_table_tmp.tsv | grep -v "^Resistance") > resfinder_table.tsv
rm resfinder_table_tmp.tsv
# Run argAnnot
conda run -n argnorm argnorm resfinder -i resfinder_table.tsv --db resfinder -o resfinder_table.norm.tsv
```

- We cath the files `ResFinder_results_tab.txt` in file `resfinder_tab.tsv` (we will not use )

```bash
touch tmpresfinder.tsv
cat resfinder/*.tsv >> tmpresfinder.tsv
awk '{if (prev && $0 ~ "No hit found") {prev=""; next} if (prev) print prev; prev=$0}' tmpresfinder.tsv > tmpresfinder1.tsv
awk -F"\t" '
NF == 1 { last = $0; next }
NF > 1 { print $0 "\t" last }' tmpresfinder1.tsv | grep -v "Resistance gene" > tmpresfinder2.tsv
grep "Resistance gene" tmpresfinder.tsv |head -n1 | awk '{ print $0 "\tARG class" }' > resfinder_tab.tsv
cat tmpresfinder2.tsv >> resfinder.tsv
rm tmpresfinder1.tsv
rm tmpresfinder2.tsv
rm tmpresfinder.tsv
```

## Retrieve the unigenes flagged as ARGs and complement the ARO annotation

Run the `retrieve_unigenes.R` file, the output should have:
- genes_prot_dna.csv

Use that file in the next section to retrieve abundances and gene lengths

## Retrieve the lengths and abundances 

- Retrieve the names of unigenes identified as ARGs in any of the tools
- `genes_prot_dna.csv` has that list of names 
- Export the file to the server 

Process the `genes_prot_dna.csv` file in the server and extract the lengths 

```bash
cd /work/microbiome/users/juan/arg_compare
# process the csv file and produce a txt file instead
awk -F"," '{print $1}' genes_prot_dna.csv | tail -n +2 |sed 's/\r$//' |sed 's/\t.*//' |sed 's/\r//' | sort | uniq > genes_prot_dna.txt
rm genes_prot_dna.csv
```

### Lengths of GMGC10 sequences

In the server get the lenghts for all GMGC10 unigenes
```bash
cd /work/microbiome/users/juan/arg_compare
# conda activate which seqkit 
# head -n10 /work/microbiome/global_data_spire/GMGC10.data/GMGC10.95nr.fna | seqkit fx2tab --length --name -i 
seqkit fx2tab --length --name -i  /work/microbiome/global_data_spire/GMGC10.data/GMGC10.95nr.fna > length_fna.tsv
seqkit fx2tab --length --name -i  /work/microbiome/global_data_spire/GMGC10.data/GMGC10.95nr.faa > length_faa.tsv
```

Extract the lengths
```bash
# extract the lengths for the identifies genes
awk 'NR==FNR {keys[$0]; next} $1 in keys' genes_prot_dna.txt FS='\t' OFS='\t' length_fna.tsv > arg_length_fna.tsv &
awk 'NR==FNR {keys[$0]; next} $1 in keys' genes_prot_dna.txt FS='\t' OFS='\t' length_faa.tsv > arg_length_faa.tsv &
```

### Abundances 

This is a very slow process
- Extracts 100,000,000 lines at the time
- Grep to the first column of the abundances to the `genes_prot_dna` file
- Saves all the matches in `args_abundances_tmp.tsv`
- Add the header and produce `args_abundances.tsv`

```bash
cd /work/microbiome/users/juan/arg_compare/data/abundances
# open a screen
xzcat /work/microbiome/global_data_spire/GMGC10.data/GMGC10.sample-abundance.tsv.xz | awk -F'\t' '
{
    print > "partition.tsv"
}
NR % 100000000 == 0 {
    close("partition.tsv");
    system("awk -F\"\t\" \047NR==FNR {lookup[$1]; next} $1 in lookup\047 /work/microbiome/users/juan/arg_compare/genes_prot_dna.txt partition.tsv >> args_abundances_tmp.tsv");
    system("rm partition.tsv")
}
END {
    if (NR % 100000000 != 0) {
        close("partition.tsv");
        system("awk -F\"\t\" \047NR==FNR {lookup[$1]; next} $1 in lookup\047 /work/microbiome/users/juan/arg_compare/genes_prot_dna.txt partition.tsv >> args_abundances_tmp.tsv");
        system("rm partition.tsv")
    }
}' 

xzcat /work/microbiome/global_data_spire/GMGC10.data/GMGC10.sample-abundance.tsv.xz | head -n1 > args_abundances.tsv
cat args_abundances_tmp.tsv >> args_abundances.tsv
```

Retrieve all the unigene names that have abundance and see how many detected ARGs are missing abundance 

```bash
# Get all the unique unigene names from the abundance file
xzcat /work/microbiome/global_data_spire/GMGC10.data/GMGC10.sample-abundance.tsv.xz | awk -F'\t' '!seen[$1]++ { print $1 }' > tmp_unigenes_with_abundance.tsv
awk -F'\t' '{print $1}' tmp_unigenes_with_abundance.tsv | sort > unigenes_with_abundance.txt
rm tmp_unigenes_with_abundance.tsv
# sort in place (just doing it again)
sort -o unigenes_with_abundance.txt unigenes_with_abundance.txt
sort -o genes_prot_dna.txt genes_prot_dna.txt
# get detected ARGs without abundance
comm -23 genes_prot_dna.txt unigenes_with_abundance.txt > args_without_abundance.txt
```

Fetch the unigene-environment-no-abundance information for those detected ARGs with no abundance
```bash
# Retrieve the info
awk 'NR==FNR {keys[$0]; next} $1 in keys' args_without_abundance.txt FS='\t' OFS='\t' <(gunzip < /work/microbiome/global_data_spire/GMGC10.data/metadata/GMGC10.unigene-environment.tsv.gz) > unigene-environment-no-abundance.tsv
# Analyze the info

awk '{sum = 0; for (i = 2; i <= NF; i++) sum += $i; print $1"\t"$12"\t"sum}' unigene-environment-no-abundance.tsv | grep "0\t0"  |wc -l
awk '{sum = 0; for (i = 2; i <= NF; i++) sum += $i; print sum}' unigene-environment-no-abundance.tsv | sort | uniq -c 
# From 229,484 detected args, 41,024 do not have any abundance information
# Info from the unigene-environment file says those 41,024 genes are located in:
# N. gene N. envs
#  40895    0
#    128    1
#      1    2
```

## Complementing argAnnot info
- A manual curation of abricate + card is in abricate-card-manual-curation.txt
- resfinder_ARO_mapping.tsv was retrieved from argAnnot
- The curation for abricate with megares and ncbi was done as follows:

```bash
cd /work/microbiome/users/juan/arg_compare/check_missing_annot
mamba run -n rgi rgi load --card_json /work/microbiome/users/juan/rgi/card.json --local
cp   tmp.fna
mamba run -n rgi rgi main -a DIAMOND -i tmp.fna -o abricate_ncbi --local --clean -t contig 
paste <(awk -F"\t" '{print $1}' abricate_ncbi.txt |awk -F" " '{print $1}' | awk '{ sub(/_[^_]*$/, "", $1); print }') abricate_ncbi.txt  > abricate_annotation.tsv

cd /work/microbiome/users/juan/arg_compare/check_missing_annot
cp /work/microbiome/users/juan/e/abricate/db/megares/sequences tmp.fna
mamba run -n rgi rgi main -a DIAMOND -i tmp.fna -o abricate_megares --local --clean -t contig 
paste <(awk -F"\t" '{print $1}' abricate_megares.txt |awk -F" " '{print $1}' | awk '{ sub(/_[^_]*$/, "", $1); print }') abricate_megares.txt  > abricate_annotation_megares.tsv

# fARGene
mamba run -n rgi rgi main -a DIAMOND -i ../dna/fargene_predicted.fasta -o fargene_predicted_fna --local --clean -t contig --include_loose
mamba run -n rgi rgi main -a DIAMOND -i ../protein/fargene_predicted.fasta -o fargene_prot_predicted_faa --local --clean -t protein --include_loose
```

# Clustering ARGs with VSEARCH

```bash
# conda activate which seqkit 
cd /work/microbiome/users/juan/arg_compare/cluster_vsearch
zcat /work/microbiome/global_data_spire/GMGC10.data/GMGC10.95nr.fna | seqkit grep -f ../genes_prot_dna.txt > args.fa
vsearch --cluster_fast args.fa --id 0.90 --centroids centroids.fasta --uc clusters.uc --threads 4 
# conda deactivate
# conda activate which vsearch 
```


## Creating the table of abundance of ARGs

Run the file `retrieve_aros_and_abundances.R`

- `abundances_per_tool_and_unigene.rds`, large table containing sample and abundance per unigene and tool
- `abundance_per_tool.rds`, table containing sample and abundance per **selected (standardized)** ontology and tool, sum of unigene abundance within a specific ontology level
- `abundance_per_tool.csv`
- `results_tools.rds`, list with detection information per tool, complemented with AROs

## Commands to remember 
```bash
python -m resfinder -o . -l 0.6 -t 0.8 --acquired -ifq test*
snakemake -s Snakefile_protein -p -j 80 --use-conda
snakemake -s Snakefile_dna -n -p -j 80 --use-conda
snakemake -s Snakefile_dna -p -j 80 --use-conda

export CGE_RESFINDER_RESGENE_DB="/work/microbiome/users/juan/resfinder_databases/resfinder_db"
export CGE_RESFINDER_RESPOINT_DB="/work/microbiome/users/juan/resfinder_databases/pointfinder_db"
export CGE_DISINFINDER_DB="/work/microbiome/users/juan/resfinder_databases/disinfinder_db"
```


## Counting total number of unigenes in habitats of interest

Ensure that [GMGC10.data](https://git.embl.org/coelho/GMGC10.data) is available
in the current directory (in particular the two files
`metadata/metadata/GMGC10.sample.meta.tsv.gz` and
`GMGC10.sample-abundance.tsv.xz` are needed; you may need to run `git-annex
get` first).


```bash
pixi exec -s python -s pandas python count-unigenes_in_interesting_habitats.py
```
