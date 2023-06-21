
library(readr)
library(tidyverse)

files <- list.files("/user_data/ahd/EPS_PIPELINE/data/metatranscriptomics/RSEM_new/results", pattern = "*.genes.results")

for (file in files) {
    temp <- read_tsv(paste0("/user_data/ahd/EPS_PIPELINE/data/metatranscriptomics/RSEM_new/results", file)) %>% 
    mutate(MAG = substring(gene_id, 1, nchar(gene_id) - 6)) %>% group_by(MAG) %>% 
    summarise(gene_id = gene_id, MAG_TPM = ((effective_count/effective_length)/(sum(effective_count/effective_length)))*1000000) %>%
    write_tsv(paste0("/user_data/ahd/EPS_PIPELINE/data/metatranscriptomics/RSEM_new/results/MAG_TPM_files/", file))
    print(file)
}