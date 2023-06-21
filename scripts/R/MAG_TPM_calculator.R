
library(readr)
library(tidyverse)

files <- list.files("/mnt/ahd/EPS_PIPELINE/data/metatranscriptomics/RSEM_new/results", pattern = "*.genes.results")

for (file in files) {
    temp <- read_tsv(paste0("/mnt/ahd/EPS_PIPELINE/data/metatranscriptomics/RSEM_new/results/", file)) %>%
    filter(effective_length > 0) %>%
    mutate(MAG = substring(gene_id, 1, nchar(gene_id) - 6), RKM = expected_count/effective_length) %>% group_by(MAG) %>% 
    summarise(gene_id = gene_id, TPM = TPM, MAG_TPM = ((expected_count/effective_length)/sum(RKM))*1000000) %>%
    write_tsv(temp, paste0("/mnt/ahd/EPS_PIPELINE/data/metatranscriptomics/RSEM_new/results/MAG_TPM_files/", file))
    print(file)
}
