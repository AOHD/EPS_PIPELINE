library(data.table)
library(readr)
library(tidyverse)

files <- list.files("/mnt/ahd/EPS_PIPELINE/data/metatranscriptomics/RSEM_new/results", pattern = "*.genes.results")

for (file in files) {
    temp <- read_tsv(paste0("/mnt/ahd/EPS_PIPELINE/data/metatranscriptomics/RSEM_new/results/", file)) %>% 
    mutate(MAG = substring(gene_id, 1, nchar(gene_id) - 6),
           RKM = ifelse(effective_length == 0, 0, expected_count/effective_length)) %>% 
    group_by(MAG) %>% 
    reframe(gene_id = gene_id, TPM = TPM, MAG_TPM = ((expected_count/effective_length)/(sum(RKM)))*1000000) %>%
    write_tsv(paste0("/mnt/ahd/EPS_PIPELINE/data/metatranscriptomics/RSEM_new/MAG_TPM_files/", file))
    print(file)
}

for (file in files) {
  temp <- fread(paste0("/mnt/ahd/EPS_PIPELINE/data/metatranscriptomics/RSEM_new/results/", file))
  temp[, MAG := substring(gene_id, 1, nchar(gene_id) - 6)]
  temp[, RKM := ifelse(effective_length == 0, 0, expected_count / effective_length)]
  temp[, MAG_TPM := ((expected_count / effective_length) / sum(RKM)) * 1000000, by = MAG]
  fwrite(temp, paste0("/mnt/ahd/EPS_PIPELINE/data/metatranscriptomics/RSEM_new/MAG_TPM_files/", file))
  print(file)
}
