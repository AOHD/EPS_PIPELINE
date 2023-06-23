library(data.table)
library(readr)

files <- list.files("/user_data/ahd/EPS_PIPELINE/data/metatranscriptomics/RSEM_new/results", pattern = "*.genes.results")

for (file in files) {
  temp <- fread(paste0("/user_data/ahd/EPS_PIPELINE/data/metatranscriptomics/RSEM_new/results/", file))
  temp[, MAG := substring(gene_id, 1, nchar(gene_id) - 6)]
  temp[, RKM := ifelse(effective_length == 0, 0, expected_count / effective_length)]
  temp[, MAG_TPM := (RKM / sum(RKM)) * 1000000, by = MAG]
  write.table(temp, paste0("/user_data/ahd/EPS_PIPELINE/data/metatranscriptomics/RSEM_new/MAG_TPM_files/", file, ".tsv"), sep = "\t", row.names = F, quote = F)
  print(file)
}
