library( data.table)
library(tidyverse)
#-------------------------------------------------------------------------
#TPM table                                                              |
#-------------------------------------------------------------------------

# For TPM with annotation, detectM data

dt_TPM <- data.table(
  "geneannotation" = as.character()
)
for (file in list.files("data/metatranscriptomics/detectm_geneTPMannot/")) {
  
  temp <- data.table::fread(paste0("data/metatranscriptomics/detectm_geneTPMannot/",file),
                            sep = "\t",
                            header = TRUE
  ) %>% mutate(geneannotation = paste(gene,annotation, sep = "ÆØÅ")) %>% select(geneannotation, TPM)
  
  temp <- setnames(temp, "TPM", file)
  
  dt_TPM <- merge.data.table(dt_TPM, temp, by = "geneannotation", all = TRUE)
  print(file)
  
}


dt_TPM <- dt_TPM %>% separate(geneannotation, into = c("gene", "annotation"), sep = "ÆØÅ")

dt_TPM[is.na(dt_TPM)] <- 0


dt_TPM <- as.data.frame(dt_TPM)


write_tsv(dt_TPM,"data/metatranscriptomics/detectM_TPM_summarised.tsv")


# For TPM with annotation, RSEM data

dt_TPM <- data.table::fread(paste0("data/metatranscriptomics/RSEM/results/LIB-Glomicave-0001_resequenced_mRNA.genes.results"),
                            sep = "\t",
                            header = TRUE
) %>% select(gene_id) %>% column_to_rownames(gene_id)

for (file in list.files("data/metatranscriptomics/RSEM/results/", pattern = "*.genes.results")) {
  
  temp <- data.table::fread(paste0("data/metatranscriptomics/RSEM/results/",file),
                            sep = "\t",
                            header = TRUE
  ) %>% select(gene_id, TPM)
  
  temp <- setnames(temp, "TPM", file)
  
  dt_TPM <- merge.data.table(dt_TPM, temp, by = "gene_id")
  print(file)
  
}


dt_TPM <- dt_TPM %>% separate(geneannotation, into = c("gene", "annotation"), sep = "ÆØÅ")

dt_TPM[is.na(dt_TPM)] <- 0


dt_TPM <- as.data.frame(dt_TPM)


write_tsv(dt_TPM,"data/metatranscriptomics/detectM_TPM_summarised.tsv")

#-------------------------------------------------------------------------
#Forward reads table                                                    |
#-------------------------------------------------------------------------

dt_forward <- data.table(
  "gene" = as.character()
)

for (file in list.files("data/metatranscriptomics/detectm_geneforward/")) {
  
  temp <- data.table::fread(paste0("data/metatranscriptomics/detectm_geneforward/",file),
                            sep = "\t",
                            header = TRUE
  )
  temp <- setnames(temp, "reverse_read_count", file)
  
  dt_forward <- merge.data.table(dt_forward, temp, by = "gene", all = TRUE)
  print(file)
  
}
dt_forward[is.na(dt_forward)] <- 0

write_tsv(dt_forward,"data/metatranscriptomics/detectM_forwardgene_summarised.tsv")
