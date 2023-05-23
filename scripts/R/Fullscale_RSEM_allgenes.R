library(tidyverse)
library(ggplot2)
library(readxl)
library(viridis)
library(patchwork)
library(data.table)
library(scales)
library(stringr)
`%ni%` <- Negate(`%in%`)

# Preparing data
database_stack <- data.frame(
  Target_label = as.character(),
  Psiblast = as.character()
)

for (file in list.files("/mnt/ahd/EPS_PIPELINE/output_proximity_filtration/psi_operon_full")) {
  temp <- read.csv2(paste0("/mnt/ahd/EPS_PIPELINE/output_proximity_filtration/psi_operon_full/",file), sep = "\t") %>% 
    select(Psiblast, Target_label, Query_label, Function, operon) %>% mutate(Psiblast = str_sub(file, end = -5))
  
  database_stack <- database_stack %>% bind_rows(temp)
  
}

database_stack <- database_stack %>% rename(operonNO = operon, gene_id = Target_label, operon = Psiblast, gene_annotation = Query_label) %>% mutate(EPS = "x")

metadata <- read_xlsx("/mnt/ahd/EPS_PIPELINE/data/metatranscriptomics/metadata_R.xlsx", sheet = 2) %>% 
  group_by(Notes) %>% mutate(id = row_number())

# Loading RSEM data
# data_RSEM <- data.table::fread("//wsl$/Ubuntu/user_data/ahd/EPS_PIPELINE/data/metatranscriptomics/RSEM_TPM_summarised.tsv")
# 
# saveRDS(data_RSEM, file = "//wsl$/Ubuntu/user_data/ahd/EPS_PIPELINE/data/metatranscriptomics/RSEM_TPM_summarised.rds")

data_RSEM <- readRDS(file = "/mnt/ahd/EPS_PIPELINE/data/metatranscriptomics/RSEM_TPM_summarised.rds") %>%
  select(gene_id, `LIB-Glomicave-0189_mRNA.genes.results`:`LIB-Glomicave-0151_resequenced_mRNA.genes.results`) %>%
  left_join(database_stack, by = "gene_id") %>% filter(!is.na(EPS))

saveRDS(data_RSEM, file = "/mnt/ahd/EPS_PIPELINE/data/metatranscriptomics/RSEM_TPM_EPS_Fullscale_All.rds")

data_RSEM <- readRDS("/mnt/ahd/EPS_PIPELINE/data/metatranscriptomics/RSEM_TPM_EPS_Fullscale_All.rds")


table_RSEM <- data_RSEM %>% rename(
  SA_Glomicave_0151 = 40,
  SA_Glomicave_0152 = 39,
  SA_Glomicave_0153 = 38,
  SA_Glomicave_0154 = 37,
  SA_Glomicave_0155 = 36,
  SA_Glomicave_0156 = 35,
  SA_Glomicave_0157 = 34,
  SA_Glomicave_0158 = 33,
  SA_Glomicave_0159 = 32,
  SA_Glomicave_0160 = 31,
  SA_Glomicave_0161 = 30,
  SA_Glomicave_0162 = 29,
  SA_Glomicave_0163 = 28,
  SA_Glomicave_0164 = 27,
  SA_Glomicave_0165 = 26,
  SA_Glomicave_0166 = 25,
  SA_Glomicave_0167 = 24,
  SA_Glomicave_0168 = 23,
  SA_Glomicave_0169 = 22,
  SA_Glomicave_0170 = 21,
  SA_Glomicave_0171 = 20,
  SA_Glomicave_0172 = 19,
  SA_Glomicave_0173 = 18,
  SA_Glomicave_0174 = 17,
  SA_Glomicave_0175 = 16,
  SA_Glomicave_0176 = 15,
  SA_Glomicave_0177 = 14,
  SA_Glomicave_0178 = 13,
  SA_Glomicave_0179 = 12,
  SA_Glomicave_0180 = 11,
  SA_Glomicave_0181 = 10,
  SA_Glomicave_0182 = 9,
  SA_Glomicave_0183 = 8,
  SA_Glomicave_0184 = 7,
  SA_Glomicave_0185 = 6,
  SA_Glomicave_0186 = 5,
  SA_Glomicave_0187 = 4,
  SA_Glomicave_0188 = 3,
  SA_Glomicave_0189 = 2)

gather_RSEM <- table_RSEM %>% gather(2:40,
                                     key = "Sample_ID",
                                     value = "TPM") %>%
  mutate(Sample_ID = gsub("_","-",Sample_ID)) %>% mutate(Function = as.character(Function))

abundance <- readRDS("/mnt/ahd/EPS_PIPELINE/data/metatranscriptomics/abundance_metaG.rds")

EPS_table_RSEM <- gather_RSEM %>% left_join(metadata, by = "Sample_ID") %>% mutate(Function = if_else(is.na(Function), "NA", Function)) %>%
  mutate(
    Function = str_replace(Function, "MOD", "Modification"),
    Function = str_replace(Function, "PE", "Polymerization & Export"),
    Function = str_replace(Function, "GT", "Glycosyl Transferase"),
    Function = str_replace(Function, "ABC", "ABC transporter"),
    Function = str_replace(Function, "SY", "Synthase"),
    Function = str_replace(Function, "PS", "Precursor Synthesis"),
    Function = str_replace(Function, "REG", "Regulation"),
    Function = str_replace(Function, "DEG", "Degradation"),
    Function = str_replace(Function, "TRANS", "Transport"),
    Function = str_replace(Function, "POL", "Polymerisation"),
    Function = str_replace(Function, "POLTRANS", "Polymerisation and Transport"),
    Function = str_replace(Function, "PRE", "Precursor synthesis"),
    Function = str_replace(Function, "NA", "Unknown function")) %>%
  mutate(
    Notes = str_replace(Notes, "Full-scale measurements aerobic stage", "Aerobic stage"),
    Notes = str_replace(Notes, "Full-scale measurements sidestream tank 1", "Sidestream tank 1"),
    Notes = str_replace(Notes, "Full-scale measurements sidestream tank 2", "Sidestream tank 2"),
    Notes = str_replace(Notes, "Full-scale measurements return sludge", "Return sludge"),
    Notes = str_replace(Notes, "Full-scale measurements anoxic stage", "Anoxic stage")) %>% rename(Target_label = gene_id) %>%
  mutate(MAG_id = str_sub(Target_label, end = -7)) %>%
  left_join(abundance) %>%
  mutate(rel_TPM = TPM/`Relative abundance (%)`) %>%
  rename(`Processing tank` = Notes) %>%
  mutate(operonNO = paste0(operon, "_", operonNO))

saveRDS(EPS_table_RSEM, "/mnt/ahd/EPS_PIPELINE/data/metatranscriptomics/EPS_table_Fullscale_RSEM_All.rds")

EPS_table_RSEM_All <- readRDS("/mnt/ahd/EPS_PIPELINE/data/metatranscriptomics/EPS_table_Fullscale_RSEM_All.rds")



