library(tidyverse)
library(ggplot2)
library(readxl)
library(viridis)
library(patchwork)
library(data.table)
`%ni%` <- Negate(`%in%`)

#Preparing data
database_stack <- data.frame(
  Target_label = as.character(),
  Psiblast = as.character()
)

for (file in list.files("/mnt/ahd/EPS_PIPELINE/output_proximity_filtration/psi_operon_full")) {
  temp <- read.csv2(paste0("/mnt/ahd/EPS_PIPELINE/output_proximity_filtration/psi_operon_full/",file), sep = "\t") %>% filter(!is.na(Query_label)) %>% select(Psiblast, Target_label, Query_label, Function)

  database_stack <- database_stack %>% bind_rows(temp)

}

database_stack <- database_stack %>% filter(!is.na(Psiblast)) %>% rename(gene_id = Target_label, operon = Psiblast, gene_annotation = Query_label)

metadata <- read_xlsx("/mnt/ahd/EPS_PIPELINE/data/metatranscriptomics/metadata_R.xlsx")

# Loading RSEM data
data_RSEM <- data.table::fread("//wsl$/Ubuntu/mnt/ahd/EPS_PIPELINE/data/metatranscriptomics/RSEM_TPM_summarised.tsv")

saveRDS(data_RSEM, file = "//wsl$/Ubuntu/mnt/ahd/EPS_PIPELINE/data/metatranscriptomics/RSEM_TPM_summarised.rds")

 data_RSEM <- readRDS(file = "/mnt/ahd/data/metatranscriptomics/RSEM_TPM_summarised.rds") %>%
   select(gene_id, `LIB-Glomicave-0060_mRNA.genes.results`:`LIB-Glomicave-0001_resequenced_mRNA.genes.results`) %>%
   left_join(database_stack, by = "gene_id") %>% filter(!is.na(operon))

 saveRDS(data_RSEM, file = "/mnt/ahd/data/metatranscriptomics/RSEM_TPM_EPS.rds")

 data_RSEM <- readRDS(file = "/mnt/ahd/EPS_PIPELINE/data/metatranscriptomics/RSEM_TPM_summarised.rds") %>%
   select(gene_id, `LIB-Glomicave-0110_combined_mRNA.genes.results`:`LIB-Glomicave-0061_mRNA.genes.results`) %>%
   left_join(database_stack, by = "gene_id") %>% filter(!is.na(operon))
 ##

saveRDS(data_RSEM, file = "/mnt/ahd/EPS_PIPELINE/data/metatranscriptomics/RSEM_TPM_EPS_DPAO.rds")

data_RSEM <- readRDS(file = "/mnt/ahd/EPS_PIPELINE/data/metatranscriptomics/RSEM_TPM_EPS.rds")

table_RSEM <- data_RSEM %>% rename(SA_Glomicave_0001 = 61,
                                   SA_Glomicave_0002 = 60,
                                   SA_Glomicave_0003 = 59,
                                   SA_Glomicave_0004 = 58,
                                   SA_Glomicave_0005 = 57,
                                   SA_Glomicave_0006 = 56,
                                   SA_Glomicave_0007 = 55,
                                   SA_Glomicave_0008 = 54,
                                   SA_Glomicave_0009 = 53,
                                   SA_Glomicave_0010 = 52,
                                   SA_Glomicave_0011 = 51,
                                   SA_Glomicave_0012 = 50,
                                   SA_Glomicave_0013 = 49,
                                   SA_Glomicave_0014 = 48,
                                   SA_Glomicave_0015 = 47,
                                   SA_Glomicave_0016 = 46,
                                   SA_Glomicave_0017 = 45,
                                   SA_Glomicave_0018 = 44,
                                   SA_Glomicave_0019 = 43,
                                   SA_Glomicave_0020 = 42,
                                   SA_Glomicave_0021 = 41,
                                   SA_Glomicave_0022 = 40,
                                   SA_Glomicave_0023 = 39,
                                   SA_Glomicave_0024 = 38,
                                   SA_Glomicave_0025 = 37,
                                   SA_Glomicave_0026 = 36,
                                   SA_Glomicave_0027 = 35,
                                   SA_Glomicave_0028 = 34,
                                   SA_Glomicave_0029 = 33,
                                   SA_Glomicave_0030 = 32,
                                   SA_Glomicave_0031 = 31,
                                   SA_Glomicave_0032 = 30,
                                   SA_Glomicave_0033 = 29,
                                   SA_Glomicave_0034 = 28,
                                   SA_Glomicave_0035 = 27,
                                   SA_Glomicave_0036 = 26,
                                   SA_Glomicave_0037 = 25,
                                   SA_Glomicave_0038 = 24,
                                   SA_Glomicave_0039 = 23,
                                   SA_Glomicave_0040 = 22,
                                   SA_Glomicave_0041 = 21,
                                   SA_Glomicave_0042 = 20,
                                   SA_Glomicave_0043 = 19,
                                   SA_Glomicave_0044 = 18,
                                   SA_Glomicave_0045 = 17,
                                   SA_Glomicave_0046 = 16,
                                   SA_Glomicave_0047 = 15,
                                   SA_Glomicave_0048 = 14,
                                   SA_Glomicave_0049 = 13,
                                   SA_Glomicave_0050 = 12,
                                   SA_Glomicave_0051 = 11,
                                   SA_Glomicave_0052 = 10,
                                   SA_Glomicave_0053 = 9,
                                   SA_Glomicave_0054 = 8,
                                   SA_Glomicave_0055 = 7,
                                   SA_Glomicave_0056 = 6,
                                   SA_Glomicave_0057 = 5,
                                   SA_Glomicave_0058 = 4,
                                   SA_Glomicave_0059 = 3,
                                   SA_Glomicave_0060 = 2)


table_RSEM <- data_RSEM %>% rename(SA_Glomicave_0061 = 51,
                                   SA_Glomicave_0062 = 50,
                                   SA_Glomicave_0063 = 49,
                                   SA_Glomicave_0064 = 48,
                                   SA_Glomicave_0065 = 47,
                                   SA_Glomicave_0066 = 46,
                                   SA_Glomicave_0067 = 45,
                                   SA_Glomicave_0068 = 44,
                                   SA_Glomicave_0069 = 43,
                                   SA_Glomicave_0070 = 42,
                                   SA_Glomicave_0071 = 41,
                                   SA_Glomicave_0072 = 40,
                                   SA_Glomicave_0073 = 39,
                                   SA_Glomicave_0074 = 38,
                                   SA_Glomicave_0075 = 37,
                                   SA_Glomicave_0076 = 36,
                                   SA_Glomicave_0077 = 35,
                                   SA_Glomicave_0078 = 34,
                                   SA_Glomicave_0079 = 33,
                                   SA_Glomicave_0080 = 32,
                                   SA_Glomicave_0081 = 31,
                                   SA_Glomicave_0082 = 30,
                                   SA_Glomicave_0083 = 29,
                                   SA_Glomicave_0084 = 28,
                                   SA_Glomicave_0085 = 27,
                                   SA_Glomicave_0086 = 26,
                                   SA_Glomicave_0087 = 25,
                                   SA_Glomicave_0088 = 24,
                                   SA_Glomicave_0089 = 23,
                                   SA_Glomicave_0090 = 22,
                                   SA_Glomicave_0091 = 21,
                                   SA_Glomicave_0092 = 20,
                                   SA_Glomicave_0093 = 19,
                                   SA_Glomicave_0094 = 18,
                                   SA_Glomicave_0095 = 17,
                                   SA_Glomicave_0096 = 16,
                                   SA_Glomicave_0097 = 15,
                                   SA_Glomicave_0098 = 14,
                                   SA_Glomicave_0099 = 13,
                                   SA_Glomicave_0100 = 12,
                                   SA_Glomicave_0101 = 11,
                                   SA_Glomicave_0102 = 10,
                                   SA_Glomicave_0103 = 9,
                                   SA_Glomicave_0104 = 8,
                                   SA_Glomicave_0105 = 7,
                                   SA_Glomicave_0106 = 6,
                                   SA_Glomicave_0107 = 5,
                                   SA_Glomicave_0108 = 4,
                                   SA_Glomicave_0109 = 3,
                                   SA_Glomicave_0110 = 2)

gather_RSEM <- table_RSEM %>% gather(2:61,
                                     key = "Sample_ID",
                                     value = "TPM") %>%
  mutate(Sample_ID = gsub("_","-",Sample_ID)) %>% mutate(Function = as.character(Function))


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
    Function = str_replace(Function, "NA", "Unknown function"))

saveRDS(EPS_table_RSEM, "/mnt/ahd/EPS_PIPELINE/data/raw/EPS_table_PAO_RSEM.rds")

EPS_table_RSEM <- readRDS(file = "/mnt/ahd/EPS_PIPELINE/data/raw/EPS_table_PAO_RSEM.rds")

##Plotting data functions

expressionPlot <- function(EPS_table = EPS_table_RSEM, facet = "Control", replicate = "A", legend = "none") {
  ggplot(EPS_table %>% filter(Replicate == paste(facet, replicate, sep = "_"), operon %ni% c("HA_streptococcus", "NulO_merged", "stewartan")), aes(Time_point, TPM)) +
    theme_bw() +
    ggtitle(paste(facet, replicate, sep = " ")) +
    geom_col(aes(fill = Function)) +
    facet_wrap(~operon) +
    theme(legend.position=legend) +
    ylim(0,12.5)
}

expressionPlot1 <- function(EPS_table = EPS_table_RSEM, facet = "Control", replicate = "A", legend = "none") {
  ggplot(EPS_table %>% filter(Replicate == paste(facet, replicate, sep = "_"), operon %in% c("HA_streptococcus","stewartan")), aes(Time_point, TPM)) +
    theme_bw() +
    ggtitle(paste(facet, replicate, sep = " ")) +
    geom_col(aes(fill = Function)) +
    facet_wrap(~operon) +
    theme(legend.position=legend) +
    ylim(0,160)
}


