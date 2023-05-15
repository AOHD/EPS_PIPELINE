library(DESeq2)
library(tximport)
library(readr)
library(here)
library(readxl)
library(data.table)
library(tidyverse)
library(patchwork)
library(scales)
library(viridis)
`%ni%` <- Negate(`%in%`)

#Preparing dataset with EPS-relevant genes
# database_stack <- data.frame(
#   Target_label = as.character(),
#   Psiblast = as.character()
# )
# 
# for (file in list.files("/mnt/ahd/EPS_PIPELINE/output_proximity_filtration/psi_operon_full")) {
#   temp <- read.csv2(paste0("/mnt/ahd/EPS_PIPELINE/output_proximity_filtration/psi_operon_full/",file), sep = "\t") %>% filter(!is.na(Query_label)) %>% select(Psiblast, Target_label, Query_label, Function)
# 
#   database_stack <- database_stack %>% bind_rows(temp)
# 
# }
# 
# database_stack <- database_stack %>% filter(!is.na(Psiblast)) %>% rename(gene_id = Target_label, operon = Psiblast, gene_annotation = Query_label)
# 
# metadata <- read_xlsx("/mnt/ahd/EPS_PIPELINE/data/metatranscriptomics/metadata_R.xlsx")
# 
# #Load dds object
# dds <- readRDS("/mnt/ahd/EPS_PIPELINE/data/metatranscriptomics/DESeq2_EPS_DPAO.rds")
# 
# 
# 
# #Extract and wrangle normalised counts
# 
# dds_counts <-as.data.frame(counts(dds, normalized=T))
# 
# dds_counts <- dds_counts[rowSums(dds_counts[])>0,]
# 
# dds_counts$gene_id <- row.names(dds_counts)
# 
# dds_counts <- dds_counts %>% left_join(database_stack, by = "gene_id")  %>% filter(!is.na(operon)) %>% select(gene_id, everything())
# 
# table_dds_counts <- dds_counts %>% rename(SA_Glomicave_0061 = 2,
#                                           SA_Glomicave_0062 = 3,
#                                           SA_Glomicave_0063 = 4,
#                                           SA_Glomicave_0064 = 5,
#                                           SA_Glomicave_0065 = 6,
#                                           SA_Glomicave_0066 = 7,
#                                           SA_Glomicave_0067 = 8,
#                                           SA_Glomicave_0068 = 9,
#                                           SA_Glomicave_0069 = 10,
#                                           SA_Glomicave_0070 = 11,
#                                           SA_Glomicave_0071 = 12,
#                                           SA_Glomicave_0072 = 13,
#                                           SA_Glomicave_0073 = 14,
#                                           SA_Glomicave_0074 = 15,
#                                           SA_Glomicave_0075 = 16,
#                                           SA_Glomicave_0076 = 17,
#                                           SA_Glomicave_0077 = 18,
#                                           SA_Glomicave_0078 = 19,
#                                           SA_Glomicave_0079 = 20,
#                                           SA_Glomicave_0080 = 21,
#                                           SA_Glomicave_0081 = 22,
#                                           SA_Glomicave_0082 = 23,
#                                           SA_Glomicave_0083 = 24,
#                                           SA_Glomicave_0084 = 25,
#                                           SA_Glomicave_0085 = 26,
#                                           SA_Glomicave_0086 = 27,
#                                           SA_Glomicave_0087 = 28,
#                                           SA_Glomicave_0088 = 29,
#                                           SA_Glomicave_0089 = 30,
#                                           SA_Glomicave_0090 = 31,
#                                           SA_Glomicave_0091 = 32,
#                                           SA_Glomicave_0092 = 33,
#                                           SA_Glomicave_0093 = 34,
#                                           SA_Glomicave_0094 = 35,
#                                           SA_Glomicave_0095 = 36,
#                                           SA_Glomicave_0096 = 37,
#                                           SA_Glomicave_0097 = 38,
#                                           SA_Glomicave_0098 = 39,
#                                           SA_Glomicave_0099 = 40,
#                                           SA_Glomicave_0100 = 41,
#                                           SA_Glomicave_0101 = 42,
#                                           SA_Glomicave_0102 = 43,
#                                           SA_Glomicave_0103 = 44,
#                                           SA_Glomicave_0104 = 45,
#                                           SA_Glomicave_0105 = 46,
#                                           SA_Glomicave_0106 = 47,
#                                           SA_Glomicave_0107 = 48,
#                                           SA_Glomicave_0108 = 49,
#                                           SA_Glomicave_0109 = 50,
#                                           SA_Glomicave_0110 = 51)
# 
# 
# table_dds_counts <- dds_counts %>% rename(SA_Glomicave_0001 = 2,
#                                          SA_Glomicave_0002 = 3,
#                                          SA_Glomicave_0003 = 4,
#                                          SA_Glomicave_0004 = 5,
#                                          SA_Glomicave_0005 = 6,
#                                          SA_Glomicave_0006 = 7,
#                                          SA_Glomicave_0007 = 8,
#                                          SA_Glomicave_0008 = 9,
#                                          SA_Glomicave_0009 = 10,
#                                          SA_Glomicave_0010 = 11,
#                                          SA_Glomicave_0011 = 12,
#                                          SA_Glomicave_0012 = 13,
#                                          SA_Glomicave_0013 = 14,
#                                          SA_Glomicave_0014 = 15,
#                                          SA_Glomicave_0015 = 16,
#                                          SA_Glomicave_0016 = 17,
#                                          SA_Glomicave_0017 = 18,
#                                          SA_Glomicave_0018 = 19,
#                                          SA_Glomicave_0019 = 20,
#                                          SA_Glomicave_0020 = 21,
#                                          SA_Glomicave_0021 = 22,
#                                          SA_Glomicave_0022 = 23,
#                                          SA_Glomicave_0023 = 24,
#                                          SA_Glomicave_0024 = 25,
#                                          SA_Glomicave_0025 = 26,
#                                          SA_Glomicave_0026 = 27,
#                                          SA_Glomicave_0027 = 28,
#                                          SA_Glomicave_0028 = 29,
#                                          SA_Glomicave_0029 = 30,
#                                          SA_Glomicave_0030 = 31,
#                                          SA_Glomicave_0031 = 32,
#                                          SA_Glomicave_0032 = 33,
#                                          SA_Glomicave_0033 = 34,
#                                          SA_Glomicave_0034 = 35,
#                                          SA_Glomicave_0035 = 36,
#                                          SA_Glomicave_0036 = 37,
#                                          SA_Glomicave_0037 = 38,
#                                          SA_Glomicave_0038 = 39,
#                                          SA_Glomicave_0039 = 40,
#                                          SA_Glomicave_0040 = 41,
#                                          SA_Glomicave_0041 = 42,
#                                          SA_Glomicave_0042 = 43,
#                                          SA_Glomicave_0043 = 44,
#                                          SA_Glomicave_0044 = 45,
#                                          SA_Glomicave_0045 = 46,
#                                          SA_Glomicave_0046 = 47,
#                                          SA_Glomicave_0047 = 48,
#                                          SA_Glomicave_0048 = 49,
#                                          SA_Glomicave_0049 = 50,
#                                          SA_Glomicave_0050 = 51,
#                                          SA_Glomicave_0051 = 52,
#                                          SA_Glomicave_0052 = 53,
#                                          SA_Glomicave_0053 = 54,
#                                          SA_Glomicave_0054 = 55,
#                                          SA_Glomicave_0055 = 56,
#                                          SA_Glomicave_0056 = 57,
#                                          SA_Glomicave_0057 = 58,
#                                          SA_Glomicave_0058 = 59,
#                                          SA_Glomicave_0059 = 60,
#                                          SA_Glomicave_0060 = 61)
# gather_dds_counts <- table_dds_counts %>% gather(2:51,
#                                            key = "Sample_ID",
#                                            value = "counts") %>%
#   mutate(Sample_ID = gsub("_","-",Sample_ID)) %>% mutate(Function = as.character(Function))
# 
# 
# EPS_table_dds <- gather_dds_counts %>% left_join(metadata, by = "Sample_ID") %>% mutate(Function = if_else(is.na(Function), "NA", Function)) %>%
#   mutate(
#     Function = str_replace(Function, "MOD", "Modification"),
#     Function = str_replace(Function, "PE", "Polymerization & Export"),
#     Function = str_replace(Function, "GT", "Glycosyl Transferase"),
#     Function = str_replace(Function, "ABC", "ABC transporter"),
#     Function = str_replace(Function, "SY", "Synthase"),
#     Function = str_replace(Function, "PS", "Precursor Synthesis"),
#     Function = str_replace(Function, "REG", "Regulation"),
#     Function = str_replace(Function, "DEG", "Degradation"),
#     Function = str_replace(Function, "TRANS", "Transport"),
#     Function = str_replace(Function, "POL", "Polymerisation"),
#     Function = str_replace(Function, "POLTRANS", "Polymerisation and Transport"),
#     Function = str_replace(Function, "PRE", "Precursor synthesis"),
#     Function = str_replace(Function, "NA", "Unknown function"),
#     N_added = str_replace(N_added, "NA", "Control"))
# 
# saveRDS(EPS_table_dds, "/mnt/ahd/EPS_PIPELINE/data/raw/EPS_table_DPAO_dds.rds")

EPS_table_dds <- readRDS("/mnt/ahd/EPS_PIPELINE/data/raw/EPS_table_DPAO_dds.rds")

##Plotting data functions

expressionPlot <- function(EPS_table = EPS_table_dds, facet = "Control", replicate = "A", legend = "none", ylim = 160) {
  ggplot(EPS_table %>% filter(Replicate == paste(facet, replicate, sep = "_"), operon %ni% c("HA_streptococcus", "NulO_merged", "stewartan")), aes(Time_point, counts)) +
    theme_bw() +
    ggtitle(paste(facet, replicate, sep = " ")) +
    geom_col(aes(fill = Function)) +
    facet_wrap(~operon) +
    theme(legend.position=legend) +
    ylim(0,ylim)
}

expressionPlot1 <- function(EPS_table = EPS_table_dds, facet = "Control", replicate = "A", legend = "none", ylim = 1000) {
  ggplot(EPS_table %>% filter(Replicate == paste(facet, replicate, sep = "_"), operon %in% c("HA_streptococcus", "stewartan")), aes(Time_point, counts)) +
    theme_bw() +
    ggtitle(paste(facet, replicate, sep = " ")) +
    geom_col(aes(fill = Function)) +
    facet_wrap(~operon) +
    theme(legend.position=legend) +
    ylim(0,ylim)
}




