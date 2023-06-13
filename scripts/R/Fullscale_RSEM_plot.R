## Creates an .rds object with all the genes in detected EPS gene clusters, not just the ones detected by the queries

library(tidyverse)
library(ggplot2)
library(readxl)
library(viridis)
library(patchwork)
library(data.table)
library(scales)
library(ggbeeswarm)
`%ni%` <- Negate(`%in%`)

# Preparing data
database_stack <- data.frame(
  Target_label = as.character(),
  Psiblast = as.character()
)

for (file in list.files("/user_data/ahd/EPS_PIPELINE/output_proximity_filtration/psi_operon_full")) {
  temp <- read.csv2(paste0("/user_data/ahd/EPS_PIPELINE/output_proximity_filtration/psi_operon_full/",file), sep = "\t") %>% filter(!is.na(Query_label)) %>% 
    select(Psiblast, Target_label, Query_label, Function, operon) %>% mutate(Psiblast = str_sub(file, end = -5))
  
  database_stack <- database_stack %>% bind_rows(temp)
  
}

database_stack <- database_stack %>% filter(!is.na(Psiblast)) %>% rename(operonNO = operon, gene_id = Target_label, operon = Psiblast, gene_annotation = Query_label)

metadata <- read_xlsx("/mnt/ahd/EPS_PIPELINE/data/metatranscriptomics/metadata_R.xlsx", sheet = 2) %>% 
  group_by(Notes) %>% mutate(id = row_number())

# Loading RSEM data
#data_RSEM <- data.table::fread("/user_data/ahd/EPS_PIPELINE/data/metatranscriptomics/RSEM_TPM_summarised.tsv")

#saveRDS(data_RSEM, file = "/user_data/ahd/EPS_PIPELINE/data/metatranscriptomics/RSEM_TPM_summarised.rds")

data_RSEM <- readRDS(file = "/user_data/ahd/EPS_PIPELINE/data/metatranscriptomics/RSEM_TPM_summarised.rds") %>%
   select(gene_id, `LIB-Glomicave-0189_mRNA.genes.results`:`LIB-Glomicave-0151_resequenced_mRNA.genes.results`)
 
saveRDS(data_RSEM, file = "/user_data/ahd/EPS_PIPELINE/data/metatranscriptomics/RSEM_TPM_EPS_Fullscale.rds")


data_RSEM <- readRDS("/user_data/ahd/EPS_PIPELINE/data/metatranscriptomics/RSEM_TPM_EPS_Fullscale.rds") %>%
  left_join(database_stack, by = "gene_id") %>% filter(!is.na(operon))

saveRDS(data_RSEM, "/user_data/ahd/EPS_PIPELINE/data/metatranscriptomics/RSEM_TPM_EPS_Fullscale_filtered.rds")

data_RSEM <- readRDS("/mnt/ahd/EPS_PIPELINE/data/metatranscriptomics/RSEM_TPM_EPS_Fullscale_filtered.rds")

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
    Function = str_replace(Function, "NA", "Unknown function"),
    operon = str_replace(operon, "alginate", "Alginate"),
    operon = str_replace(operon, "celluloseI", "Cellulose I"),
    operon = str_replace(operon, "celluloseII", "Cellulose II"),
    operon = str_replace(operon, "celluloseIII", "Cellulose III"),
    operon = str_replace(operon, "cellulose_Ac", "Acetylated cellulose"),
    operon = str_replace(operon, "cellulose_NA", "Unclassified cellulose"),
    operon = str_replace(operon, "HA_Pasteurella", "HA (pmHAS)"),
    operon = str_replace(operon, "HA_streptococcus", "HA (has)"),
    operon = str_replace(operon, "NulO_merged", "NulO"),
    operon = str_replace(operon, "pel_merged", "Pel"),
    operon = str_replace(operon, "pnag_pga", "PNAG (pga)"),
    operon = str_replace(operon, "B_subtilis_EPS", "*B. subtilis* EPS"),
    operon = str_replace(operon, "pnag_ica", "PNAG (ica)"),
    operon = str_replace(operon, "xanthan", "Xanthan"),
    operon = str_replace(operon, "psl", "Psl"),
    operon = str_replace(operon, "curdlan", "Curdlan"),
    operon = str_replace(operon, "diutan", "Diutan"),
    operon = str_replace(operon, "succinoglycan", "Succinoglycan"),
    operon = str_replace(operon, "gellan2", "Gellan 2"),
    operon = str_replace(operon, "burkholderia_eps", "*Burkholderia* EPS"),
    operon = str_replace(operon, "amylovoran", "Amylovoran"),
    operon = str_replace(operon, "ColA", "Colanic Acid"),
    operon = str_replace(operon, "salecan", "Salecan"),
    operon = str_replace(operon, "stewartan", "Stewartan"),
    operon = str_replace(operon, "vps", "*Vibrio* EPS"),
    operon = str_replace(operon, "rhizobium_eps", "Rhizobium EPS"),
    operon = str_replace(operon, "gellan1", "Gellan1"),
    operon = str_replace(operon, "acetan", "Acetan"),
    operon = str_replace(operon, "s88", "S88"),
    operon = str_replace(operon, "galactoglucan", "Galactoglucan"),
    operon = str_replace(operon, "levan", "Levan"),
    operon = str_replace(operon, "methanolan", "Methanolan"),
    operon = str_replace(operon, "synechan", "Synechan"),
    operon = str_replace(operon, "B_fragilis_PS_A", "*B. fragilis* PS A"),
    operon = str_replace(operon, "B_fragilis_PS_B", "*B. fragilis* PS B"),
    operon = str_replace(operon, "B_pseudomallei_EPS", "*B. pseudomallei* PS"),
    operon = str_replace(operon, "cepacian", "Cepacian"),
    operon = str_replace(operon, "E_faecalis_PS", "*E. faecalis* PS"),
    operon = str_replace(operon, "emulsan", "Emulsan"),
    operon = str_replace(operon, "EPS273", "EPS273"),
    operon = str_replace(operon, "GG", "GG"),
    operon = str_replace(operon, "glucorhamnan", "Glucorhamnan"),
    operon = str_replace(operon, "L_johnsonii_ATCC_33200_EPS_A", "*L. johnsonii* PS A"),
    operon = str_replace(operon, "L_johnsonii_ATCC_11506_EPS_B", "*L. johnsonii* PS B"),
    operon = str_replace(operon, "L_johnsonii_ATCC_2767_EPS_C", "*L. johnsonii* PS C"),
    operon = str_replace(operon, "L_lactis_EPS", "*L. lactis* PS"),
    operon = str_replace(operon, "L_plantarum_HePS", "*L. plantarum* PS"),
    operon = str_replace(operon, "phosphonoglycan", "Phosphonoglycan")) %>% rename(`Sample#` = id) %>%
  mutate(
    Notes = str_replace(Notes, "Full-scale measurements aerobic stage", "Aerobic stage"),
    Notes = str_replace(Notes, "Full-scale measurements sidestream tank 1", "Sidestream tank 1"),
    Notes = str_replace(Notes, "Full-scale measurements sidestream tank 2", "Sidestream tank 2"),
    Notes = str_replace(Notes, "Full-scale measurements return sludge", "Return sludge"),
    Notes = str_replace(Notes, "Full-scale measurements anoxic stage", "Anoxic stage")) %>%
  mutate(MAG_id = str_sub(gene_id, end = -7)) %>%
  left_join(abundance) %>%
  mutate(rel_TPM = TPM/`Relative abundance (%)`) %>%
  rename(`Process tank` = Notes)



saveRDS(EPS_table_RSEM, "/mnt/ahd/EPS_PIPELINE/data/metatranscriptomics/EPS_table_Fullscale_RSEM.rds")

EPS_table_RSEM <- readRDS("/mnt/ahd/EPS_PIPELINE/data/metatranscriptomics/EPS_table_Fullscale_RSEM.rds")



## Summarise by taking the sum of all instances of a given operon per Process tank per Sample#

EPS_table_RSEM_summarized <- EPS_table_RSEM %>%
  group_by(`Sample#`, `Process tank`, operon) %>%
  summarise(TPM = sum(TPM)) %>% ungroup()

EPS_table_RSEM_summarized_rel <- EPS_table_RSEM %>%
  group_by(`Sample#`, `Process tank`, operon) %>%
  summarise(rel_TPM = sum(rel_TPM)) %>% ungroup()


EPS_table_RSEM_summarized$operon <- factor(EPS_table_RSEM_summarized$operon) %>% fct_reorder(EPS_table_RSEM_summarized$TPM, .desc = TRUE)

EPS_table_RSEM_summarized_rel$operon <- factor(EPS_table_RSEM_summarized_rel$operon) %>% fct_reorder(EPS_table_RSEM_summarized_rel$rel_TPM, .desc = TRUE)

EPS_table_RSEM_summarized$`Process tank` <- factor(EPS_table_RSEM_summarized$`Process tank`,levels = 
                                                        c("Anoxic stage", "Aerobic stage", "Return sludge", "Sidestream tank 1", "Sidestream tank 2"))

EPS_table_RSEM_summarized_rel$`Process tank` <- factor(EPS_table_RSEM_summarized_rel$`Process tank`,levels = 
                                                        c("Anoxic stage", "Aerobic stage", "Return sludge", "Sidestream tank 1", "Sidestream tank 2"))

## Bar plot for presentation

EPS_table_RSEM_summarized_barplot <- EPS_table_RSEM_summarized %>% filter(`Process tank` == "Aerobic stage") %>%
  group_by(operon) %>%
  summarise(
    mean_TPM = mean(TPM),
    sd_TPM = sd(TPM),
    n_TPM = n(),
    se_TPM = sd_TPM/sqrt(n_TPM)
  )

