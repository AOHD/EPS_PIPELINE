library(tidyverse)
library(ggplot2)
library(readxl)
library(viridis)
library(patchwork)
library(data.table)
library(scales)
library(ggbeeswarm)
library(ggtext)
`%ni%` <- Negate(`%in%`)


# Loading RSEM data
data_RSEM_counts <- data.table::fread("/user_data/ahd/EPS_PIPELINE/data/metatranscriptomics/RSEM_counts_summarised.tsv")


saveRDS(data_RSEM_counts, file = "/user_data/ahd/EPS_PIPELINE/data/metatranscriptomics/RSEM_counts_summarised.rds")

data_RSEM_counts <- readRDS(file = "/user_data/ahd/EPS_PIPELINE/data/metatranscriptomics/RSEM_counts_summarised.rds") %>%
 select(gene_id, `LIB-Glomicave-0189_mRNA.genes.results`:`LIB-Glomicave-0151_resequenced_mRNA.genes.results`)

saveRDS(data_RSEM_counts, file = "/user_data/ahd/EPS_PIPELINE/data/metatranscriptomics/RSEM_counts_EPS_Fullscale.rds")

metadata <- read_xlsx("/user_data/ahd/EPS_PIPELINE/data/metatranscriptomics/metadata_R.xlsx", sheet = 2) %>% 
  group_by(Notes) %>% mutate(id = row_number())

data_RSEM_counts <- readRDS(file = "/user_data/ahd/EPS_PIPELINE/data/metatranscriptomics/RSEM_counts_EPS_Fullscale.rds")

data_RSEM_TPM <- readRDS("/user_data/ahd/EPS_PIPELINE/data/metatranscriptomics/RSEM_TPM_EPS_Fullscale.rds")

table_RSEM_counts <- data_RSEM_counts %>% rename(
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

table_RSEM_TPM <- data_RSEM_TPM %>% rename(
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

gather_RSEM_counts <- table_RSEM_counts %>% gather(2:40,
                                     key = "Sample_ID",
                                     value = "expected_count") %>%
  mutate(Sample_ID = gsub("_","-",Sample_ID))

gather_RSEM_TPM <- table_RSEM_TPM %>% gather(2:40,
                                            key = "Sample_ID",
                                            value = "TPM") %>%
  mutate(Sample_ID = gsub("_","-",Sample_ID))


EPS_table_RSEM_counts <- gather_RSEM_counts %>% left_join(metadata, by = "Sample_ID")
EPS_table_RSEM_TPM <- gather_RSEM_TPM %>% left_join(metadata, by = "Sample_ID")


EPS_table_RSEM_counts <- EPS_table_RSEM_counts %>% rename(`Sample#` = id, `Processing tank` = Notes) %>% 
  mutate(MAG_id = str_sub(gene_id, end = -7))

EPS_table_RSEM_TPM <- EPS_table_RSEM_TPM %>% rename(`Sample#` = id, `Processing tank` = Notes) %>% 
  mutate(MAG_id = str_sub(gene_id, end = -7))

saveRDS(EPS_table_RSEM_counts, file = "/user_data/ahd/EPS_PIPELINE/data/metatranscriptomics/EPS_table_RSEM_counts.rds")
saveRDS(EPS_table_RSEM_TPM, file = "/user_data/ahd/EPS_PIPELINE/data/metatranscriptomics/EPS_table_RSEM_TPM.rds")

EPS_table_RSEM_counts <- readRDS("/user_data/ahd/EPS_PIPELINE/data/metatranscriptomics/EPS_table_RSEM_counts.rds")
EPS_table_RSEM_TPM <- readRDS("/user_data/ahd/EPS_PIPELINE/data/metatranscriptomics/EPS_table_RSEM_counts.rds")


EPS_table_RSEM_summarized_counts <- EPS_table_RSEM_counts %>%
  filter(MAG_id %in% c("Ega_18-Q3-R5-49_MAXAC.001", "AalW_18-Q3-R10-53_BAT3C.524", "Lyne_18-Q3-R50-59_MAXAC.006", "AalE_18-Q3-R2-46_BATAC.251")) %>%
  group_by(`Sample#`, `Processing tank`, MAG_id) %>%
  summarise(expected_count = sum(expected_count)) %>% ungroup()

EPS_table_RSEM_summarized_TPM <- EPS_table_RSEM_TPM %>%
  filter(MAG_id %in% c("Ega_18-Q3-R5-49_MAXAC.001", "AalW_18-Q3-R10-53_BAT3C.524", "Lyne_18-Q3-R50-59_MAXAC.006", "AalE_18-Q3-R2-46_BATAC.251")) %>%
  group_by(`Sample#`, `Processing tank`, MAG_id) %>%
  summarise(expected_count = sum(expected_count)) %>% ungroup()

saveRDS(EPS_table_RSEM_summarized_counts, file = "/user_data/ahd/EPS_PIPELINE/data/metatranscriptomics/EPS_table_RSEM_counts_summarized_selected.rds")
saveRDS(EPS_table_RSEM_summarized_TPM, file = "/user_data/ahd/EPS_PIPELINE/data/metatranscriptomics/EPS_table_RSEM_TPM_summarized_selected.rds")


EPS_table_RSEM_summarized_counts <- readRDS("/user_data/ahd/EPS_PIPELINE/data/metatranscriptomics/EPS_table_RSEM_counts_summarized_selected.rds") %>%
  mutate(
    `Processing tank` = str_replace(`Processing tank`, "Full-scale measurements aerobic stage", "Aerobic stage"),
    `Processing tank` = str_replace(`Processing tank`, "Full-scale measurements sidestream tank 1", "Sidestream tank 1"),
    `Processing tank` = str_replace(`Processing tank`, "Full-scale measurements sidestream tank 2", "Sidestream tank 2"),
    `Processing tank` = str_replace(`Processing tank`, "Full-scale measurements return sludge", "Return sludge"),
    `Processing tank` = str_replace(`Processing tank`, "Full-scale measurements anoxic stage", "Anoxic stage")) %>%
    group_by(`Processing tank`, MAG_id) %>% summarise(
      mean = mean(expected_count),
      sd = sd(expected_count),
      n = n(),
      se = sd/sqrt(n)
    ) %>%
    mutate(MAG_id = ifelse(MAG_id == "Ega_18-Q3-R5-49_MAXAC.001", "*Ca.* P. baldrii", MAG_id),
           MAG_id = ifelse(MAG_id == "AalW_18-Q3-R10-53_BAT3C.524", "*Ca.* P. hodrii", MAG_id),
           MAG_id = ifelse(MAG_id == "Lyne_18-Q3-R50-59_MAXAC.006", "*Ca.* M. subdominans", MAG_id),
           MAG_id = ifelse(MAG_id == "AalE_18-Q3-R2-46_BATAC.251", "midas_g_461 midas_s_461", MAG_id))

EPS_table_RSEM_summarized_TPM <- readRDS("/user_data/ahd/EPS_PIPELINE/data/metatranscriptomics/EPS_table_RSEM_TPM_summarized_selected.rds") %>%
  mutate(
    `Processing tank` = str_replace(`Processing tank`, "Full-scale measurements aerobic stage", "Aerobic stage"),
    `Processing tank` = str_replace(`Processing tank`, "Full-scale measurements sidestream tank 1", "Sidestream tank 1"),
    `Processing tank` = str_replace(`Processing tank`, "Full-scale measurements sidestream tank 2", "Sidestream tank 2"),
    `Processing tank` = str_replace(`Processing tank`, "Full-scale measurements return sludge", "Return sludge"),
    `Processing tank` = str_replace(`Processing tank`, "Full-scale measurements anoxic stage", "Anoxic stage")) %>%
  group_by(`Processing tank`, MAG_id) %>% summarise(
    mean = mean(expected_count),
    sd = sd(expected_count),
    n = n(),
    se = sd/sqrt(n)
  ) %>%
  mutate(MAG_id = ifelse(MAG_id == "Ega_18-Q3-R5-49_MAXAC.001", "*Ca.* P. baldrii", MAG_id),
         MAG_id = ifelse(MAG_id == "AalW_18-Q3-R10-53_BAT3C.524", "*Ca.* P. hodrii", MAG_id),
         MAG_id = ifelse(MAG_id == "Lyne_18-Q3-R50-59_MAXAC.006", "*Ca.* M. subdominans", MAG_id),
         MAG_id = ifelse(MAG_id == "AalE_18-Q3-R2-46_BATAC.251", "midas_g_461 midas_s_461", MAG_id))


summarized_counts_Fullscale <- ggplot(EPS_table_RSEM_summarized, aes(`Processing tank`, mean, color = MAG_id)) +
  geom_point(size = 3) +
  xlab("") + ylab("expected_count")+
  theme_bw() +
  geom_errorbar(
    aes(ymin = mean - se, ymax = mean + se, linewidth = 0.5), width = 0.2, show.legend = FALSE
  ) +
  scale_y_continuous(breaks = pretty_breaks(n=4)) +
  theme(text=element_text(size=30, colour = "black"),
        axis.text.x = element_text(size = 22, colour = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(size = 22, colour = "black"),
        axis.title = element_text(size = 26, face = "bold", colour = "black"),
        legend.text = element_markdown(size = 26, colour = "black"),
        legend.title = element_text(size = 26, face ="bold", colour = "black"),
        legend.position = "right") +
  guides(color = guide_legend(override.aes = list(size = 5),
  ))


ggsave("/user_data/ahd/EPS_PIPELINE/figures/expression/expression_RSEM_Fullscale_population_sum.pdf", summarized_counts_Fullscale, limitsize = FALSE, width = 20, height = 20, dpi = 300)

