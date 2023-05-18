library(tidyverse)
library(ggplot2)
library(readxl)
library(viridis)
library(patchwork)
library(data.table)
library(scales)
library(ggbeeswarm)
library(ggtext)
library("ggnewscale")
`%ni%` <- Negate(`%in%`)


metadata <- read_xlsx("/user_data/ahd/EPS_PIPELINE/data/metatranscriptomics/metadata_R.xlsx", sheet = 2) %>% 
  group_by(Notes) %>% mutate(id = row_number())


data_RSEM_TPM <- readRDS("/user_data/ahd/EPS_PIPELINE/data/metatranscriptomics/RSEM_TPM_EPS_Fullscale.rds")


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

gather_RSEM_TPM <- table_RSEM_TPM %>% gather(2:40,
                                            key = "Sample_ID",
                                            value = "TPM") %>%
  mutate(Sample_ID = gsub("_","-",Sample_ID))


EPS_table_RSEM_counts <- gather_RSEM_counts %>% left_join(metadata, by = "Sample_ID")
EPS_table_RSEM_TPM <- gather_RSEM_TPM %>% left_join(metadata, by = "Sample_ID")

EPS_table_RSEM_TPM <- EPS_table_RSEM_TPM %>% rename(`Sample#` = id, `Processing tank` = Notes) %>% 
  mutate(MAG_id = str_sub(gene_id, end = -7))

saveRDS(EPS_table_RSEM_TPM, file = "/user_data/ahd/EPS_PIPELINE/data/metatranscriptomics/EPS_table_RSEM_TPM.rds")


EPS_table_RSEM_summarized_TPM <- EPS_table_RSEM_TPM %>%
  filter(MAG_id %in% c("Ega_18-Q3-R5-49_MAXAC.001", "AalW_18-Q3-R10-53_BAT3C.524", "Lyne_18-Q3-R50-59_MAXAC.006", "AalE_18-Q3-R2-46_BATAC.251")) %>%
  group_by(`Sample#`, `Processing tank`, MAG_id) %>%
  summarise(TPM = sum(TPM)) %>% ungroup()

saveRDS(EPS_table_RSEM_summarized_TPM, file = "/user_data/ahd/EPS_PIPELINE/data/metatranscriptomics/EPS_table_RSEM_TPM_summarized_selected.rds")

EPS_table_RSEM_summarized_TPM <- readRDS("/mnt/ahd/EPS_PIPELINE/data/metatranscriptomics/EPS_table_RSEM_TPM_summarized_selected.rds")
EPS_table_RSEM <- readRDS("/mnt/ahd/EPS_PIPELINE/data/metatranscriptomics/EPS_table_Fullscale_RSEM.rds")


EPS_table_RSEM_summarized_TPM <- EPS_table_RSEM_summarized_TPM %>%
  group_by(`Processing tank`, MAG_id) %>% summarise(
    mean_TPM = mean(TPM),
    sd_TPM = sd(TPM),
    n_TPM = n(),
    se_TPM = sd_TPM/sqrt(n_TPM)
  ) %>%
  mutate(MAG_id = ifelse(MAG_id == "Ega_18-Q3-R5-49_MAXAC.001", "*Ca.* P. baldrii", MAG_id),
         MAG_id = ifelse(MAG_id == "AalW_18-Q3-R10-53_BAT3C.524", "*Ca.* P. hodrii", MAG_id),
         MAG_id = ifelse(MAG_id == "Lyne_18-Q3-R50-59_MAXAC.006", "*Ca.* M. subdominans", MAG_id),
         MAG_id = ifelse(MAG_id == "AalE_18-Q3-R2-46_BATAC.251", "midas_g_461 midas_s_461", MAG_id)) %>% 
  mutate(`Processing tank` = ifelse(`Processing tank` == "Full-scale measurements anoxic stage", "Anoxic", `Processing tank`),
                                 `Processing tank` = ifelse(`Processing tank` == "Full-scale measurements aerobic stage", "Aerobic", `Processing tank`),
                                 `Processing tank` = ifelse(`Processing tank` == "Full-scale measurements return sludge", "RS", `Processing tank`),
                                 `Processing tank` = ifelse(`Processing tank` == "Full-scale measurements sidestream tank 1", "ST1", `Processing tank`),
                                 `Processing tank` = ifelse(`Processing tank` == "Full-scale measurements sidestream tank 2", "ST2", `Processing tank`))
EPS_table_RSEM_summarized_TPM$`Processing tank` <- factor(EPS_table_RSEM_summarized_TPM$`Processing tank`,levels = 
                                             c("Anoxic", "Aerobic", "RS", "ST1", "ST2"))


EPS_table_operons <- EPS_table_RSEM %>% filter(MAG_id %in% c("Ega_18-Q3-R5-49_MAXAC.001", "AalW_18-Q3-R10-53_BAT3C.524", "Lyne_18-Q3-R50-59_MAXAC.006", "AalE_18-Q3-R2-46_BATAC.251")) %>% 
  group_by(`Sample#`,`Processing tank`, MAG_id, operon) %>% 
  summarise(TPM = sum(TPM)) %>% ungroup() %>%
  group_by(`Processing tank`, MAG_id, operon) %>% 
  summarise(
    mean_TPM_operons = mean(TPM),
    sd_TPM_operons = sd(TPM),
    n_TPM_operons = n(),
    se_TPM_operons = sd_TPM_operons/sqrt(n_TPM_operons)
  ) %>%
  mutate(MAG_id = ifelse(MAG_id == "Ega_18-Q3-R5-49_MAXAC.001", "*Ca.* P. baldrii", MAG_id),
         MAG_id = ifelse(MAG_id == "AalW_18-Q3-R10-53_BAT3C.524", "*Ca.* P. hodrii", MAG_id),
         MAG_id = ifelse(MAG_id == "Lyne_18-Q3-R50-59_MAXAC.006", "*Ca.* M. subdominans", MAG_id),
         MAG_id = ifelse(MAG_id == "AalE_18-Q3-R2-46_BATAC.251", "midas_g_461 midas_s_461", MAG_id)) %>%
  mutate(`Processing tank` = ifelse(`Processing tank` == "Anoxic stage", "Anoxic", `Processing tank`),
         `Processing tank` = ifelse(`Processing tank` == "Aerobic stage", "Aerobic", `Processing tank`),
         `Processing tank` = ifelse(`Processing tank` == "Return sludge", "RS", `Processing tank`),
         `Processing tank` = ifelse(`Processing tank` == "Sidestream tank 1", "ST1", `Processing tank`),
         `Processing tank` = ifelse(`Processing tank` == "Sidestream tank 2", "ST2", `Processing tank`))

EPS_table_operons$`Processing tank` <- factor(EPS_table_operons$`Processing tank`,levels = 
                                             c("Anoxic", "Aerobic", "RS", "ST1", "ST2"))

operon_plot <- ggplot(EPS_table_operons, aes(`Processing tank`, mean_TPM_operons, color = MAG_id, shape = operon)) +
  geom_pointrange(size = 1, aes(ymin = mean_TPM_operons - se_TPM_operons, ymax = mean_TPM_operons + se_TPM_operons)) +
  geom_line(aes(group = interaction(MAG_id, operon))) +
  xlab("") + ylab("TPM") +
  theme_bw() +
  scale_y_continuous(trans = "log10") +
  theme(text=element_text(size=30, colour = "black"),
        axis.text.x = element_text(size = 22, colour = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(size = 22, colour = "black"),
        axis.title = element_blank(),
        legend.text = element_markdown(size = 26, colour = "black"),
        legend.title = element_text(size = 26, face ="bold", colour = "black"),
        legend.position = "right") + labs(color = "exoPS", shape = "MAG")

species_plot <- ggplot(EPS_table_RSEM_summarized_TPM, aes(`Processing tank`, mean_TPM, color = MAG_id)) +
  geom_pointrange(size = 1, aes(ymin = mean_TPM - se_TPM, ymax = mean_TPM + se_TPM)) +
  geom_line(aes(group = MAG_id)) +
  xlab("") + ylab("TPM") +
  theme_bw() +
  scale_y_continuous(trans = "log10") +
  theme(text=element_text(size=30, colour = "black"),
        axis.text.x = element_text(size = 22, colour = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(size = 22, colour = "black"),
        axis.title = element_text(size = 26, face = "bold", colour = "black"),
        legend.text = element_markdown(size = 26, colour = "black"),
        legend.title = element_text(size = 26, face ="bold", colour = "black"),
        legend.position = "none")

species_operon_plot <- species_plot + operon_plot

ggsave("/mnt/ahd/EPS_PIPELINE/figures/expression/expression_RSEM_Fullscale_species_operon.png", species_operon_plot, limitsize = FALSE, width = 30, height = 20, dpi = 300)

#Operon expression normalized by MAG TPM


EPS_table_normalized <- left_join(EPS_table_operons, EPS_table_RSEM_summarized_TPM, by = c("Processing tank", "MAG_id")) %>% 
  select(`Processing tank`, MAG_id, operon, mean_TPM_operons, mean_TPM) %>%
  mutate(normalized_TPM = mean_TPM_operons/mean_TPM)


normalized_plot <- ggplot(EPS_table_normalized, aes(`Processing tank`, normalized_TPM, color = MAG_id, shape = operon)) +
  geom_point(size = 1) +
  geom_line(aes(group = interaction(MAG_id, operon))) +
  xlab("") + ylab("TPM") +
  theme_bw() +
  scale_y_continuous(trans = "log10") +
  theme(text=element_text(size=30, colour = "black"),
        axis.text.x = element_text(size = 22, colour = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(size = 22, colour = "black"),
        axis.title = element_blank(),
        legend.text = element_markdown(size = 26, colour = "black"),
        legend.title = element_text(size = 26, face ="bold", colour = "black"),
        legend.position = "right") + labs(color = "exoPS", shape = "MAG")


