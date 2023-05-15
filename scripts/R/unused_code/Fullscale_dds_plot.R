library(DESeq2)
library(tximport)
library(readr)
library(here)
library(readxl)
library(data.table)
library(tidyverse)
library(patchwork)
library(scales)
library(stringr)
library(viridis)
library(forcats)
library(ggbeeswarm)
`%ni%` <- Negate(`%in%`)
#Preparing dataset with EPS-relevant genes
 database_stack <- data.frame(
   Target_label = as.character(),
   Psiblast = as.character()
 )
 
for (file in list.files("/mnt/ahd/EPS_PIPELINE/output_proximity_filtration/psi_operon_full")) {
  temp <- read.csv2(paste0("/mnt/ahd/EPS_PIPELINE/output_proximity_filtration/psi_operon_full/",file), sep = "\t") %>% filter(!is.na(Query_label)) %>% 
    select(Psiblast, Target_label, Query_label, Function) %>% mutate(Psiblast = str_sub(file, end = -5))
  
  database_stack <- database_stack %>% bind_rows(temp)
  
}

database_stack <- database_stack %>% filter(!is.na(Psiblast)) %>% rename(gene_id = Target_label, operon = Psiblast, gene_annotation = Query_label)
 

metadata <- read_xlsx("/mnt/ahd/EPS_PIPELINE/data/metatranscriptomics/metadata_R.xlsx", sheet = 2) %>% 
  group_by(Notes) %>% mutate(id = row_number())


 #Load dds object
 dds <- readRDS("/mnt/ahd/EPS_PIPELINE/data/metatranscriptomics/DESeq2_EPS_Fullscale.rds")
 
 
 
 #Extract and wrangle normalised counts
 
 dds_counts <-as.data.frame(counts(dds, normalized=T))
 
 dds_counts <- dds_counts[rowSums(dds_counts[])>0,]
 
 dds_counts$gene_id <- row.names(dds_counts)
 
 dds_counts <- dds_counts %>% left_join(database_stack, by = "gene_id")  %>% filter(!is.na(operon)) %>% select(gene_id, everything())
 
 table_dds_counts <- dds_counts %>% rename(
                                   SA_Glomicave_0151 = 2,
                                   SA_Glomicave_0152 = 3,
                                   SA_Glomicave_0153 = 4,
                                   SA_Glomicave_0154 = 5,
                                   SA_Glomicave_0155 = 6,
                                   SA_Glomicave_0156 = 7,
                                   SA_Glomicave_0157 = 8,
                                   SA_Glomicave_0158 = 9,
                                   SA_Glomicave_0159 = 10,
                                   SA_Glomicave_0160 = 11,
                                   SA_Glomicave_0161 = 12,
                                   SA_Glomicave_0162 = 13,
                                   SA_Glomicave_0163 = 14,
                                   SA_Glomicave_0164 = 15,
                                   SA_Glomicave_0165 = 16,
                                   SA_Glomicave_0166 = 17,
                                   SA_Glomicave_0167 = 18,
                                   SA_Glomicave_0168 = 19,
                                   SA_Glomicave_0169 = 20,
                                   SA_Glomicave_0170 = 21,
                                   SA_Glomicave_0171 = 22,
                                   SA_Glomicave_0172 = 23,
                                   SA_Glomicave_0173 = 24,
                                   SA_Glomicave_0174 = 25,
                                   SA_Glomicave_0175 = 26,
                                   SA_Glomicave_0176 = 27,
                                   SA_Glomicave_0177 = 28,
                                   SA_Glomicave_0178 = 29,
                                   SA_Glomicave_0179 = 30,
                                   SA_Glomicave_0180 = 31,
                                   SA_Glomicave_0181 = 32,
                                   SA_Glomicave_0182 = 33,
                                   SA_Glomicave_0183 = 34,
                                   SA_Glomicave_0184 = 35,
                                   SA_Glomicave_0185 = 36,
                                   SA_Glomicave_0186 = 37,
                                   SA_Glomicave_0187 = 38,
                                   SA_Glomicave_0188 = 39,
                                   SA_Glomicave_0189 = 40)

 gather_dds_counts <- table_dds_counts %>% gather(2:40,
                                            key = "Sample_ID",
                                            value = "counts") %>%
   mutate(Sample_ID = gsub("_","-",Sample_ID)) %>% mutate(Function = as.character(Function))
 
 
 EPS_table_dds <- gather_dds_counts %>% left_join(metadata, by = "Sample_ID") %>% mutate(Function = if_else(is.na(Function), "NA", Function)) %>%
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
    operon = str_replace(operon, "pnag_eps", "PNAG (eps)"),
    operon = str_replace(operon, "pnag_ica", "PNAG (ica)"),
    operon = str_replace(operon, "xanthan", "Xanthan"),
    operon = str_replace(operon, "psl", "Psl"),
    operon = str_replace(operon, "curdlan", "Curdlan"),
    operon = str_replace(operon, "diutan", "Diutan"),
    operon = str_replace(operon, "succinoglycan", "Succinoglycan"),
    operon = str_replace(operon, "gellan2", "Gellan2"),
    operon = str_replace(operon, "burkholderia_eps", "Burkholderia EPS"),
    operon = str_replace(operon, "amylovoran", "Amylovoran"),
    operon = str_replace(operon, "ColA", "Colanic Acid"),
    operon = str_replace(operon, "salecan", "Salecan"),
    operon = str_replace(operon, "stewartan", "Stewartan"),
    operon = str_replace(operon, "vps", "Vibrio EPS"),
    operon = str_replace(operon, "rhizobium_eps", "Rhizobium EPS"),
    operon = str_replace(operon, "gellan1", "Gellan1"),
    operon = str_replace(operon, "acetan", "Acetan"),
    operon = str_replace(operon, "s88", "s88"),
    operon = str_replace(operon, "levan", "Levan"),
    operon = str_replace(operon, "methanolan", "Methanolan"),
    operon = str_replace(operon, "synechan", "Synechan")) %>% rename(`Sample#` = id) %>%
   mutate(
     Notes = str_replace(Notes, "Full-scale measurements aerobic stage", "Aerobic stage"),
     Notes = str_replace(Notes, "Full-scale measurements sidestream tank 1", "Sidestream tank 1"),
     Notes = str_replace(Notes, "Full-scale measurements sidestream tank 2", "Sidestream tank 2"),
     Notes = str_replace(Notes, "Full-scale measurements return sludge", "Return sludge"),
     Notes = str_replace(Notes, "Full-scale measurements anoxic stage", "Anoxic stage"))
 
 saveRDS(EPS_table_dds, "/mnt/ahd/EPS_PIPELINE/data/raw/EPS_table_Fullscale_dds.rds")

EPS_table_dds <- readRDS("/mnt/ahd/EPS_PIPELINE/data/raw/EPS_table_Fullscale_dds.rds")


##Plotting data functions

expressionPlot <- function(EPS_table = EPS_table_dds, experiment = "Anoxic stage", legend = "none", n = 5) {
  ggplot(EPS_table %>% filter(Notes == experiment, operon %ni% c("HA (has)", "NulO", "Stewartan")), aes(`Sample#`, counts)) +
    theme_bw() +
    ggtitle(experiment) +
    geom_col(aes(fill = Function)) +
    facet_wrap(~operon) +
    theme(legend.position=legend) +
    ylim(0,150) +
    scale_x_continuous(breaks = 1:n) +
    theme(axis.text.x = element_text(size=6))
}



expressionPlot1 <- function(EPS_table = EPS_table_dds, experiment = "Anoxic stage", legend = "none", n = 5) {
  ggplot(EPS_table %>% filter(Notes == experiment, operon %in% c("HA (has)", "Stewartan")), aes(`Sample#`, counts)) +
    theme_bw() +
    ggtitle(experiment) +
    geom_col(aes(fill = Function)) +
    facet_wrap(~operon) +
    theme(legend.position=legend) +
    ylim(0,1500) +
    scale_x_continuous(breaks = 1:n) +
    theme(axis.text.x = element_text(size=6))
}


##New plot with error bars and collapsing

EPS_table_dds_summarized <- EPS_table_dds %>%
  filter(counts != 0) %>%
  group_by(`Sample#`, Notes, operon) %>%
  summarise(counts = mean(counts)) %>% ungroup() %>%
  rename(`Processing tank` = Notes)
  # group_by(`Notes`, operon) %>% 
  # summarise(mean_counts = mean(counts),
  #           sd_counts = sd(counts))

# EPS_table_dds_summarized$lower <- EPS_table_dds_summarized$mean_counts - EPS_table_dds_summarized$sd_counts
# EPS_table_dds_summarized$upper <- EPS_table_dds_summarized$mean_counts + EPS_table_dds_summarized$sd_counts
# 
EPS_table_dds_summarized$operon <- factor(EPS_table_dds_summarized$operon) %>% fct_reorder(EPS_table_dds_summarized$counts, .desc = TRUE)


summarized_dds_Fullscale <- ggplot(EPS_table_dds_summarized %>% filter(operon %ni% c("NulO")), aes(`Processing tank`, counts, color = `Processing tank`)) +
  geom_beeswarm(cex = 3) +
  facet_wrap(~operon, scales = "free") +
  xlab("") + ylab("Mean normalized count")+
  theme_bw() +
  ylim(c(0,NA)) +
  scale_fill_brewer(palette="Set2") +
  theme(axis.text.x = element_blank()) +
  theme(text=element_text(size=12, colour = "black"),
        axis.text = element_text(size = 6, colour = "black"),
        axis.title = element_text(size = 6, face = "bold", colour = "black"),
        legend.text = element_text(size = 6, colour = "black"),
        legend.title = element_text(size = 6, face ="bold", colour = "black"))


ggsave("/mnt/ahd/EPS_PIPELINE/figures/expression/expression_dds_Fullscale_summarized.png", summarized_dds_Fullscale, limitsize = FALSE, width = 14, height = 10, dpi = 600)

