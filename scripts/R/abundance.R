library(dplyr)
library(ggplot2)
library(SummarizedExperiment)
library(ggpmisc)

abundance_metaT <- read_tsv("/mnt/ahd/EPS_PIPELINE/data/metatranscriptomics/relative_abundance_metaTs.txt") %>%
  select(1, `HQ_sp_692_concat.fa.LIB-Glomicave-0151_resequenced_mRNA_R1.fq Relative Abundance (%)`:`HQ_sp_692_concat.fa.LIB-Glomicave-0189_mRNA_R1.fq Relative Abundance (%)`) %>%
  filter(Genome != "unmapped") %>%
  dplyr::rename(
    `SA-Glomicave-0151` = 2,
    `SA-Glomicave-0152` = 3,
    `SA-Glomicave-0153` = 4,
    `SA-Glomicave-0154` = 5,
    `SA-Glomicave-0155` = 6,
    `SA-Glomicave-0156` = 7,
    `SA-Glomicave-0157` = 8,
    `SA-Glomicave-0158` = 9,
    `SA-Glomicave-0159` = 10,
    `SA-Glomicave-0160` = 11,
    `SA-Glomicave-0161` = 12,
    `SA-Glomicave-0162` = 13,
    `SA-Glomicave-0163` = 14,
    `SA-Glomicave-0164` = 15,
    `SA-Glomicave-0165` = 16,
    `SA-Glomicave-0166` = 17,
    `SA-Glomicave-0167` = 18,
    `SA-Glomicave-0168` = 19,
    `SA-Glomicave-0169` = 20,
    `SA-Glomicave-0170` = 21,
    `SA-Glomicave-0171` = 22,
    `SA-Glomicave-0172` = 23,
    `SA-Glomicave-0173` = 24,
    `SA-Glomicave-0174` = 25,
    `SA-Glomicave-0175` = 26,
    `SA-Glomicave-0176` = 27,
    `SA-Glomicave-0177` = 28,
    `SA-Glomicave-0178` = 29,
    `SA-Glomicave-0179` = 30,
    `SA-Glomicave-0180` = 31,
    `SA-Glomicave-0181` = 32,
    `SA-Glomicave-0182` = 33,
    `SA-Glomicave-0183` = 34,
    `SA-Glomicave-0184` = 35,
    `SA-Glomicave-0185` = 36,
    `SA-Glomicave-0186` = 37,
    `SA-Glomicave-0187` = 38,
    `SA-Glomicave-0188` = 39,
    `SA-Glomicave-0189` = 40,
    MAG_id = 1) %>%
pivot_longer(cols = 2:40, names_to = "Sample_ID", values_to = "Relative metatranscriptomic abundance (%)") %>% 
  mutate(`Relative metatranscriptomic abundance (%)` = ifelse(`Relative metatranscriptomic abundance (%)` == 0, 0.0000003, `Relative metatranscriptomic abundance (%)`))

saveRDS(abundance_metaT, file = "/mnt/ahd/EPS_PIPELINE/data/metatranscriptomics/abundance_metaT.rds")

abundance_metaG <- read_tsv("/user_data/ahd/EPS_PIPELINE/data/metatranscriptomics/relative_abundance_metagenome.tsv") %>% 
  dplyr::rename(MAG_id = Genome, `Relative abundance (%)` = 2) %>% filter(MAG_id != "unmapped")

saveRDS(abundance_metaG, file = "/user_data/ahd/EPS_PIPELINE/data/metatranscriptomics/abundance_metaG.rds")

phylum <- read_tsv("./data/raw/magstats.tsv") %>% 
  mutate(
    phylum =   str_extract(midas4_tax, "p:[^,]*"),
    phylum =   str_remove(phylum, "p:"),
    
    class =  str_extract(midas4_tax, "c:[^,]*"),
    class =  str_remove(class, "c:"),
    class =  str_remove(class, ";"),
    
    order =  str_extract(midas4_tax, "o:[^,]*"),
    order =  str_remove(order, "o:"),
    order =  str_remove(order, ";"),
    
    family =  str_extract(midas4_tax, "f:[^,]*"),
    family =  str_remove(family, "f:"),
    family =  str_remove(family, ";"),
    
    genus =  str_extract(midas4_tax, "g:[^,]*"),
    genus =  str_remove(genus, "g:"),
    genus =  str_remove(genus, ";"),
    
    species = str_extract(midas4_tax, "g:[^;]*"),
    # species =  str_remove(genus, "g:"),
    # species = str_remove(species, "s:"),
    # species = str_remove(species, ";")
    # 
  ) %>% 
  dplyr::rename(ID = bin) %>%
  select(ID, phylum, class, order, family, genus, species) %>% 
  setNames(c("MAG_id", "phylum", "class", "order", "family", "genus", "species"))



###
# Comparing all taxa with a heatmap
###

abundance <- read_tsv("/mnt/ahd/EPS_PIPELINE/data/metatranscriptomics/relative_abundance_metagenome.tsv") %>% 
  dplyr::rename(MAG_id = Genome, Metagenome = 2) %>% filter(MAG_id != "unmapped") %>% left_join(abundance_metaT) %>% left_join(phylum) %>%
  pivot_longer(cols = 2:41, values_to = "Relative Abundance (%)", names_to = "Sample_ID") %>%
  mutate(Sample_ID = ifelse(Sample_ID != "Metagenome", substr(Sample_ID, nchar(Sample_ID)-2, nchar(Sample_ID)), Sample_ID)) %>%
  mutate(Sample_ID = ifelse(Sample_ID != "Metagenome", paste0("Sample ", Sample_ID), Sample_ID))


abundance$class <- factor(abundance$class) %>% fct_reorder(abundance$`Relative Abundance (%)`, .desc = TRUE)
abundance$order <- factor(abundance$order) %>% fct_reorder(abundance$`Relative Abundance (%)`, .desc = TRUE)
abundance$genus <- factor(abundance$genus) %>% fct_reorder(abundance$`Relative Abundance (%)`, .desc = TRUE)
abundance$species <- factor(abundance$species) %>% fct_reorder(abundance$`Relative Abundance (%)`, .desc = TRUE)


abundances <- ggplot(abundance, aes(Sample_ID, species)) +
  geom_tile(aes(fill =`Relative Abundance (%)`)) +
  scale_fill_gradientn(trans = "log", colors = c("black", "orange", "red"), 
                       na.value = "black", breaks = c(0.00001, 0.001, 0.05, 3), name = "Log(Relative Abundance (%))") +
  theme(axis.text.y = element_text(size = 1),
        # axis.ticks.y = element_blank(),
        # axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_text(size = 10),
        axis.text.x = element_text(size = 3)
  ) +
  ylab("Species")

ggsave("/mnt/ahd/EPS_PIPELINE/figures/abundances/abundance_overview_species.png", abundances, limitsize = FALSE, width = 14, height = 10, dpi = 1000)

###
# Comparing Metagenome with the rest using scatterplots
###

abundance2 <- read_tsv("/mnt/ahd/EPS_PIPELINE/data/metatranscriptomics/relative_abundance_metagenome.tsv") %>% 
  dplyr::rename(MAG_id = Genome, Metagenome = 2) %>% filter(MAG_id != "unmapped") %>% left_join(abundance_metaT) %>%
  pivot_longer(cols = 3:41, values_to = "Relative Abundance (%)", names_to = "Sample_ID")

correlation_seperate <-ggplot(abundance2, aes(Metagenome, `Relative Abundance (%)`)) +
  geom_point() +
  facet_wrap(~Sample_ID) +
  scale_x_continuous(trans = "log", breaks = c(0.000001, 0.001, 3)) +
  scale_y_continuous(trans = "log", breaks = c(0.000001, 0.001, 3)) +
  ylab("Log(Relative abundance (%) - Samples)") +
  xlab("Log(Relative abundance (%) - Metagenome)") +
  stat_poly_line() +
  stat_poly_eq(aes(label = paste(after_stat(rr.label), sep = "*\", \"*"))) +
  geom_point()

abundance3 <- abundance2 %>% dplyr::group_by(MAG_id, Metagenome) %>% dplyr::summarise(`Relative Abundance (%)` = mean(`Relative Abundance (%)`))

correlation_together <- ggplot(abundance3, aes(Metagenome, `Relative Abundance (%)`)) +
  geom_point() +
  scale_x_continuous(trans = "log", breaks = c(0.000001, 0.001, 3)) +
  scale_y_continuous(trans = "log", breaks = c(0.000001, 0.001, 3)) +
  ylab("Log(Relative abundance (%)) - Mean of Samples") +
  xlab("Log(Relative abundance (%)) - Metagenome") +
  stat_poly_line() +
  stat_poly_eq(aes(label = paste(after_stat(eq.label),
                                 after_stat(rr.label), sep = "*\", \"*"))) +
  geom_point()

ggsave("/mnt/ahd/EPS_PIPELINE/figures/abundances/rel_abun_correlation_seperate_log.png", correlation_seperate, limitsize = FALSE, width = 14, height = 10, dpi = 600)

###
# Comparing Samples with each other
###

# Comparing Sample 151 with the rest

abundance4 <- read_tsv("/mnt/ahd/EPS_PIPELINE/data/metatranscriptomics/relative_abundance_metagenome.tsv") %>% 
  dplyr::rename(MAG_id = Genome, Metagenome = 2) %>% filter(MAG_id != "unmapped") %>% left_join(abundance_metaT) %>%
  select(-Metagenome) %>% pivot_longer(cols = 3:40, values_to = "Relative Abundance (%)", names_to = "Sample_ID")

correlation_seperate <-ggplot(abundance4, aes(`SA-Glomicave-0151`, `Relative Abundance (%)`)) +
  geom_point() +
  facet_wrap(~Sample_ID) +
  scale_x_continuous(trans = "identity", breaks = c(0.000001, 0.001, 3)) +
  scale_y_continuous(trans = "identity", breaks = c(0.000001, 0.001, 3)) +
  ylab("Relative abundance (%) - Samples") +
  xlab("Relative abundance (%) - Sample 151") +
  stat_poly_line() +
  stat_poly_eq(aes(label = paste(after_stat(rr.label), sep = "*\", \"*"))) +
  geom_point()

ggsave("/mnt/ahd/EPS_PIPELINE/figures/abundances/rel_abun_correlation_seperate_0151.png", correlation_seperate, limitsize = FALSE, width = 14, height = 10, dpi = 600)



#Comparing Sample 175 with the rest

abundance5 <- read_tsv("/mnt/ahd/EPS_PIPELINE/data/metatranscriptomics/relative_abundance_metagenome.tsv") %>% 
  dplyr::rename(MAG_id = Genome, Metagenome = 2) %>% filter(MAG_id != "unmapped") %>% left_join(abundance_metaT) %>%
  select(-Metagenome) %>% pivot_longer(cols = c(2:25, 27:40), values_to = "Relative Abundance (%)", names_to = "Sample_ID")

correlation_seperate <-ggplot(abundance5, aes(`SA-Glomicave-0175`, `Relative Abundance (%)`)) +
  geom_point() +
  facet_wrap(~Sample_ID) +
  scale_x_continuous(trans = "identity", breaks = c(0.000001, 0.001, 3)) +
  scale_y_continuous(trans = "identity", breaks = c(0.000001, 0.001, 3)) +
  ylab("Relative abundance (%) - Samples") +
  xlab("Relative abundance (%) - Sample 175") +
  stat_poly_line() +
  stat_poly_eq(aes(label = paste(after_stat(rr.label), sep = "*\", \"*"))) +
  geom_point()

ggsave("/mnt/ahd/EPS_PIPELINE/figures/abundances/rel_abun_correlation_seperate_0175.png", correlation_seperate, limitsize = FALSE, width = 14, height = 10, dpi = 600)



