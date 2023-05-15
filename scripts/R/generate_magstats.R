library("tidyverse")

#--------------------------------------------------------------------------------------------
#  Loading, cleaning, filtering and saving HQ-MAG statistics data from MGP1000 and glomicave |
#--------------------------------------------------------------------------------------------

HQlist <- read_xlsx("data/raw/metaTs_to_692MAGs.xlsx", sheet = 1) %>% select(Genome)

glomistats = read.csv(file = "data/raw/20221102_glomicave_214_HQ_MAGs_w_MiDAS4_tax.txt", header = TRUE, sep = "\t") %>% 
  mutate(name = gsub("bin.", "glom_", name)) %>% filter(name %in% HQlist$Genome) %>% select(name, midas4_tax)

MGPstats = read.csv(file = "data/raw/MiDAS4_gtdb207_MGP1000_MAG_statistics_table.txt", header = TRUE, sep = "\t", skip = 1) %>% 
  rename(name = MAG, midas4_tax = MiDAS4_8_1Tax) %>% select(name, midas4_tax) %>% mutate(name = gsub(".fa", "", name)) %>% filter(name %in% HQlist$Genome)

magstats <- glomistats %>%
  bind_rows(MGPstats) %>% rename(bin = name)
write_tsv(magstats, file="data/raw/magstats.tsv")


