library(tidyverse)
library(ggplot2)

midas4tax <- read.csv2("./data/raw/magstats.tsv", sep = "\t") %>%
  separate(midas4_tax, sep = ",", into = c("domain", "phylum","class","order","family","genus","species"))

genus <- ggplot(midas4tax %>% count(phylum), aes(y=reorder(phylum, n), x="", fill=n))+
  geom_tile(show.legend = FALSE)+
  scale_fill_gradientn(colors = c("white", "orange", "red"), limits = c(0,250), name = "Rel. abund.")+
  geom_text(aes(label = round(n, 1))) +
  theme_minimal()
genus
  


midas41000 <- read.csv2("//wsl$/Ubuntu/mnt/server/user_data/ahd/EPS_PIPELINE/data/raw/MiDAS4_gtdb207_MGP1000_MAG_statistics_table.txt", sep = "\t", header = T)%>% 
  separate(X.18, sep = ";", into = c("domain", "phylum","class","order","family","genus","species")) 