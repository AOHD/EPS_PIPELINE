library("tidyverse")

#-----------------------------------------------------------------
#        Defining the gffRead function for loading gff data       
#-----------------------------------------------------------------
gffRead <- function(gffFile){
  gff <- read.delim(gffFile, header=F, comment.char="#")
  colnames(gff) = c("seqname", "source", "feature", "start", "end",
                    "score", "strand", "frame", "attributes")
  gff <- gff[complete.cases(gff), ]
  
  gff <- gff %>%
    separate(attributes, c("prokkaID", "attributes"), sep = ";", extra = "merge") %>%
    separate(prokkaID, c("ID1", "ProkkaNO"), sep = -5)
  
  
  return(gff)}

#-----------------------------------------------------------------
#          Using the gffRead function to import gff files         
#-----------------------------------------------------------------
gff <- tibble(
  file = list.files("data/raw/reduced_gff"),
  location = list.files("data/raw/reduced_gff", full.names = TRUE),
  ID = str_remove(file, ".gff"),
  df = map(location, gffRead)
  ) %>% 
  unnest(df)

#-----------------------------------------------------------------
#  Removing columns redundant to psiblast data in gff dataframe   
#-----------------------------------------------------------------
gff <-  gff %>% select(!c("score", "feature", "file", "location", 
                          "frame", "attributes",  "source", "ID1"))

write_tsv(gff, file="data/raw/gff.tsv")


