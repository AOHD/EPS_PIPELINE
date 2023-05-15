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

#Importing RSEM data with tximport

samples <- list.files("user_data/ahd/EPS_PIPELINE/data/metatranscriptomics/RSEM/results/", pattern = "*.genes.results")
samples


files <- file.path("user_data/ahd/EPS_PIPELINE/data/metatranscriptomics/RSEM/results", samples)
files
names(files) <- paste0("sample ", 1:186)
files
files_cut <- files[61:110]
files_cut

txi.rsem <- tximport(files_cut, type = "rsem", txIn = FALSE, txOut = FALSE)
head(txi.rsem$counts)


#Filtering out genes with no counts in any sample
txi.rsem$abundance <- txi.rsem$abundance[apply(txi.rsem$length,
                                1,
                                function(row) all(row !=0 )),]

txi.rsem$counts <- txi.rsem$counts[apply(txi.rsem$length,
                             1,
                             function(row) all(row !=0 )),]

txi.rsem$length <- txi.rsem$length[apply(txi.rsem$length,
                             1,
                             function(row) all(row !=0 )),]


#Importing EPS gene annotation database
database_stack <- data.frame(
  Target_label = as.character(),
  Psiblast = as.character()
)

for (file in list.files("user_data/ahd/EPS_PIPELINE/output_proximity_filtration/psi_operon_full")) {
  temp <- read.csv2(paste0("user_data/ahd/EPS_PIPELINE/output_proximity_filtration/psi_operon_full/",file), sep = "\t") %>% filter(!is.na(Query_label)) %>% select(Psiblast, Target_label, Query_label, Function)
  
  database_stack <- database_stack %>% bind_rows(temp)
  
}

database_stack <- database_stack %>% filter(!is.na(Psiblast)) %>% rename(gene_id = Target_label, operon = Psiblast, gene_annotation = Query_label)


#Joining RSEM data with EPS data, significantly reducing the number of genes and thus speeding up the analysis
txi.rsem$abundance <- txi.rsem$abundance[rownames(txi.rsem$abundance) %in% database_stack$gene_id,]

txi.rsem$counts <- txi.rsem$counts[rownames(txi.rsem$counts) %in% database_stack$gene_id,]

txi.rsem$length <- txi.rsem$length[row.names(txi.rsem$length) %in% database_stack$gene_id,]

txi.rsem <- readRDS("/user_data/ahd/EPS_PIPELINE/data/metatranscriptomics/txi_rsem_DPAO_EPS.rds")


#Loading metadata
metadata <- read_xlsx("/user_data/ahd/EPS_PIPELINE/data/metatranscriptomics/metadata_R.xlsx")

metadata <- metadata[61:110,]

rownames(metadata) <- colnames(txi.rsem$counts)

#Importing tximport object into DESeq2
dds <- DESeqDataSetFromTximport(txi.rsem, metadata, ~Time_point + N_added)
dds$Substrate_added <- relevel( dds$N_added, "NA" )
# dds is now ready for DESeq() see DESeq2 vignette



#Performing DESq2	command (computationaly expensive, unless you've joined the RSEM data with your own gene annotation database)
dds <- DESeq(dds)

saveRDS(dds, "/user_data/ahd/EPS_PIPELINE/data/metatranscriptomics/DESeq2_EPS_DPAO.rds")

#Now you can load the dds object and perform further analysis
dds <- readRDS("/user_data/ahd/data/metatranscriptomics/DESeq2_EPS.rds")

res <- results(dds, contrast = c("N_added", "NA", "NO3"))

summary(res)

# Transform count data using the variance stablilizing transform
deseq2VST <-  varianceStabilizingTransformation(dds)

# Convert the DESeq transformed object to a data frame
deseq2VST <- assay(deseq2VST)
deseq2VST <- as.data.frame(deseq2VST)
deseq2VST$Gene <- rownames(deseq2VST)
head(deseq2VST)

# Keep only the significantly differentiated genes where the fold-change was at least 3
sigGenes <- rownames(deseq2ResDF[deseq2ResDF$padj <= .05,])
deseq2VST <- deseq2VST[deseq2VST$Gene %in% sigGenes,]



#Exploring results


res <- results(dds, contrast = c("Substrate_added", "Control", "Arginine and aspartate"))

deseq2ResDF <- as.data.frame(res) %>% mutate(gene_id = rownames(deseq2ResDF)) %>% left_join(database_stack, by = "gene_id")

