library(DESeq2)
library(tximport)
library(tximportData)
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
# dir <- system.file("extdata", package = "tximportData")
# list.files(dir)
# 
# samples <- list.files("data/metatranscriptomics/RSEM/results/", pattern = "*.genes.results")
# samples
# 
# 
# files <- file.path("./data/metatranscriptomics/RSEM/results", samples)
# files
# names(files) <- paste0("sample ", 1:186)
# files
# files_cut <- files[1:60]
# files_cut
# 
# txi.rsem <- tximport(files_cut, type = "rsem", txIn = FALSE, txOut = FALSE)
# head(txi.rsem$counts)
# 
# saveRDS(txi.rsem, "//wsl$/Ubuntu/mnt/ahd/EPS_PIPELINE/data/metatranscriptomics/txi_rsem.rds")
# 
# txi.rsem <- readRDS("/mnt/ahd/EPS_PIPELINE/data/metatranscriptomics/txi_rsem.rds")
# 
#  
# #Filtering out genes with no counts in any sample
# txi.rsem$abundance <- txi.rsem$abundance[apply(txi.rsem$length,
#                                 1,
#                                 function(row) all(row !=0 )),]
# 
# txi.rsem$counts <- txi.rsem$counts[apply(txi.rsem$length,
#                              1,
#                              function(row) all(row !=0 )),]
# 
# txi.rsem$length <- txi.rsem$length[apply(txi.rsem$length,
#                              1,
#                              function(row) all(row !=0 )),]
# 
# 
# #Joining RSEM data with EPS data, significantly reducing the number of genes and thus speeding up the analysis
# txi.rsem$abundance <- txi.rsem$abundance[rownames(txi.rsem$abundance) %in% database_stack$gene_id,]
# 
# txi.rsem$counts <- txi.rsem$counts[rownames(txi.rsem$counts) %in% database_stack$gene_id,]
# 
# txi.rsem$length <- txi.rsem$length[row.names(txi.rsem$length) %in% database_stack$gene_id,]
# 
# saveRDS(txi.rsem, "/user_data/ahd/EPS_PIPELINE/data/metatranscriptomics/txi_rsem.rds")

txi.rsem <- readRDS("/user_data/ahd/EPS_PIPELINE/data/metatranscriptomics/txi_rsem_PAO_EPS.rds")

#Loading metadata
metadata <- read_xlsx("user_data/ahd/EPS_PIPELINE/data/metatranscriptomics/metadata_R.xlsx")

metadata <- metadata[1:60,]

rownames(metadata) <- colnames(txi.rsem$counts)

#Importing tximport object into DESeq2
dds <- DESeqDataSetFromTximport(txi.rsem, metadata, ~Time_point + Substrate_added)
dds$Substrate_added <- relevel( dds$Substrate_added, "Control" )
# dds is now ready for DESeq() see DESeq2 vignette

#Performing DESq2	command (computationaly expensive, unless you've joined the RSEM data with your own gene annotation database)
dds <- DESeq(dds)

saveRDS(dds, "/user_data/ahd/EPS_PIPELINE/data/metatranscriptomics/DESeq2_EPS_PAO.rds")

#Now you can load the dds object and perform further analysis
dds <- readRDS("/mnt/ahd/data/metatranscriptomics/DESeq2_EPS_PAO.rds")

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

