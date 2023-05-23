library(tidyverse)
library(treeio)
library(ggtree)
library(ggtreeExtra)
library(ggplot2)
library(ggnewscale)
library(ggrepel)
library(scales)
library(readxl)
library(viridis)
library(patchwork)
library(data.table)
library(aplot)
library(ggplotify)
`%ni%` <- Negate(`%in%`)


##################
# Full-scale 
##################

# source("./scripts/R/HQ_MAG_tree.R")
# saveRDS(tree, file = "/mnt/ahd/EPS_PIPELINE/data/metatranscriptomics/HQtree.rds")

tree <- readRDS("/mnt/ahd/EPS_PIPELINE/data/raw/HQtree.rds")
EPS_table_RSEM <- readRDS("/mnt/ahd/EPS_PIPELINE/data/metatranscriptomics/EPS_table_Fullscale_RSEM.rds") %>% 
  mutate(`Processing tank` = ifelse(`Processing tank` == "Anoxic stage", "Anoxic", `Processing tank`),
         `Processing tank` = ifelse(`Processing tank` == "Aerobic stage", "Aerobic", `Processing tank`),
         `Processing tank` = ifelse(`Processing tank` == "Return sludge", "RS", `Processing tank`),
         `Processing tank` = ifelse(`Processing tank` == "Sidestream tank 1", "ST1", `Processing tank`),
         `Processing tank` = ifelse(`Processing tank` == "Sidestream tank 2", "ST2", `Processing tank`))

EPS_table_RSEM$`Processing tank` <- factor(EPS_table_RSEM$`Processing tank`,levels = 
                                                        c("Anoxic", "Aerobic", "RS", "ST1", "ST2"))



##################
# Full-scale, several experiments
##################
expression_phylogeny_plural <- function(eps, legend = "right", limits = c(1, NA),
                                 treelayout = "rectangular", trans = "identity", 
                                 width = 2, offset = 2, title = eps) {

  # EPS_table_filt <- EPS_table_RSEM %>% filter(operon %in% c(eps)) %>% mutate(label = substring(gene_id, 1, nchar(gene_id) - 6)) %>% 
  #   group_by(label,`Processing tank`) %>% summarise(TPM = mean(TPM)) %>% filter(any(TPM >= limits[1])) %>% ungroup() %>% 
  #   pivot_wider(names_from = `Processing tank`, values_from = TPM) %>% as.data.frame()
  # 
  
  EPS_table_filt <- EPS_table_RSEM %>% filter(operon %in% c(eps)) %>% mutate(label = substring(gene_id, 1, nchar(gene_id) - 6)) %>% 
    group_by(`Relative abundance (%)`, label, `Sample#`,`Processing tank`) %>% summarise(TPM = sum(TPM)) %>% filter(any(TPM >= limits[1])) %>% ungroup() %>% 
    group_by(`Relative abundance (%)`, label, `Processing tank`) %>% summarise(TPM = mean(TPM)) %>%
    pivot_wider(names_from = `Processing tank`, values_from = TPM) %>% as.data.frame()  
  
  
  rownames(EPS_table_filt) <- EPS_table_filt$label 
  
  EPS_table_filt <- EPS_table_filt %>% select(-label)
  
  tree_tib_relevant <- tree %>% 
    as_tibble() %>% filter(label %in% rownames(EPS_table_filt))
  
  tree_tib_irrelevant <- tree %>% 
    as_tibble() %>% filter(label %ni% tree_tib_relevant$label)
  
  tree_filt <- tree %>% drop.tip(tree_tib_irrelevant$label)
  tree_filt <- tree_filt %>% as_tibble() %>%
    mutate(
      genus = str_replace_all(genus, 
                              pattern = "Ca_", 
                              replacement = "\\*Ca.\\* ") ,
      species = str_replace_all(species,
                                pattern = "Ca_*_",
                                replacement = ""),
      species = str_replace_all(species,
                                pattern = "Microthrix_",
                                replacement = ""),
      species = str_replace_all(species,
                                pattern = "Nitrospira_",
                                replacement = ""),
      species = str_replace_all(species,
                                pattern = "Sphingopyxis_",
                                replacement = ""),
      species = str_replace_all(species,
                                pattern = "Phreatobacter_",
                                replacement = ""),
      species = str_replace_all(species,
                                pattern = "Haliscomenobacter_",
                                replacement = ""),
      species = str_replace_all(species,
                                pattern = "Tabrizicola_",
                                replacement = ""),
      species = str_replace_all(species,
                                pattern = "Dechloromonas_",
                                replacement = ""),
      species = str_replace_all(species,
                                pattern = "Amarolinea_",
                                replacement = ""),
      species = str_replace_all(species,
                                pattern = "Methylophosphatis_",
                                replacement = ""),
      species = str_replace_all(species,
                                pattern = "Accumulibacter_",
                                replacement = ""),
      genus = ifelse(label %in% c("AalE_18-Q3-R2-46_BAT3C.188", "AalW_18-Q3-R10-53_BAT3C.524", "Bjer_18-Q3-R1-45_BAT3C.93",
                                  "Ega_18-Q3-R5-49_MAXAC.001", "EsbW_18-Q3-R4-48_BAT3C.295", "Hirt_18-Q3-R61-65_BAT3C.386",
                                  "Hjor_18-Q3-R7-51_BAT3C.81_sub", "Hjor_18-Q3-R7-51_MAXAC.088", "Mari_18-Q3-R65-66_BAT3C.41",
                                  "Ribe_18-Q3-R11-54_MAXAC.001"), "*Ca.* Phosphoribacter", genus),
      species = ifelse(label %in% c("EsbW_18-Q3-R4-48_BAT3C.295", "AalW_18-Q3-R10-53_BAT3C.524", "Bjer_18-Q3-R1-45_BAT3C.93"), "baldrii", species),
      species = ifelse(label %in% c("AalE_18-Q3-R2-46_BAT3C.188", "Ega_18-Q3-R5-49_MAXAC.001", "Ribe_18-Q3-R11-54_MAXAC.001"), "hodrii", species),
      species = ifelse(label %in% c("Hirt_18-Q3-R61-65_BAT3C.386"), "Pbr3", species),
      species = ifelse(label %in% c("Hjor_18-Q3-R7-51_MAXAC.088"), "Pbr4", species),
      species = ifelse(label %in% c("Hjor_18-Q3-R7-51_BAT3C.81_sub"), "Pbr5", species),
      species = ifelse(label %in% c("Mari_18-Q3-R65-66_BAT3C.41"), "Pbr6", species),
      genus = ifelse(species == "midas_s_45", "*Ca.* Lutibacillus", genus),
      species = ifelse(species == "midas_s_45", "vidarii", species),
      genus = ifelse(str_detect(genus, "Ca.|midas"), genus, paste0("*",genus,"*")),
      species = ifelse(str_detect(species, "midas") | str_detect(genus, "Ca."), species, paste0("*",species,"*"))
    ) %>%
    as.treedata()

  rel_abun_tree <- EPS_table_RSEM %>% 
    mutate(label = substring(gene_id, 1, nchar(gene_id) - 6)) %>% ungroup() %>%
    select(label, `Relative abundance (%)`) %>% filter(label %in% tree_tib_relevant$label) %>% group_by(label, `Relative abundance (%)`) %>%
    distinct(label)
  
  
  tree_filt_tib <- tree_filt %>% as_tibble() %>% left_join(rel_abun_tree, by = "label") %>% select(-`Relative abundance (%)`) %>%
    rename(label_original = label) 
  
  tree_filt_tib$species <- with(tree_filt_tib, make.unique(as.character(species)))
  
  tree_filt_tib <- tree_filt_tib %>% mutate(label = paste(genus, species)) %>% as.treedata()
  
  tree_plot <- ggtree(tree_filt_tib, layout = treelayout, branch.length = "none") +
    # geom_tiplab(
    #   mapping = aes(node = node, label = paste(label)),
    #   geom = "text",
    #   size = 3
    # ) +
    xlim_tree(1) +
    ggtitle(title) +
    theme(
      plot.title = element_markdown(size = 40)
    )
  
  rel_abun_tree <- rel_abun_tree %>% select(-`Relative abundance (%)`) %>% 
    rename(label_original = label) %>% left_join(as_tibble(tree_filt_tib), by = "label_original")
  
  
  EPS_table_filt_plot <- EPS_table_RSEM %>% filter(operon %in% c(eps)) %>% 
    mutate(label_original = substring(gene_id, 1, nchar(gene_id) - 6)) %>% left_join(as_tibble(tree_filt_tib), by = "label_original") %>%
    group_by(`Relative abundance (%)`, label, `Sample#`,`Processing tank`) %>% summarise(TPM = sum(TPM)) %>% filter(any(TPM >= limits[1])) %>% ungroup() %>% 
    group_by(`Relative abundance (%)`, label, `Processing tank`) %>% summarise(TPM = mean(TPM)) %>% ungroup() %>% filter(!is.na(label))
  
  
  
  abun_plot <- ggplot(rel_abun_tree, aes(x = "Rel. abund.", y = label)) +
    geom_point(aes(size = `Relative abundance (%)`)) + 
    theme_bw() +
    scale_size_continuous(range = c(1, 15)) +
    theme(axis.line = element_blank(),
          panel.border = element_blank(),
          axis.text.y = element_markdown(colour = "Black", size = 30),
          axis.ticks.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          plot.title = element_text(size = 35),
          axis.text.x = element_blank(),
          legend.text = element_text(size = 30),
          legend.title = element_markdown(size = 35)
    )
  
  EPS_plot <- ggplot(EPS_table_filt_plot, aes(`Processing tank`, label, size = TPM)) +
    geom_point(color = "orange")+
    theme_bw() +
    scale_size_continuous(range = c(1, 15)) +
    theme(axis.line = element_blank(),
          panel.border = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          plot.title = element_text(size = 35),
          axis.text.x = element_text(size = 30, color = "Black", angle = 30),
          legend.text = element_text(size = 30),
          legend.title = element_markdown(size = 35)
          )
  
  MAG_TPM_plot <- ggplot(EPS_table_filt_plot, aes(`Processing tank`, label, size = TPM/`Relative abundance (%)`)) +
    geom_point(color = "red")+
    theme_bw() +
    scale_size_continuous(range = c(1, 20)) +
    theme(axis.line = element_blank(),
          panel.border = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          plot.title = element_text(size = 35),
          axis.text.x = element_text(size = 30, color = "Black", angle = 30),
          legend.text = element_text(size = 30),
          legend.title = element_text(size = 35)
    ) +
    guides(size=guide_legend(title="MAG-normalized TPM"))
  
  heattree <- EPS_plot %>% insert_left(abun_plot, width = 0.2) %>% insert_left(tree_plot, width = 1) %>% insert_right(MAG_TPM_plot)
  
  return(heattree)
  
}

HA <- expression_phylogeny_plural(eps = "HA (has)", title = "HA (has)", offset = 1)

ggsave("/mnt/ahd/EPS_PIPELINE/figures/expression_connection/HA_has_new.pdf", HA, limitsize = FALSE, width = 30, height = 18, dpi = 300)


Glucorhamnan <- expression_phylogeny_plural(eps = "Glucorhamnan", title = "Glucorhamnan", offset = 1)

ggsave("/mnt/ahd/EPS_PIPELINE/figures/expression_connection/Glucorhamnan_new.pdf", Glucorhamnan, limitsize = FALSE, width = 30, height = 18, dpi = 300)


Stewartan <- expression_phylogeny_plural(eps = "Stewartan", title = "Stewartan", offset = 1)

ggsave("/mnt/ahd/EPS_PIPELINE/figures/expression_connection/Stewartan_new.pdf", Stewartan, limitsize = FALSE, width = 30, height = 18, dpi = 300)


celluloseNA <- expression_phylogeny_plural(eps = "Unclassified cellulose", title = "Unclassified\ncellulose", offset = 2)

ggsave("/mnt/ahd/EPS_PIPELINE/figures/expression_connection/cellulose_NA_new.pdf", celluloseNA, limitsize = FALSE, width = 30, height = 18, dpi = 300)


Levan <- expression_phylogeny_plural(eps = "Levan", title = "Levan", offset = 1)

ggsave("/mnt/ahd/EPS_PIPELINE/figures/expression_connection/Levan_new.pdf", Levan, limitsize = FALSE, width = 30, height = 18, dpi = 300)


Alginate <- expression_phylogeny_plural(eps = "Alginate", title = "Alginate", offset = 0.01, limits = c(0.1, NA))

ggsave("/mnt/ahd/EPS_PIPELINE/figures/expression_connection/Alginate_new.pdf", Alginate, limitsize = FALSE, width = 30, height = 18, dpi = 300)


E_faecalis_PS <- expression_phylogeny_plural(eps = "*E. faecalis* PS", title = "*E. faecalis* PS")

ggsave("/mnt/ahd/EPS_PIPELINE/figures/expression_connection/E_faecalis_PS_new.pdf", E_faecalis_PS, limitsize = FALSE, width = 30, height = 18, dpi = 300)


pmHAS <- expression_phylogeny_plural(eps = "HA (pmHAS)", title = "HA (pmHAS)")

ggsave("/mnt/ahd/EPS_PIPELINE/figures/expression_connection/pmHAS_new.pdf", pmHAS, limitsize = FALSE, width = 30, height = 18, dpi = 300)

Gellan <- expression_phylogeny_plural(eps = "Gellan 2", title = "Gellan")

ggsave("/mnt/ahd/EPS_PIPELINE/figures/expression_connection/gellan_new.pdf", Gellan, limitsize = FALSE, width = 30, height = 18, dpi = 300)

Galactoglucan <- expression_phylogeny_plural(eps = "Galactoglucan", title = "Galactoglucan")

ggsave("/mnt/ahd/EPS_PIPELINE/figures/expression_connection/Galactoglucan_new.pdf", Galactoglucan, limitsize = FALSE, width = 30, height = 18, dpi = 300)

cellulose_Ace <- expression_phylogeny_plural(eps = "Acetylated cellulose", title = "Acetylated cellulose", offset = 0.01)

ggsave("/mnt/ahd/EPS_PIPELINE/figures/expression_connection/cellulose_Ace_new.pdf", cellulose_Ace, limitsize = FALSE, width = 30, height = 18, dpi = 300)

Diutan <- expression_phylogeny_plural(eps = "Diutan", title = "Diutan", offset = 0.01)

ggsave("/mnt/ahd/EPS_PIPELINE/figures/expression_connection/Diutan_new.pdf", Diutan, limitsize = FALSE, width = 30, height = 18, dpi = 300)

Pel <- expression_phylogeny_plural(eps = "Pel", title = "Pel", offset = 0.01)

ggsave("/mnt/ahd/EPS_PIPELINE/figures/expression_connection/Pel_new.pdf", Pel, limitsize = FALSE, width = 30, height = 18, dpi = 300)

CelluloseII <- expression_phylogeny_plural(eps = "Cellulose II", title = "Cellulose II", offset = 0.01)

ggsave("/mnt/ahd/EPS_PIPELINE/figures/expression_connection/CelluloseII_new.pdf", CelluloseII, limitsize = FALSE, width = 30, height = 18, dpi = 300)

S88 <- expression_phylogeny_plural(eps = "S88", title = "S88", offset = 0.01)

ggsave("/mnt/ahd/EPS_PIPELINE/figures/expression_connection/S88_new.pdf", S88, limitsize = FALSE, width = 30, height = 18, dpi = 300)

# CelluloseIII <- expression_phylogeny_plural(eps = "Cellulose III", title = "Cellulose III", offset = 0.01)
# 
# ggsave("/mnt/ahd/EPS_PIPELINE/figures/expression_connection/CelluloseIII_new.pdf", CelluloseIII, limitsize = FALSE, width = 30, height = 18, dpi = 300)

Sphingans <- expression_phylogeny_plural(eps = c("S88", "Diutan", "Gellan 2"), title = "Sphingans\n(S88, Diutan, Gellan)", offset = 0.01)

ggsave("/mnt/ahd/EPS_PIPELINE/figures/expression_connection/Sphingans_new.pdf", Sphingans, limitsize = FALSE, width = 30, height = 18, dpi = 300)




expression_phylogeny <- function(eps, experiment = "Anoxic stage", legend = "right", 
                                 treelayout = "rectangular", trans = "identity", width = 2, offset = 0.5, lowlim = 1) {
  
  EPS_table_filt <- EPS_table_RSEM %>% filter(operon %in% c(eps), `Processing tank` == experiment) %>% mutate(label = substring(gene_id, 1, nchar(gene_id) - 6)) %>% 
    group_by(label, gene_annotation) %>% summarise(TPM = mean(TPM)) %>% filter(any(TPM >= lowlim)) %>% ungroup() %>% 
    pivot_wider(names_from = gene_annotation, values_from = TPM) %>% as.data.frame()
  
  rownames(EPS_table_filt) <- EPS_table_filt$label 
  
  EPS_table_filt <- EPS_table_filt %>% select(-label)
  
  tree_tib_relevant <- tree %>% 
    as_tibble() %>% filter(label %in% rownames(EPS_table_filt))
  
  tree_tib_irrelevant <- tree %>% 
    as_tibble() %>% filter(label %ni% tree_tib_relevant$label)
  
  tree_filt <- tree %>% drop.tip(tree_tib_irrelevant$label)
  
  rel_abun <- EPS_table_RSEM %>% 
    mutate(label = substring(gene_id, 1, nchar(gene_id) - 6)) %>% ungroup() %>%
    select(label, `Relative abundance (%)`) %>% filter(label %in% tree_tib_relevant$label) %>% group_by(label, `Relative abundance (%)`) %>%
    distinct(label)
  
  tree_filt_tib <- tree_filt %>% as_tibble() %>% left_join(rel_abun, by = "label") %>% as.treedata()
  
  
  tree_plot <- ggtree(tree_filt_tib, layout = treelayout) +
    geom_tiplab( 
      mapping = aes(node = node, label = paste(label, genus)),
      geom = "text",
      size = 3
    ) + 
    geom_tippoint(aes(size = `Relative abundance (%)`))
  
  heattree <-gheatmap(tree_plot,EPS_table_filt, offset = offset, colnames_angle = -45, colnames_offset_y = 0.4, width = width, hjust = 0) +
    ggtree::vexpand(.1, -1)+
    scale_fill_gradientn(trans = trans, colors = c("black", "orange", "red"), limits = c(lowlim, NA), name = "TPM") +
    ggtitle(paste(experiment, eps))+
    theme(legend.position = legend)
  return(heattree)
}

stewartan_single <- expression_phylogeny(eps = "Stewartan", offset = 1)

alginate_single <- expression_phylogeny(eps = "Alginate", offset = 0.01)


HA_has_single <- expression_phylogeny(eps = "HA (has)", offset = 1)

ggsave("/mnt/ahd/EPS_PIPELINE/figures/expression_connection/stewartan_single.pdf", stewartan_single, limitsize = FALSE, width = 30, height = 18, dpi = 600)

ggsave("/mnt/ahd/EPS_PIPELINE/figures/expression_connection/HA_has_single.pdf", HA_has_single, limitsize = FALSE, width = 30, height = 18, dpi = 600)

cellulose_NA_single <- expression_phylogeny(eps = "Unclassified cellulose", offset = 1)

ggsave("/mnt/ahd/EPS_PIPELINE/figures/expression_connection/cellulose_NA_single.pdf", cellulose_NA_single, limitsize = FALSE, width = 30, height = 18, dpi = 600)



# expression_phylogeny_plural <- function(eps, experiment = "Anoxic stage", legend = "right", limits = c(1, 25),
#                                  treelayout = "rectangular", trans = "identity", width = 2, offset = 2) {
# 
#   EPS_table_filt <- EPS_table_RSEM %>% filter(operon %in% c(eps), `Processing tank` == experiment) %>% mutate(label = substring(gene_id, 1, nchar(gene_id) - 6)) %>%
#     group_by(label, gene_annotation) %>% summarise(TPM = mean(TPM)) %>% filter(any(TPM >= limits[1])) %>% pivot_wider(names_from = gene_annotation, values_from = TPM) %>% as.data.frame()
# 
#   rownames(EPS_table_filt) <- EPS_table_filt$label
# 
#   EPS_table_filt <- EPS_table_filt %>% select(-label)
# 
# 
#   tree_tib_relevant <- tree %>%
#     as_tibble() %>% filter(label %in% rownames(EPS_table_filt))
# 
#   tree_tib_irrelevant <- tree %>%
#     as_tibble() %>% filter(label %ni% tree_tib_relevant$label)
# 
#   tree_filt <- tree %>% drop.tip(tree_tib_irrelevant$label)
# 
#   rel_abun <- EPS_table_RSEM %>% 
#     mutate(label = substring(gene_id, 1, nchar(gene_id) - 6)) %>% ungroup() %>%
#     select(label, `Relative abundance (%)`) %>% filter(label %in% tree_tib_relevant$label) %>% group_by(label, `Relative abundance (%)`) %>%
#     distinct(label)
#   
#   tree_filt_tib <- tree_filt %>% as_tibble() %>% left_join(rel_abun, by = "label") %>% as.treedata()
#   
#   
#   tree_plot <- ggtree(tree_filt_tib, layout = treelayout) +
#     geom_tiplab( 
#       mapping = aes(node = node, label = paste(label, genus)),
#       geom = "text",
#       size = 5
#     ) + 
#     geom_tippoint(aes(size = `Relative abundance (%)`))
#   
# 
#   EPS_table_plot <- EPS_table_RSEM %>% filter(operon %in% c(eps), `Processing tank` == experiment) %>% mutate(label = substring(gene_id, 1, nchar(gene_id) - 6)) %>%
#     group_by(label, gene_annotation) %>% filter(any(TPM >= limits[1])) %>% summarise(TPM = mean(TPM))
# 
#   heattree <- ggplot(EPS_table_plot, aes(gene_annotation, label)) +
#     geom_tile(aes(fill = TPM)) +
#     ggtitle(experiment)+
#     scale_fill_gradientn(trans = trans, colors = c("black", "orange", "red"), limits = limits, name = "TPM", na.value = "transparent") +
#     theme(axis.text.y = element_blank(),
#           axis.ticks.y = element_blank(),
#           axis.title.y = element_blank(),
#           axis.title.x = element_blank(),
#           plot.title = element_text(size = 10),
#           axis.text.x = element_text(size = 6)
#     )
# 
#   heattree <- heattree %>% insert_left(tree_plot, width = 2.5)
# 
#   for (item in list) {
# 
#     temp <- EPS_table_RSEM %>% filter(operon %in% c(eps), `Processing tank` == item) %>% mutate(label = substring(gene_id, 1, nchar(gene_id) - 6)) %>%
#       group_by(label, gene_annotation) %>% filter(any(TPM >= limits[1])) %>% summarise(TPM = mean(TPM))
# 
#     temp_plot <- ggplot(temp, aes(gene_annotation, label)) +
#       ggtitle(item, )+
#       geom_tile(aes(fill = TPM)) +
#       scale_fill_gradientn(trans = trans, colors = c("black", "orange", "red"), limits = limits, name = "TPM", na.value = "transparent") +
#       theme(axis.text.y = element_blank(),
#             axis.ticks.y = element_blank(),
#             axis.title.y = element_blank(),
#             axis.title.x = element_blank(),
#             plot.title = element_text(size = 10),
#             axis.text.x = element_text(size = 6)
#       )
# 
#     heattree <- heattree %>% insert_right(temp_plot)
# 
#   }
# 
#   return(heattree)
# }
# 
# 
# list <- c("Aerobic stage", "Sidestream tank 1", "Sidestream tank 2",
#           "Return sludge")
# 
# stewartan_plural <- expression_phylogeny_plural(eps = "Burkholderia EPS")
# 
# HA_has_plural <- expression_phylogeny_plural(eps = "HA (has)")
# 
# 
# ggsave("/mnt/ahd/EPS_PIPELINE/figures/expression_connection/stewartan_plural.pdf", stewartan_plural, limitsize = FALSE, width = 30, height = 18, dpi = 600)
# 
# ggsave("/mnt/ahd/EPS_PIPELINE/figures/expression_connection/HA_has_plural.pdf", HA_has_plural, limitsize = FALSE, width = 30, height = 18, dpi = 600)
# 
# expression_phylogeny_plural(eps = "Levan")

# # source("./scripts/R/HQ_MAG_tree.R")
# # saveRDS(tree, file = "/mnt/ahd/EPS_PIPELINE/data/metatranscriptomics/HQtree.rds")
# tree <- readRDS("/mnt/ahd/EPS_PIPELINE/data/raw/HQtree.rds")
# 
# EPS_table_RSEM <- readRDS("/mnt/ahd/EPS_PIPELINE/data/metatranscriptomics/EPS_table_Fullscale_RSEM.rds")
# 
# expression_phylogeny <- function(eps, experiment = "Control", replicate = "A", legend = "right", 
#                                  treelayout = "rectangular", trans = "identity",width = 2, offset = 0.5, time = "T0") {
# 
# EPS_table_filt <- EPS_table_RSEM %>% filter(operon %in% c(eps), Replicate == paste(experiment, replicate, sep = "_"), Time_point == time) %>% mutate(label = substring(gene_id, 1, nchar(gene_id) - 6)) %>% 
#   group_by(label, gene_annotation) %>% summarise(TPM = max(TPM)) %>% pivot_wider(names_from = gene_annotation, values_from = TPM) %>% as.data.frame()
# 
# rownames(EPS_table_filt) <- EPS_table_filt$label 
# 
# EPS_table_filt <- EPS_table_filt %>% select(-label)
# 
# 
# tree_tib_relevant <- tree %>% 
#   as_tibble() %>% filter(label %in% rownames(EPS_table_filt))
# 
# tree_tib_irrelevant <- tree %>% 
#   as_tibble() %>% filter(label %ni% tree_tib_relevant$label)
# 
# tree_filt <- tree %>% drop.tip(tree_tib_irrelevant$label)
# tree_filt_tib <- tree_filt %>% as_tibble()
# 
# 
# tree_plot <- ggtree(tree_filt, layout = treelayout) +
#   geom_tiplab( 
#     mapping = aes(node = node, label = paste(phylum,genus)),
#     geom = "text",
#     size = 3
#   )
# 
# heattree <-gheatmap(tree_plot,EPS_table_filt, offset = offset, colnames_angle = -45, colnames_offset_y = 0.4, width = width, hjust = 0) +
#   ggtree::vexpand(.1, -1)+
#   scale_fill_gradientn(trans = trans, colors = c("black", "orange", "red"), limits = c(10, NA), name = "TPM") +
#   ggtitle(paste(eps,"Substrate:",experiment, sep = " "))+
#   theme(legend.position = legend)
# return(heattree)
# }
# 
# 
# 
# 
# expression_phylogeny("stewartan", time = "T-1", experiment = exp, offset = 3) + 
#   expression_phylogeny("stewartan", time = "T0", experiment = exp) + 
#   expression_phylogeny("stewartan", time = "T1", experiment = exp) + 
#   expression_phylogeny("stewartan", time = "T2", experiment = exp) + 
#   expression_phylogeny("stewartan", legend = "right", time = "T3", experiment = exp)
# 
# expression_phylogeny("HA_streptococcus", time = "T-1", experiment = exp, offset = 3) + 
#   expression_phylogeny("HA_streptococcus", time = "T0", experiment = exp) + 
#   expression_phylogeny("HA_streptococcus", time = "T1", experiment = exp) + 
#   expression_phylogeny("HA_streptococcus", time = "T2", experiment = exp) + 
#   expression_phylogeny("HA_streptococcus", legend = "right", time = "T3", experiment = exp)


