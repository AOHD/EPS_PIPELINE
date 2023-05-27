library("data.table")
library("ggtree")
library("ggtreeExtra")
library("ggnewscale")
library("readxl")
library("treeio")
library("tidyverse")
library("glue")
library("shadowtext")
`%ni%` <- Negate(`%in%`)

##------------------------------------------------------------------------
##  Import percentage of identified genes in each HQ-MAG for all queries  
##------------------------------------------------------------------------
# File with metadata on the query genes, mostly for total genes in query
query_metadata. <- excel_sheets("./data/raw/Query_figur.xlsx") %>%
  sapply(function(X) read_xlsx("./data/raw/Query_figur.xlsx", sheet = X, skip = 1), USE.NAMES = T) %>% 
  lapply(as.data.frame) %>% 
  `[`(!(names(.) %in% c("Abbreveations", "HA_S_pyogenes"))) %>%
  rbindlist(fill = TRUE) %>% 
  mutate(Polysaccheride = ifelse(Polysaccheride == "S88","s88",Polysaccheride),
         Psiblast = ifelse(Psiblast == "S88","s88",Psiblast))

# Percentage of percent identity filtrated genes in each HQ-MAG
psi_perc_filt <- list.files("./output_proximity_filtration/psi_percID_filt/") %>%
  map(function(query){
    name_query <- tools::file_path_sans_ext(query)
    fread(paste0("./output_proximity_filtration/psi_percID_filt/", query)) %>% 
      group_by(ID) %>% 
      summarise(genes_percent = length(unique(Query_label))) %>% 
      mutate(
        genes_percent = 100 * genes_percent/(query_metadata. %>% filter(Psiblast == name_query) %>% nrow)
      ) %>% 
      setNames(c("ID", name_query))
  })



psi_perc_filt <- psi_perc_filt %>%
  reduce(full_join, by = "ID") %>%
  tibble::column_to_rownames(var = "ID")



colnames(psi_perc_filt) <- psi_perc_filt %>% colnames() %>%
  str_replace("alginate", "Alginate") %>% 
  str_replace("HA_Pasteurella", "HA (pmHAS)") %>% 
  str_replace("HA_streptococcus", "HA (has)") %>% 
  str_replace("NulO_merged", "NulO") %>% 
  str_replace("pel_merged", "Pel") %>% 
  str_replace("pnag_pga", "PNAG (pga)") %>% 
  str_replace("B_subtilis_EPS", "B. subtilis EPS") %>% 
  str_replace("pnag_ica", "PNAG (ica)") %>% 
  str_replace("xanthan", "Xanthan") %>% 
  str_replace("psl", "Psl") %>% 
  str_replace("curdlan", "Curdlan") %>% 
  str_replace("diutan", "Diutan") %>% 
  str_replace("gellan2", "Gellan2") %>%
  str_replace("burkholderia_eps", "Burkholderia EPS") %>%
  str_replace("amylovoran", "Amylovoran") %>%
  str_replace("ColA", "Colanic Acid") %>%
  str_replace("salecan", "Salecan") %>%
  str_replace("stewartan", "Stewartan") %>%
  str_replace("vps", "Vibrio polysaccharide") %>%
  str_replace("rhizobium_eps", "Rhizobium EPS") %>%
  str_replace("gellan1", "Gellan1") %>%
  str_replace("acetan", "Acetan") %>%
  str_replace("s88", "s88") %>%
  str_replace("levan", "Levan") %>% 
  str_replace("synechan", "Synechan") %>%
  str_replace("methanolan", "Methanolan") %>%
  str_replace("celluloseI", "Cellulose I") %>% 
  str_replace("celluloseII", "Cellulose II") %>%
  str_replace("celluloseIII", "Cellulose III") %>%
  str_replace("cellulose_Ac", "Acetylated cellulose") %>%
  str_replace("cellulose_NA", "Unclassified cellulose") %>%
  str_replace("succinoglycan", "Succinoglycan") %>%
  str_replace("galactoglucan", "Galactoglucan") %>%
  str_replace("B_fragilis_PS_A", "B. fragilis PS A") %>%
  str_replace("B_fragilis_PS_B", "B. fragilis PS B") %>%
  str_replace("B_pseudomallei_EPS", "B. pseudomallei PS") %>%
  str_replace("cepacian", "Cepacian") %>%
  str_replace("E_faecalis_PS", "E. faecalis PS") %>%
  str_replace("emulsan", "Emulsan") %>%
  str_replace("EPS273", "EPS273") %>%
  str_replace("GG", "GG") %>%
  str_replace("glucorhamnan", "Glucorhamnan") %>%
  str_replace("L_johnsonii_ATCC_33200_EPS_A", "L. johnsonii PS A") %>%
  str_replace("L_johnsonii_ATCC_11506_EPS_B", "L. johnsonii PS B") %>%
  str_replace("L_johnsonii_ATCC_2767_EPS_C", "L. johnsonii PS C") %>%
  str_replace("L_lactis_EPS", "L. lactis PS") %>%
  str_replace("L_plantarum_HePS", "L. plantarum PS") %>%
  str_replace("phosphonoglycan", "Phosphonoglycan")


# Percentage of proximity filtrated genes in each HQ-MAG
psi_proxi_filt <- list.files("./output_proximity_filtration/psi_proxi_filt/") %>%
  map(.f = function(query){
    name_query <- tools::file_path_sans_ext(query)
    fread(paste0("./output_proximity_filtration/psi_proxi_filt/", query)) %>% 
      group_by(ID) %>% 
      summarise(genes_percent = length(unique(Query_label))) %>% 
      mutate(
        genes_percent = 100 * genes_percent/(query_metadata. %>% filter(Psiblast == name_query) %>% nrow)
      ) %>% 
      setNames(c("ID", name_query))
  }
  )

psi_proxi_filt <- psi_proxi_filt %>%
  reduce(full_join, by = "ID") %>% 
  tibble::column_to_rownames(var = "ID")

colnames(psi_proxi_filt) <- psi_proxi_filt %>% colnames() %>% 
  str_replace("alginate", "Alginate") %>% 
  str_replace("HA_Pasteurella", "HA (pmHAS)") %>% 
  str_replace("HA_streptococcus", "HA (has)") %>% 
  str_replace("NulO_merged", "NulO") %>% 
  str_replace("pel_merged", "Pel") %>% 
  str_replace("pnag_pga", "PNAG (pga)") %>% 
  str_replace("B_subtilis_EPS", "B. subtilis EPS") %>% 
  str_replace("pnag_ica", "PNAG (ica)") %>% 
  str_replace("psl", "Psl") %>%
  str_replace("diutan", "Diutan") %>%
  str_replace("xanthan", "Xanthan") %>%
  str_replace("gellan2", "Gellan2") %>%
  str_replace("burkholderia_eps", "Burkholderia EPS") %>%
  str_replace("amylovoran", "Amylovoran") %>%
  str_replace("ColA", "Colanic Acid") %>%
  str_replace("salecan", "Salecan") %>%
  str_replace("stewartan", "Stewartan") %>%
  str_replace("vps", "Vibrio polysaccharide") %>%
  str_replace("acetan", "Acetan") %>%
  str_replace("s88", "s88") %>% 
  str_replace("synechan", "Synechan") %>% 
  str_replace("levan", "Levan") %>%
  str_replace("methanolan", "Methanolan") %>%
  str_replace("celluloseI", "Cellulose I") %>% 
  str_replace("celluloseII", "Cellulose II") %>%
  str_replace("celluloseIII", "Cellulose III") %>%
  str_replace("cellulose_Ac", "Acetylated cellulose") %>%
  str_replace("cellulose_NA", "Unclassified cellulose") %>%
  str_replace("succinoglycan", "Succinoglycan") %>%
  str_replace("galactoglucan", "Galactoglucan") %>%
  str_replace("B_subtilis_EPS", "B. subtilis EPS") %>%
  str_replace("B_fragilis_PS_A", "B. fragilis PS A") %>%
  str_replace("B_fragilis_PS_B", "B. fragilis PS B") %>%
  str_replace("B_pseudomallei_EPS", "B. pseudomallei PS") %>%
  str_replace("cepacian", "Cepacian") %>%
  str_replace("E_faecalis_PS", "E. faecalis PS") %>%
  str_replace("emulsan", "Emulsan") %>%
  str_replace("EPS273", "EPS273") %>%
  str_replace("GG", "GG") %>%
  str_replace("glucorhamnan", "Glucorhamnan") %>%
  str_replace("L_johnsonii_ATCC_33200_EPS_A", "L. johnsonii PS A") %>%
  str_replace("L_johnsonii_ATCC_11506_EPS_B", "L. johnsonii PS B") %>%
  str_replace("L_johnsonii_ATCC_2767_EPS_C", "L. johnsonii PS C") %>%
  str_replace("L_lactis_EPS", "L. lactis PS") %>%
  str_replace("L_plantarum_HePS", "L. plantarum PS") %>%
  str_replace("phosphonoglycan", "Phosphonoglycan")

original_operons <- c("Alginate", "Cellulose I", "Cellulose II", "Curdlan", "Diutan", "HA (has)", "HA (pmHAS)", 
                      "Pel", "PNAG (eps)", "PNAG (ica)", "PNAG (pga)", "Psl", "s88", "Xanthan", "Succinoglycan","Salecan")


psi_proxi_filt_all <- psi_proxi_filt %>%
  mutate(ID = row.names(.)) %>%
  pivot_longer(-ID) %>% filter(name != c("NulO"))

psi_perc_filt_all <- psi_perc_filt %>%
  mutate(ID = row.names(.)) %>%
  pivot_longer(-ID) %>% filter(name %ni% c("gellan", "NulO"))

psi_proxi_filt_all <- psi_proxi_filt_all %>%
  rbind(
    psi_perc_filt_all %>%
      group_by(name) %>%
      # head(n=16) %>%
      mutate(value = NA) %>%
      filter(!(name %in% unique(psi_proxi_filt_all$name))) %>%
      ungroup()
  ) %>%
  mutate(
    Pathway = case_when(
      name == "Levan" ~ "Sucrase-dependent",
      name %in% c("Cellulose I", "Cellulose II", "Cellulose III", "Acetylated cellulose", "Unclassified cellulose", "Alginate",
                  "Curdlan", "HA (pmHAS)", "HA (has)", "Pel", "PNAG (ica)", "PNAG (pga)") ~ "Synthase-dependent",
      name == "Phosphonoglycan" ~ "PEP-mutase dependent",
      TRUE ~ "Wzx/Wzy-dependent"
    ),
    color = case_when(
      Pathway == "Sucrase-dependent" ~ "Blue",
      Pathway == "Synthase-dependent" ~ "Green",
      Pathway == "Wzx/Wzy-dependent" ~ "Red",
      Pathway == "PEP-mutase dependent" ~ "Purple",
      TRUE ~ "Black"
    )
  )



psi_proxi_filt_all <- psi_proxi_filt_all %>% group_by(Pathway) %>% arrange(name, .by_group = TRUE) 


mycols <- as.vector(levels(psi_proxi_filt_all$color))

# psi_proxi_filt_all$name <- factor(psi_proxi_filt_all$name,
#                                        levels=c(unique(psi_proxi_filt_all$name)))
# 
# psi_perc_filt_all$name <- factor(psi_perc_filt_all$name,
#                                       levels=c(unique(psi_perc_filt_all$name)))
##


# psi_proxi_filt_original <- psi_proxi_filt %>% 
#   mutate(ID = row.names(.)) %>% 
#   pivot_longer(-ID) %>% filter(name != c("NulO")) %>% filter(name %in% original_operons)
# 
# psi_perc_filt_original <- psi_perc_filt %>% 
#   mutate(ID = row.names(.)) %>% 
#   pivot_longer(-ID) %>% filter(name %ni% c("gellan", "NulO")) %>% filter(name %in% original_operons)
# 
# psi_proxi_filt_original <- psi_proxi_filt_original %>%
#   rbind(
#     psi_perc_filt_original %>%
#       group_by(name) %>%
#       # head(n=16) %>%
#       mutate(value = NA) %>%
#       filter(!(name %in% unique(psi_proxi_filt_original$name))) %>% 
#       ungroup()
#   )
# 
# psi_proxi_filt_original$name <- factor(psi_proxi_filt_original$name, 
#                               levels=c(unique(psi_proxi_filt_original$name)))
# 
# psi_perc_filt_original$name <- factor(psi_perc_filt_original$name, 
#                              levels=c(unique(psi_perc_filt_original$name)))
# 
# 
# 
# 
# psi_proxi_filt_new <- psi_proxi_filt %>% 
#   mutate(ID = row.names(.)) %>% 
#   pivot_longer(-ID) %>% filter(name != c("NulO")) %>% filter(name %ni% original_operons)
# 
# psi_perc_filt_new <- psi_perc_filt %>% 
#   mutate(ID = row.names(.)) %>% 
#   pivot_longer(-ID) %>% filter(name %ni% c("gellan", "NulO")) %>% filter(name %ni% original_operons)
# 
# 
# psi_proxi_filt_new <- psi_proxi_filt_new %>%
#   rbind(
#     psi_perc_filt_new %>%
#       group_by(name) %>%
#       # head(n=16) %>%
#       mutate(value = NA) %>%
#       filter(!(name %in% unique(psi_proxi_filt_new$name))) %>% 
#       ungroup()
# )
# 
# psi_proxi_filt_new$name <- factor(psi_proxi_filt_new$name, 
#                               levels=c(unique(psi_proxi_filt_new$name)))
# 
# psi_perc_filt_new$name <- factor(psi_perc_filt_new$name, 
#                               levels=c(unique(psi_perc_filt_new$name)))

##---------------------------------------------------------------
##            Import meta information for each HQ-MAG            
##---------------------------------------------------------------
# Import phylum of HQ-MAGs                    
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
    
    species = str_extract(midas4_tax, "s:[^,]*"),
    species = str_remove(species, "s:"),
    species = str_remove(species, ";")
    
  ) %>% 
  rename(ID = bin) %>%
  select(ID, phylum, class, order, family, genus, species) %>% 
  setNames(c("label", "phylum", "class", "order", "family", "genus", "species")) %>%
  filter(!is.na(label), phylum != "Halobacterota")


##---------------------------------------------------------------
##                         Creating tree                         
##---------------------------------------------------------------
# Creating tibble with tree information
# tree_tib_old <- read.tree("./tree_results/gtdb_31102022/HQMAG_31102022.tree") %>% 
#   as_tibble() 
# 
# tree_tib_relevant <- read.tree("./tree_results/gtdb_16112022/infer/gtdbtk.bac120.decorated.tree") %>% 
#   as_tibble() %>% filter(grepl("_BAT",label)|grepl("_MAX",label)|grepl("glom_",label))
# 
# tree_tib_irrelevant <- read.tree("./tree_results/gtdb_16112022/infer/gtdbtk.bac120.decorated.tree") %>% 
#   as_tibble() %>% filter(label %ni% tree_tib_relevant$label)


tree_tib <- read.tree("./tree_results/gtdb_16112022/infer/gtdbtk.bac120.decorated_chloroflex.tree") %>% #drop.tip(tree_tib_irrelevant$label) %>%
  as_tibble() %>% #mutate(label = gsub('\'', '', label)) %>%
  full_join(phylum, by = "label") %>%  
  filter(!is.na(node)) %>% 
  as.treedata() %>% as_tibble() %>% 
  # Add phylum to intermediate nodes (only those with unique offspring phylum)
  mutate(
    offspring_phylum = map(node, function(x) offspring(., x) %>% pull(phylum) %>% unique %>% `[`(!is.na(.))),
    parent_offspring_phylum = map(parent, function(x) {
      offspring(., x) %>% 
        pull(phylum) %>% 
        unique() %>% 
        `[`(!is.na(.))
    }),
    phylum = case_when(
      offspring_phylum %>% map(length) == 0 ~ phylum,
      offspring_phylum %>% map(length) == 1 ~ as.character(offspring_phylum),
      TRUE ~ "NA"),
    phylum = ifelse(phylum == "NA", NA, phylum),
    phylum_ancestor = case_when(
      parent_offspring_phylum %>% map(length) > 1 & offspring_phylum %>% map(length) == 1 ~ as.character(offspring_phylum),
      TRUE ~ "NA",
    )
  ) %>% 
  select(-offspring_phylum)


# Selecting ancestors
tree_tib_ancestors <- tree_tib %>% filter(phylum_ancestor != "NA") %>% group_by(phylum_ancestor) %>% filter(parent == min(parent))

ancestor_list <- list()

ancestor_list[tree_tib_ancestors$phylum_ancestor] <- tree_tib_ancestors$parent

tree_tib1 <- tree_tib

for (Ancestor in names(ancestor_list)) {
  tree_tib1 <- tree_tib1 %>% mutate(
    phylum = ifelse(node == ancestor_list[[Ancestor]] & phylum %ni% c("Chloroflexi", "Patescibacteria", "Acidobacteriota", "Planctomycetota"),
                    Ancestor, phylum),
    phylum_ancestor = ifelse(node == ancestor_list[[Ancestor]] & phylum %ni% c("Chloroflexi", "Patescibacteria","Acidobacteriota", "Planctomycetota"),
                             Ancestor, phylum_ancestor)
  )
}


# Creating treedata object
tree <- tree_tib1 %>%
  as.treedata() 

##---------------------------------------------------------------
##                         Plotting tree                         
##---------------------------------------------------------------

# Tree with inspiration from https://yulab-smu.top/treedata-book/chapter10.html
perc_and_proxi_fruit_layers <- function(proxi, title){
  
  colors <- c(
    "#ffc885",
    "#f29d57",
    "#e56e33",
    "#d7301d"
  )
    fruit_list <- geom_fruit_list(
      geom_fruit(
        data = proxi, 
        geom = geom_tile,
        pwidth = 0.8*length(unique(proxi$name)),
        mapping = aes(y = ID, x = name, fill = value, width = 20, height = 2),
        axis.params = list(axis = "x", text.size = 20, text.angle = 45, hjust = 0, vjust = 1, 
                           title = title, title.size = 20),               # axis.params = list(axis = "x", text.size = 2.5, text.angle = -90, hjust = 0, vjust = 0.5),
        grid.params = list(vline = FALSE, color = "black", alpha = 1, size = 0.6),
        offset = 2
      ),
      scale_fill_stepsn(
        colors = "orange",
        na.value = "transparent",  
        n.breaks = 6,
        breaks = waiver(),
        limits = c(20, 100)
      )
    )
    return(fruit_list)
}


# rectangular

plot_operon_HQ_mag_original <- function(proxi = psi_proxi_filt_all, savename = "original", phylumNO = 10, layout = "rectangular"){
  phylum_displayed <- phylum$phylum %>% table %>% `[`(order(.)) %>% tail(phylumNO) %>% names()
  
  tree_plot <- ggtree(
    tree, 
    layout = layout, 
    branch.length = "none",
    lwd = 0.5
  ) +
    # geom_hilight(
    #   data = filter(tree_tib1, phylum_ancestor != "NA") %>% group_by(phylum_ancestor) %>% filter(node == min(node)) %>%
    #     rename(Phylum_Ancestor = phylum_ancestor) %>%
    #     filter(phylum %in% phylum_displayed), 
    #   aes(node = node), fill = "gray",
    #   extendto = 4.5, alpha = 0.4) +
    geom_cladelab(
      data = filter(tree_tib1, phylum_ancestor != "NA") %>% group_by(phylum_ancestor) %>% filter(node == min(node)) %>%
        filter(phylum %in% phylum_displayed) %>% rename(Phylum = phylum), 
      mapping = aes(node = node, label = phylum),
      angle = 0, barsize = 2,
      geom = "text", fontsize = 18, horizontal = TRUE,
      align = TRUE, offset.text = 2
    ) +
    scale_color_manual(guide = "none") +
    # perc_and_proxi_fruit_layers(proxi %>% filter(Pathway == "Synthase-dependent"), title = "Synthase-dependent") +
    perc_and_proxi_fruit_layers(proxi %>% filter(Pathway == "Wzx/Wzy-dependent"), title = "Wzx/Wzy-dependent") +
    # perc_and_proxi_fruit_layers(proxi %>% filter(Pathway == "Sucrase-dependent"), title = "Sucrase-dependent") +
    ggtree::hexpand(.3, -1) +
    ggtree::hexpand(.3, 1) +
    ggtree::vexpand(.1, 1) +
     scale_y_reverse() + 
     coord_flip(clip = "off")  +
    theme(legend.position="none")
  
  ggsave(paste0("./figures/trees/HQ_MAG_tree_fan_", savename,".pdf"), width =40, height = 60, limitsize = FALSE, dpi = 150,
         plot = tree_plot
  )
}

##---------------------------------------------------------------
##                         Trees of genes in each system                         
##---------------------------------------------------------------
# Tree with inspiration from https://yulab-smu.top/treedata-book/chapter10.html
perc_and_proxi_fruit_layers_genes <- function(perc_data, proxi_data){
  # Colors of perc filt
  perc_colors <- c(
    "#fefeff",
    "#d9dbff",
    "#aebaff",
    "#7a9bff",
    "#007eff")
  # Colors proxi filt
  proxi_colors <- c(
    "#fffefe",
    "#ffbfbf",
    "#ff7f7f",
    "#ff4040",
    "#ff0000"
  )
  
  fruit_list <- geom_fruit_list(
    geom_fruit(
      data = perc_data, 
      geom = geom_tile,
      pwidth = 0.8,
      mapping = aes(y = ID, x = name, fill = value),
      axis.params = list(axis = "x", text.size = 2, text.angle = -90, hjust = 0, vjust = 0),               # axis.params = list(axis = "x", text.size = 2.5, text.angle = -90, hjust = 0, vjust = 0.5),
      grid.params = list(vline = FALSE, color = "gray60", alpha = 1),
      offset = 0.1
    ),
    scale_fill_stepsn(
      colors = perc_colors, 
      na.value = "transparent",  
      n.breaks = 6,
      breaks = waiver(),
      limits = c(20, 50),
      guide = guide_bins(
        title = "Highest Percent Identity of Gene \n Before Proximity Filtration",
        title.position = "top",
        show.limits = TRUE,
        title.hjust = 0.5,
        title.vjust = 0.5,
        keywidth = unit(13, "mm"),
        ticks.colour = "black", 
        frame.colour = "black")
    ),
    new_scale_fill(),
    geom_fruit(
      data = proxi_data,
      geom = geom_tile,
      pwidth = 0.8,
      mapping = aes(y = ID, x = name, fill = value),
      axis.params = list(axis = "x", text.size = 2, text.angle = -90, hjust = 0, vjust = 0),               # axis.params = list(axis = "x", text.size = 2.5, text.angle = -90, hjust = 0, vjust = 0.5),
      grid.params = list(vline = FALSE, color = "gray60", alpha = 1),
      offset = 0.1
    ),
    scale_fill_stepsn(
      colors = proxi_colors,
      na.value = "transparent",
      n.breaks = 6,
      breaks = waiver(),
      limits = c(20, 50),
      guide = guide_bins(
        title = "Highest Percent Identity of Gene \n After Proximity Filtration",
        title.position = "top",
        show.limits = TRUE,
        title.hjust = 0.5,
        title.vjust = 0.5,
        keywidth = unit(13, "mm"),
        ticks.colour = "black",
        frame.colour = "black")
    )
  )
  return(fruit_list)
}

dir.create("./figures/trees/HQ_MAG_genes", showWarnings = FALSE)
plot_genes_HQ_mag <- function(eps, phylumNO = 20){
  phylum_displayed <- phylum$phylum %>% table %>% `[`(order(.)) %>% tail(phylumNO) %>% names()
  message(glue("Processing {eps}"))
  data_perc <- fread(glue("./output_proximity_filtration/psi_percID_filt/{eps}.tsv")) %>% 
    select(Query_label, ID, Percent_identity) %>%
    dplyr::rename(name = Query_label, value = Percent_identity) %>% 
    group_by(name, ID) %>% 
    filter(value == max(value)) 
  missing_genes <- query_metadata. %>% 
    filter(Psiblast %in% eps & !(Genename %in% data_perc$name)) %>% 
    pull(Genename) %>% 
    unique()
  if (length(missing_genes) != 0){
    data_perc <- data_perc %>% 
      dplyr::bind_rows(
        data.frame(
          name = missing_genes,
          ID = data_perc$ID[1],
          value = NA)
      )
  }
  
  if(file.exists(glue("./output_proximity_filtration/psi_proxi_filt/{eps}.tsv"))) {
    data_proxi <- fread(glue("./output_proximity_filtration/psi_proxi_filt/{eps}.tsv")) %>% 
      select(Query_label, ID, Percent_identity) %>% 
      dplyr::rename(name = Query_label, value = Percent_identity) %>% 
      group_by(name, ID) %>% 
      filter(value == max(value))
    data_proxi <- data_proxi %>%
      dplyr::bind_rows(
        data_perc %>%
          group_by(name) %>%
          filter( row_number() == 1 & !(name %in% unique(data_proxi$name)) ) %>%
          mutate(value = NA)
      )
    
    NA_bins <- tree_tib %>% filter(label %ni% data_proxi$ID) %>% select(ID = label)
    
    
    for (gene in unique(data_proxi$name)) {
      data_proxi <- data_proxi %>%
        dplyr::bind_rows(
          data.frame(
            name = gene,
            ID = NA_bins,
            value = NA
          )
        )
    }
    
  } else {
    data_proxi <- data.frame(matrix(ncol = 3, nrow = 0)) %>% 
      setNames(c("ID", "name", "value"))
  }
  tree_plot <- 
    ggtree(
      tree, 
      layout = "fan", 
      lwd = 0.1, 
      open.angle = 20
    ) +
    geom_hilight(
      data = filter(tree_tib1, phylum_ancestor != "NA") %>% group_by(phylum_ancestor) %>% filter(node == min(node)) %>%
        rename(Phylum_Ancestor = phylum_ancestor) %>%
        filter(phylum %in% phylum_displayed), 
      aes(node = node), 
      fill = "gray", extendto = 4.5, alpha = 0.4
    )+
    geom_cladelab(
      data = filter(tree_tib1, phylum_ancestor != "NA") %>% group_by(phylum_ancestor) %>% filter(node == min(node)) %>%
        filter(phylum %in% phylum_displayed) %>% rename(Phylum = phylum), 
      mapping = aes(node = node, label = phylum),
      angle = "auto", barsize = NA,
      geom = "text", fontsize = 2, horizontal = TRUE, 
      align = TRUE, fontface  = 0.8, offset = -1
    ) +
    scale_color_manual(guide = "none") +
    perc_and_proxi_fruit_layers_genes(perc_data = data_perc,
                                      proxi_data = data_proxi) +
    theme(legend.position = "bottom",
          plot.margin = unit(c(-20,-30,0,-30), "mm"),
          legend.margin = margin(t = -60, r = 10, l = 10))
  ggsave(glue("./figures/trees/HQ_MAG_genes/HQ_MAG_tree_fan_percid_{eps}.pdf"), width = 9, height = 10, limitsize = FALSE,
         plot = tree_plot
  )
}

