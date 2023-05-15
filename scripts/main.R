library(here)

here("/user_data/ahd/EPS_PIPELINE")
##############################################
# The main proximity filtration
##############################################
source("/user_data/ahd/EPS_PIPELINE/scripts/R/proximity_filtration.R")
# Default: min_genes=2 and perc_id=20
proximity_filtration("alginate", 
                     min_genes = 6)
proximity_filtration("psl", 
                     min_genes = 7)
proximity_filtration("pel_merged",
                     min_genes = 6)
proximity_filtration("cellulose",
                     same_database = "celluloseI",
                     min_genes = 2,
                     essential_genes = "bcsD")
proximity_filtration("cellulose",
                     same_database = "celluloseII", 
                     min_genes = 2,
                     essential_genes = c("bcsE", "bcsG"))
proximity_filtration("cellulose", 
                     same_database = "celluloseIII",
                     min_genes = 2,
                     essential_genes = c("bcsK"))
proximity_filtration("cellulose", 
                     same_database = "cellulose_Ac",
                     min_genes = 2,
                     essential_genes = c("algF", "algI", "algJ", "algX"))
proximity_filtration("cellulose", 
                     same_database = "cellulose_NA",
                     min_genes = 2,
                     exclude_gene = c("algF", "algI", "algJ", "algX", "bcsK", "bcsE", "bcsG", "bcsD"))
proximity_filtration("succinoglycan",
                     min_genes = 10)
proximity_filtration("xanthan",
                     min_genes = 6)
proximity_filtration("curdlan",
                     min_genes = 2)
proximity_filtration("pnag_pga",
                     min_genes = 3)
proximity_filtration("pnag_ica",
                     min_genes = 3)
proximity_filtration("B_subtilis_EPS", 
                     min_genes = 3)
proximity_filtration("diutan",
                     min_genes = 10)
proximity_filtration("gellan",
                     min_genes = 10)
proximity_filtration("s88", 
                     min_genes = 9)
proximity_filtration("NulO_merged", 
                     min_genes = 3)
proximity_filtration("HA_Pasteurella",
                     min_genes = 5,
                     perc_id = 20)
proximity_filtration("HA_streptococcus", 
                     min_genes = 3, 
                     perc_id = 20)
proximity_filtration("acetan", 
                     min_genes = 8)
proximity_filtration("amylovoran", 
                     min_genes = 6)
proximity_filtration("ColA", 
                     min_genes = 9)
proximity_filtration("gellan1", 
                     min_genes = 2)
proximity_filtration("gellan2", 
                     min_genes = 9)
proximity_filtration("salecan", 
                     min_genes = 6)
proximity_filtration("rhizobium_eps", 
                     min_genes = 9)
proximity_filtration("stewartan", 
                     min_genes = 6)
proximity_filtration("vps", 
                     min_genes = 9)
proximity_filtration("burkholderia_eps", 
                     min_genes = 6)
proximity_filtration("levan", 
                     min_genes = 2)
proximity_filtration("synechan", 
                     min_genes = 10)
proximity_filtration("methanolan", 
                     min_genes = 11)
proximity_filtration("galactoglucan", 
                     min_genes = 9)

proximity_filtration("B_fragilis_PS_A", 
                     min_genes = 7)
proximity_filtration("B_fragilis_PS_B", 
                     min_genes = 10)
proximity_filtration("GG", 
                     min_genes = 8)
proximity_filtration("B_pseudomallei_EPS", 
                     min_genes = 10)
proximity_filtration("phosphonoglycan", 
                     min_genes = 10)
proximity_filtration("E_faecalis_PS", 
                     min_genes = 8)
proximity_filtration("glucorhamnan", 
                     min_genes = 10)
proximity_filtration("L_plantarum_HePS", 
                     min_genes = 10)
proximity_filtration("L_johnsonii_ATCC_33200_EPS_A", 
                     min_genes = 7)
proximity_filtration("EPS273", 
                     min_genes = 8)
proximity_filtration("L_johnsonii_ATCC_11506_EPS_B", 
                     min_genes = 10)
proximity_filtration("emulsan", 
                     min_genes = 10)
proximity_filtration("L_johnsonii_ATCC_2767_EPS_C", 
                     min_genes = 10)
proximity_filtration("L_lactis_EPS", 
                     min_genes = 8)
proximity_filtration("cepacian", 
                     min_genes = 10)


##############################################
# Plotting of operons from results
# Remember to run ips before this
##############################################
source("/mnt/ahd/EPS_PIPELINE/scripts/R/plot_operon.R")

plot_operon("acetan", width = 1, plot_title = "Acetan")
plot_operon("alginate", width = 1, plot_title = "Alginate")
plot_operon("amylovoran", width = 1, plot_title = "Amylovoran")
plot_operon("cellulose", width = 1,
            same_database = "cellulose_Ac",
            query_title = "Acetylated cellulose",
            plot_title = "Acetylated cellulose")
plot_operon("cellulose", width = 1,
            same_database = "celluloseII",
            query_title = "Cellulose II",
            plot_title = "Cellulose II")
plot_operon("cellulose", width = 1,
            same_database = "celluloseIII",
            query_title = "Cellulose III",
            plot_title = "Cellulose III")
plot_operon("cellulose", width = 1,
            same_database = "cellulose_NA",
            query_title = "Unclassified cellulose",
            plot_title = "Unclassified cellulose")
plot_operon("ColA", width = 1, plot_title = "Colanic acid")
plot_operon("gellan2", width = 1, plot_title = "Gellan")
plot_operon("HA_Pasteurella", width = 1, plot_title = "HA (pmHAS)")
plot_operon("HA_streptococcus", width = 1, plot_title = "HA (has)")
plot_operon("pel_merged", width = 1, plot_title = "Pel")
plot_operon("B_subtilis_EPS", width = 1, plot_title = "B. subtilis EPS")
plot_operon("pnag_pga", width = 1, plot_title = "PNAG (pga)")
plot_operon("psl", width = 1, plot_title = "Psl")
plot_operon("s88", width = 1, plot_title = "S88")
plot_operon("salecan", width = 1, plot_title = "Salecan")
plot_operon("stewartan", width = 1, plot_title = "Stewartan")
plot_operon("vps", width = 1, plot_title = "Vibrio polysaccharide")
plot_operon("xanthan", width = 1, plot_title = "Xanthan")
plot_operon("burkholderia_eps", width = 1, plot_title = "Burkholderia EPS")
plot_operon("diutan", width  = 2, plot_title = "Diutan")
plot_operon("levan", width = 1, textsize = 1.5)
plot_operon("synechan", width = 1, plot_title = "Synechan")
plot_operon("methanolan", width = 1, plot_title = "Methanolan")
plot_operon("succinoglycan", width = 1, plot_title = "Succinoglycan")
plot_operon("E_faecalis_PS", width = 1, plot_title = "E. faecalis PS")
plot_operon("galactoglucan", width = 1, plot_title = "Galactoglucan")
plot_operon("B_fragilis_PS_A", width = 1, plot_title = "B. fragilis PS A")
plot_operon("B_fragilis_PS_B", width = 1, plot_title = "B. fragilis PS B")
plot_operon("cepacian", width = 1, plot_title = "Cepacian")
plot_operon("emulsan", width = 1, plot_title = "Emulsan")
plot_operon("GG", width = 1, plot_title = "GG")
plot_operon("L_plantarum_HePS", width = 1, plot_title = "L. plantarum PS")
plot_operon("glucorhamnan", width = 1, plot_title = "Glucorhamnan")

#Expression operon plots

# plot_operon("stewartan", width = 4, expression = TRUE, name_addon = "_expressers", query_title ="Stewartan operons with significant expression",
# mags = c("Ega_18-Q3-R5-49_MAXAC.001", "AalW_18-Q3-R10-53_BAT3C.524", "AalE_18-Q3-R2-46_BATAC.251"), limits = c(0.1, NA))

plot_operon("stewartan", width = 1.75, expression = TRUE, article_plot_domain = TRUE,
            name_addon = "_expressers_domains", query_title ="Stewartan operons with significant expression",
            mags = c("Ega_18-Q3-R5-49_MAXAC.001", "AalW_18-Q3-R10-53_BAT3C.524", "AalE_18-Q3-R2-46_BATAC.251"), limits = c(1, NA))

# plot_operon("HA_streptococcus", width = 1, expression = TRUE, name_addon = "_expressers", query_title = "HA (has) operons with significant expression",
# mags = c("Viby_18-Q3-R106-67_BAT3C.52", "Lyne_18-Q3-R50-59_MAXAC.006", "glom_0465", "Ega_18-Q3-R5-49_MAXAC.001", "glom_0297",
#          "Fred_18-Q3-R57-64_BATAC.284", "glom_0076"), limits = c(0.1, NA))

plot_operon("HA_streptococcus", width = 1.75, expression = TRUE, article_plot_domain = TRUE,
            name_addon = "_expressers_domains", query_title = "HA (has) operons with significant expression",
            mags = c("Viby_18-Q3-R106-67_BAT3C.52", "Lyne_18-Q3-R50-59_MAXAC.006", "glom_0465", "Ega_18-Q3-R5-49_MAXAC.001", "glom_0297",
                     "Fred_18-Q3-R57-64_BATAC.284", "glom_0076"), limits = c(1, NA))


# plot_operon("cellulose", same_database = "cellulose_NA", width = 1, expression = TRUE, name_addon = "_expressers_anox", query_title = "Unclassified cellulose operons with significant expression",
#             mags = c("AalW_18-Q3-R10-53_BAT3C.524", "glom_0454", "glom_0299", "Ega_18-Q3-R5-49_BAT3C.159", "glom_0288", "glom_0360"), 
#             limits = c(0.1, NA))

plot_operon("cellulose", same_database = "cellulose_NA", width = 2, expression = TRUE, article_plot_domain = TRUE,
            name_addon = "_expressers_anox_domains", query_title = "Unclassified cellulose operons with significant expression",
            mags = c("AalW_18-Q3-R10-53_BAT3C.524", "glom_0454", "glom_0299", "Ega_18-Q3-R5-49_BAT3C.159"), 
            limits = c(1, NA))


# plot_operon("cellulose", same_database = "cellulose_NA", width = 1, expression = TRUE, name_addon = "_expressers_ST1", query_title = "Unclassified cellulose operons with significant expression",
#             mags = c("AalW_18-Q3-R10-53_BAT3C.524", "glom_0454", "glom_0299", "Ega_18-Q3-R5-49_BAT3C.159", "glom_0288", "glom_0360"), 
#             limits = c(0.1, 10), Processing_tank = "Sidestream tank 1")

plot_operon("cellulose", same_database = "cellulose_NA", width = 2, expression = TRUE, article_plot_domain = TRUE,
            name_addon = "_expressers_ST1_domains", query_title = "Unclassified cellulose operons with significant expression",
            mags = c("AalW_18-Q3-R10-53_BAT3C.524", "glom_0454", "glom_0299", "Ega_18-Q3-R5-49_BAT3C.159"), 
            limits = c(1, NA), Processing_tank = "Sidestream tank 1")


# plot_operon("levan", width = 1, expression = TRUE, name_addon = "_expressers", query_title = "Levan operons with significant expression",
#             mags = c("glom_0132", "Lyne_18-Q3-R50-59_BATAC.272"), Processing_tank = "Sidestream tank 2", limits = c(0.1, 10))



# plot_operon("alginate", width = 1, expression = TRUE, name_addon = "_expressers", query_title = "Alginate operons with significant expression",
#             mags = c("glom_0322"), limits = c(0.1, 10), Processing_tank = "Aerobic stage")

plot_operon("alginate", width = 1.5, expression = TRUE, article_plot_domain = TRUE,
            name_addon = "_expressers_domains", query_title = "Alginate operons with significant expression",
            mags = c("glom_0322"), limits = c(1, NA), Processing_tank = "Aerobic stage")


# plot_operon("E_faecalis_PS", width = 1, expression = TRUE, name_addon = "_expressers", query_title = "E. faecalis PS operons with significant expression",
#             mags = c("glom_0359", "glom_0465", "Lyne_18-Q3-R50-59_MAXAC.006", "glom_0073"), limits = c(0.1, NA), Processing_tank = "Aerobic stage")

plot_operon("E_faecalis_PS", width = 2.60, expression = TRUE, article_plot_domain = TRUE,
            name_addon = "_expressers_domains", query_title = "E. faecalis PS operons with significant expression",
            mags = c("glom_0359", "glom_0465", "Lyne_18-Q3-R50-59_MAXAC.006", "glom_0073"), limits = c(1, NA), Processing_tank = "Aerobic stage")

# 
# plot_operon("cellulose", same_database = "cellulose_Ac", width = 1, expression = TRUE, name_addon = "_expressers", query_title = "Acetylated cellulose operons with significant expression",
#             mags = c("glom_0479"), 
#             limits = c(0.1, 10))


# plot_operon(c("s88"), width = 1, expression = TRUE, name_addon = "_expressers", query_title = "S88 operons with significant expression",
#             mags = c("glom_0021", "glom_0527", "Ega_18-Q3-R5-49_BATAC.204", "Aalw_18-Q3-R10-53_BATAC.245"), 
#             limits = c(0.1, 10))
# plot_operon(c("s88"), width = 2.60, expression = TRUE, article_plot_domain = TRUE,
#             name_addon = "_expressers_domains", query_title = "S88 operons with significant expression",
#             mags = c("glom_0021", "glom_0527", "Ega_18-Q3-R5-49_BATAC.204", "Aalw_18-Q3-R10-53_BATAC.245"), 
#             limits = c(1, NA))


plot_operon(c("galactoglucan"), width = 2.60, expression = TRUE, article_plot_domain = TRUE,
            name_addon = "_expressers_domains", query_title = "Galactoglucan operons with significant expression",
            mags = c("glom_0359"), 
            limits = c(1, NA))

plot_operon(c("glucorhamnan"), width = 2.60, expression = TRUE, article_plot_domain = TRUE,
            name_addon = "_expressers_domains", query_title = "Glucorhamnan gene cluster with significant expression",
            mags = c("glom_0073"), 
            limits = c(1, NA))



##############################################
# Plotting trees and metagenomic results
# Remember to create .tree file before this
##############################################
source("./scripts/R/HQ_MAG_tree.R")

plot_operon_HQ_mag_original(layout = "rectangular", psi_proxi_filt_all, savename = "rectangular_all_wzxwzynoflip")


c(
  "acetan",
  "alginate",
  "amylovoran",
  "burkholderia_eps",
  "cellulose1",
  "cellulose2",
  "ColA",
  "diutan",
  "gellan2",
  "HA_Pasteurella",
  "HA_streptococcus",
  "NulO_merged",
  "pel_merged",
  "pnag_eps",
  "pnag_pga",
  "psl",
  "s88",
  "salecan",
  "stewartan",
  "vps",
  "xanthan",
  "levan"
) %>% 
  map(plot_genes_HQ_mag)

plot_genes_HQ_mag("celluloseIII", phylumNO = 10)

# Plotting RSEM and DESeq2 expression overview and getting EPS_table



##############################################
# RSEM expression plots
##############################################

###############
# RSEM full-scale expression
###############

summarized_RSEM_Fullscale <- ggplot(EPS_table_RSEM_summarized %>% filter(operon %ni% c("NulO", "Xanthan", "Cepacian", "*B. fragilis* PS B")), aes(`Processing tank`, TPM, color = `Processing tank`)) +
  geom_beeswarm(cex = 3, size = 3) +
  facet_wrap(~operon, scales = "free", ncol = 4) +
  xlab("") + ylab("TPM")+
  theme_bw() +
  ylim(c(0,NA)) +
  # scale_color_brewer(palette="Set2") +
  scale_y_continuous(breaks = pretty_breaks(n=4)) +
  theme(strip.text=element_markdown(size=24, colour = "black"),
        axis.text.x = element_blank(),
        axis.text = element_text(size = 22, colour = "black"),
        axis.title = element_text(size = 26, face = "bold", colour = "black"),
        legend.text = element_text(size = 26, colour = "black"),
        legend.title = element_text(size = 26, face ="bold", colour = "black"),
        legend.position = "right") +
  guides(color = guide_legend(override.aes = list(size = 5),
                              ))


ggsave("/mnt/ahd/EPS_PIPELINE/figures/expression/expression_RSEM_Fullscale_summarized_sum.pdf", summarized_RSEM_Fullscale, limitsize = FALSE, width = 20, height = 20, dpi = 300)

