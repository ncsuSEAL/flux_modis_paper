#************************************************************************************
# Description: Growing season length among biome types
# Author: Xiaojie(J) Gao
# Date: 2022-03-04
#************************************************************************************
rm(list=ls())

source("Code/base.R")
source("Code/mod_base.R")

color_pal <- brewer.pal(8, "Dark2")

modis_color <- color_pal[4]
flux_color <- color_pal[3]



gsl_plot <- ggplot(north_sites_dt) +
    geom_boxplot(aes(x = IGBP, y = m_gsl, fill = modis_color), 
        position = position_nudge(x = -0.2), width = 0.3) +
    geom_boxplot(aes(x = IGBP, y = f_gsl, fill = flux_color), 
        position = position_nudge(x = 0.2), width = 0.3) +
    scale_fill_manual(name = "", values = c("#E7298A" = "#E7298A", 
        "#7570B3" = "#7570B3"), labels = c("MODIS EVI2", "FLUX GPP"), 
        guide = "legend") +
    xlab("Biome type") +
    ylab("Growing season length (day)") +
    theme_article() +
    theme(legend.position = c(0.9, 0.93)) +
    guides(fill = guide_legend(
        keyheight = 0.2,
        default.unit = "inch"
    ))

gsl_diff_plot <- ggplot(north_sites_dt) +
    geom_boxplot(aes(x = IGBP, y = f_gsl - m_gsl), width = 0.3) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    ylab("Growing season length difference (day)") +
    xlab("Biome type") +
    annotate("text", x = 0.5, y = 100, label = "Difference = FLUX - MODIS", 
        hjust = 0, vjust = 0) +
    theme_article() +
    scale_fill_identity()

compose <- gsl_plot / gsl_diff_plot + plot_annotation(tag_levels = "a")

# ggsave(file = "0_Manuscript/SI/SI_Fig_9.png", plot = compose, width = 8, 
    # height = 6, bg = "white")
ggsave(file = "Output/gsl_biomes.png", plot = compose, width = 8, 
    height = 6, bg = "white")
