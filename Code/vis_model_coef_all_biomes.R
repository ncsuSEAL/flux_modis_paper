#************************************************************************************
# Description: Model coefficients for all biome types
# Author: Xiaojie(J) Gao
# Date: 2022-03-04
#************************************************************************************
rm(list=ls())

source("Code/base.R")
source("Code/mod_base.R")
library(ggplot2)
library(patchwork)
library(egg)


color_pal <- brewer.pal(8, "Dark2")

modis_color <- color_pal[4]
flux_color <- color_pal[3]


load(file.path("Pipeline/best_model_fit_VUT.RData"))

gsl_coef_fig <- ggplot() +
    # add modis
    geom_pointrange(aes(x = IGBP, ymin = lwr, y = med, ymax = upr, 
        color = modis_color), shape = 16, alpha = 0.3, size = 0.3, 
        data = m_model$beta2, position = position_dodge2(width = 0.5)) +
    geom_pointrange(aes(x = IGBP, ymin = lwr, y = med, ymax = upr, 
        fill = modis_color), shape = 21, stroke = 1, size = 0.8, 
        data = m_model$eta2,  position = position_nudge(x = -0.35)) +
    geom_pointrange(aes(x = IGBP, ymin = lwr, y = med, ymax = upr, 
        color = flux_color), shape = 16, alpha = 0.3, size = 0.3, 
        data = f_model$beta2, position = position_dodge2(width = 0.5)) +
    geom_pointrange(aes(x = IGBP, ymin = lwr, y = med, ymax = upr, 
        fill = flux_color), shape = 21, size = 0.8, data = f_model$eta2, 
        position = position_nudge(x = -0.45)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
    ylab(expression("GSL Coefficient" ~ (g ~ C ~ m^"-2"))) +
    geom_text(aes(x = IGBP, y = rep(-20, 11), label = N), 
        data = north_sites_dt[, .N, by = "IGBP"]) +
    geom_text(aes(x = IGBP, y = rep(-15, 11), label = N), 
        data = unique(north_sites_dt[, .(site, IGBP)])[, .N, by = c("IGBP")]) +
    theme_article() +
    scale_y_continuous(breaks = seq(-40, 40, by = 10)) +
    scale_color_identity() +
    scale_fill_manual(name = "", values = c("#E7298A" = "#E7298A", 
        "#7570B3" = "#7570B3"), labels = c("MODIS EVI2", "FLUX GPP"), 
        guide = "legend") +
    theme(legend.position = c(0.9, 0.93)) +
    guides(fill = guide_legend(
        keywidth = 0.1,
        keyheight = 0.3,
        default.unit = "inch"
    ))

gsl_coef_fig

# fig: model coef for all biomes
ggsave(filename = "Output/model_coef_all_biomes.png", plot = gsl_coef_fig, 
    width = 10, height = 4, bg = "white")
