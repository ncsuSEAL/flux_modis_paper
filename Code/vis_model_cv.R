#************************************************************************************
# Description: Cross validation result figure for models
# Author: Xiaojie(J) Gao
# Date: 2022-03-04
#************************************************************************************
rm(list=ls())

source("Code/base.R")
source("Code/mod_base.R")

color_pal <- brewer.pal(8, "Dark2")
modis_color <- color_pal[4]
flux_color <- color_pal[3]


load(file.path("Pipeline/cv.RData"))

uni_yrs <- sort(unique(north_sites_dt$year))

{
    png("Output/mod_cv.png", width = 1000, height = 1000, res = 150)
    # png("0_Manuscript/SI/SI_Fig_5.png", width = 1000, height = 1000, res = 150)

    par(mar = c(3, 4, 2, 2))
    plot(uni_yrs, sqrt(f_cv_mse_1),
        type = "b", ylim = c(0, 500), col = flux_color, pch = 16,
        xlab = "Year", ylab = expression("RMSE" ~ (g ~ C ~ m^"-2")), 
        mgp = c(1.5, 0.5, 0), lwd = 1
    )
    lines(uni_yrs, sqrt(m_cv_mse_1), type = "b", col = modis_color, pch = 16, lwd = 1)
    lines(uni_yrs, sqrt(f_cv_mse_2), type = "b", col = flux_color, pch = 22, lwd = 1)
    lines(uni_yrs, sqrt(m_cv_mse_2), type = "b", col = modis_color, pch = 22, lwd = 1)
    lines(uni_yrs, sqrt(f_cv_mse_3), type = "b", col = flux_color, pch = 2, lwd = 1)
    lines(uni_yrs, sqrt(m_cv_mse_3), type = "b", col = modis_color, pch = 2, lwd = 1)
    lines(uni_yrs, sqrt(f_cv_mse_4), type = "b", col = flux_color, pch = 4, lwd = 1)
    lines(uni_yrs, sqrt(m_cv_mse_4), type = "b", col = modis_color, pch = 4, lwd = 1)
    lines(uni_yrs, sqrt(f_cv_mse_5), type = "b", col = flux_color, pch = 7, lwd = 1)
    lines(uni_yrs, sqrt(m_cv_mse_5), type = "b", col = modis_color, pch = 7, lwd = 1)

    legend("bottomright",
        ncol = 3, legend = c("Model 1", "Model 2", "Model 3", "Model 4", "Model 5", 
        "GPP based", "EVI2 based"), pch = c(16, 22, 2, 4, 7, NA, NA), lwd = 1, 
        col = c(rep("black", 5), flux_color, modis_color), bty = "n"
    )
    dev.off()
}
