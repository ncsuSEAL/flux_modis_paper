#************************************************************************************
# Description: MODIS growing season length vs. Flux growing season length.
#              The GSL are site means across years.
# Author: Xiaojie(J) Gao
# Date: 2022-03-04
#************************************************************************************
rm(list = ls())

source("Code/base.R")
source("Code/mod_base.R")


color_pal <- brewer.pal(8, "Dark2")

modis_color <- color_pal[4]
flux_color <- color_pal[3]

ever_col <- color_pal[1]
deci_col <- color_pal[2]

# compute site means across years
north_sites_dt[, ":="(m_gsl_mean = mean(m_gsl, na.rm = TRUE),
    f_gsl_mean = mean(f_gsl, na.rm = TRUE)), by = "site"]

{
    # png("0_Manuscript/SI/SI_Fig_gsl_evi_gpp.png", width = 1000, height = 1000, 
    #     res = 170)
    png("Output/modis_gsl_flux_gsl.png", width = 1000, height = 1000, res = 170)

    colorsTable <- brewer.pal(5, "Dark2")
    plot(north_sites_dt[cat == "ever", .(m_gsl_mean, f_gsl_mean)],
        mgp = c(1.5, 0.5, 0), xlab = "EVI2 GSL", ylab = "GPP GSL", pch = 16,
        col = ever_col, xlim = c(50, 250), ylim = c(50, 250)
    )
    abline(0, 1, lty = 2)
    fit_ever <- lm(f_gsl_mean ~ m_gsl_mean, data = north_sites_dt[cat == "ever"])
    # abline(fit_ever, col = ever_col, lwd = 2)
    newx <- seq(grconvertX(0, from = "nfc", to = "user"), 
        grconvertX(1, from = "nfc", to = "user"), length.out = 100)
    preds_ever <- predict(fit_ever, newdata = data.frame(m_gsl_mean = newx), 
        interval = "confidence")
    lines(newx, preds_ever[, 1], lwd = 2, col = ever_col)
    polygon(c(rev(newx), newx), c(rev(preds_ever[, 3]), preds_ever[, 2]), 
        col = Transparent(ever_col, 0.3), border = NA)
    text(50, 220, col = ever_col, labels = c(
        paste0("y = ", round(coef(fit_ever)[2], 2), "x +", 
            round(coef(fit_ever)[1], 2)),
        paste0("\n\nR2:", round(summary(fit_ever)$r.squared, 2))
    ), adj = 0)

    points(north_sites_dt[cat == "deci", .(m_gsl_mean, f_gsl_mean)], col = deci_col, 
        pch = 16)
    fit_deci <- lm(f_gsl_mean ~ m_gsl_mean, data = north_sites_dt[cat == "deci"])
    # abline(fit_deci, col = deci_col, lwd = 2)
    newx <- seq(grconvertX(0, from = "nfc", to = "user"), 
        grconvertX(1, from = "nfc", to = "user"), length.out = 100)
    preds_deci <- predict(fit_deci, newdata = data.frame(m_gsl_mean = newx), 
        interval = "confidence")
    lines(newx, preds_deci[, 1], lwd = 2, col = deci_col)
    polygon(c(rev(newx), newx), c(rev(preds_deci[, 3]), preds_deci[, 2]), 
        col = Transparent(deci_col, 0.3), border = NA)
    text(50, 200, col = deci_col, labels = c(
        paste0("y = ", round(coef(fit_deci)[2], 2), "x +", 
            round(coef(fit_deci)[1], 2)),
        paste0("\n\nR2:", round(summary(fit_deci)$r.squared, 2))
    ), adj = 0)

    legend("topleft", bty = "n", col = colorsTable[1:2], pch = 16, 
        legend = c("Evergreen", "Deciduous"))

    dev.off()
}