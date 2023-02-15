#************************************************************************************
# Description: GPPmax vs. EVI2max relationship.
# Author: Xiaojie(J) Gao
# Date: 2022-03-03
#************************************************************************************
rm(list = ls())

source("Code/base.R")
source("Code/mod_base.R")


color_pal <- brewer.pal(8, "Dark2")

ever_col <- color_pal[1]
deci_col <- color_pal[2]



# ~ Log transformation plot on the original scale
# ~~~~~~~~~~~~~~~~~
# compute site means across years
north_sites_dt[, ":="(m_EVImax_1_mean = mean(m_EVImax_1, na.rm = TRUE),
    f_gppmax_1_mean = mean(f_gppmax_1, na.rm = TRUE)), by = "site"]

{
    # png("0_Manuscript/SI/SI_Fig_6.png", width = 800, height = 800, res = 170)
    png("Output/gppmax_evi2max.png", width = 800, height = 800, res = 170)

    par(mar = c(3, 3, 1, 1))

    colorsTable <- brewer.pal(5, "Dark2")
    plot(north_sites_dt[cat == "ever", .(m_EVImax_1_mean, f_gppmax_1_mean)],
        mgp = c(1.5, 0.5, 0), xlab = "EVI2max", ylab = "GPPmax", pch = 16,
        col = ever_col, xlim = c(0.1, 0.9), ylim = c(2, 30)
    )
    fit_ever <- lm(
        log(f_gppmax_1_mean) ~ m_EVImax_1_mean, 
        data = north_sites_dt[cat == "ever"]
    )
    # abline(fit_ever, col = ever_col, lwd = 2)
    newx <- seq(grconvertX(0, from = "nfc", to = "user"), 
        grconvertX(1, from = "nfc", to = "user"), length.out = 100)
    preds_ever <- predict(fit_ever, newdata = data.frame(m_EVImax_1_mean = newx), 
        interval = "confidence")
    lines(newx, exp(preds_ever[, 1]), lwd = 2, col = ever_col)
    polygon(c(rev(newx), newx), c(rev(exp(preds_ever[, 3])), exp(preds_ever[, 2])), 
        col = Transparent(ever_col, 0.3), border = NA)
    text(0.1, 20, col = ever_col, labels = c(
        paste0("y = ", "exp(", round(coef(fit_ever)[2], 2), "x +", 
            round(coef(fit_ever)[1], 2), ")"),
        paste0("\n\nR2:", round(summary(fit_ever)$r.squared, 2))
    ), adj = 0)

    points(north_sites_dt[cat == "deci", .(m_EVImax_1_mean, f_gppmax_1_mean)], 
        col = deci_col, pch = 16)
    fit_deci <- lm(log(f_gppmax_1_mean) ~ m_EVImax_1_mean, 
        data = north_sites_dt[cat == "deci"])
    # abline(fit_deci, col = deci_col, lwd = 2)
    newx <- seq(grconvertX(0, from = "nfc", to = "user"), 
        grconvertX(1, from = "nfc", to = "user"), length.out = 100)
    preds_deci <- predict(fit_deci, newdata = data.frame(m_EVImax_1_mean = newx), 
        interval = "confidence")
    lines(newx, exp(preds_deci[, 1]), lwd = 2, col = deci_col)
    polygon(c(rev(newx), newx), c(rev(exp(preds_deci[, 3])), exp(preds_deci[, 2])), 
        col = Transparent(deci_col, 0.3), border = NA)
    text(0.1, 15, col = deci_col, labels = c(
        paste0("y = ", "exp(", round(coef(fit_deci)[2], 2), "x +", 
            round(coef(fit_deci)[1], 2), ")"),
        paste0("\n\nR2:", round(summary(fit_deci)$r.squared, 2))
    ), adj = 0)

    legend("topleft", bty = "n", col = colorsTable[1:2], pch = 16, 
        legend = c("Evergreen", "Deciduous"))

    dev.off()
}