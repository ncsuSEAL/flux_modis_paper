#*******************************************************************************
# Description: Flux and MODIS phenometrics comparison, Midgup and Midgdown only.
# 
# Author: Xiaojie Gao
# Date: 2022-03-03
#*******************************************************************************
rm(list = ls())

source("Code/base.R")
source("Code/mod_base.R")


color_pal <- brewer.pal(8, "Dark2")

ever_col <- color_pal[1]
deci_col <- color_pal[2]


# The comparion result
comp_result_all <- fread(file.path(gdir, "comparison_new.csv"))
comp_result_all$V1 <- NULL
# remove site year that flux are 0
comp_result <- comp_result_all[f_gup_1 != 0 | f_gup_2 != 0, ]

# I want to try melt and dcast here
comp_result_melt <- melt(comp_result[, c(1:30, 43), with = FALSE], 
    id = c("site", "year", "IGBP"))
comp_result_melt[, c("source", "pheno", "cycle") := tstrsplit(variable, "_")]
comp_result_dcast <- dcast(comp_result_melt, 
    site + year + pheno + cycle + IGBP ~ source, value.var = "value")

# convert to DOY
comp_result_doy <- na.omit(comp_result_dcast)[, .(
    IGBP, site, pheno, cycle, year,
    f = f - as.integer(as.Date(paste0(year, "-1-1")) - as.Date("1970-1-1")),
    m = m - as.integer(as.Date(paste0(year, "-1-1")) - as.Date("1970-1-1"))
)]

midgup_ever <- comp_result_doy[
    pheno == "midgup" & IGBP %in% c("ENF", "EBF"), 
]
midgup_deci <- comp_result_doy[
    pheno == "midgup" & !(IGBP %in% c("ENF", "EBF")), 
]

midgdown_ever <- comp_result_doy[
    pheno == "midgdown" & IGBP %in% c("ENF", "EBF"), 
]
midgdown_deci <- comp_result_doy[
    pheno == "midgdown" & !(IGBP %in% c("ENF", "EBF")),
]




{
    # png("0_Manuscript/Fig_3.png", width = 2100, height = 1200, res = 300)
    png("Output/f_m_pheno_comp_mid.png", width = 2100, height = 1200, res = 300)

    par(mgp = c(1.5, 0.5, 0), mfrow = c(1, 2), mar = c(3, 2.5, 0.8, 0.3))
    # ~ Left panel
    # ~~~~~~~~~~~~~~~~~
    # deciduous midgup
    plot(midgup_deci[, .(m, f)],
        cex = 0.5, lab = c(5, 5, 0),
        xlab = NA, ylab = NA, pch = 16, col = Transparent(deci_col, 0.3),
        xlim = range(range(midgup_deci$m), range(midgup_deci$f)), 
        ylim = range(range(midgup_deci$m), range(midgup_deci$f))
    )
    mtext("MODIS (DOY)", side = 1, line = 1.5, cex = 1)
    mtext("Flux (DOY)", side = 2, line = 1.5, cex = 1)
    abline(a = 0, b = 1, lty = 2, xpd = FALSE) # 1:1 line
    fit_deci <- lm(f ~ m, data = midgup_deci)
    abline(fit_deci, col = deci_col, lwd = 2)
    newx <- seq(grconvertX(0, from = "nfc", to = "user"), 
        grconvertX(1, from = "nfc", to = "user"), length.out = 100)
    preds_deci <- predict(fit_deci, newdata = data.frame(m = newx), 
        interval = "confidence")
    polygon(c(rev(newx), newx), c(rev(preds_deci[, 3]), preds_deci[, 2]), 
        col = Transparent(deci_col, 0.3), border = NA)

    # evergreen midgup
    points(midgup_ever[, .(m, f)],
        cex = 0.5,
        xlab = NA, ylab = NA, pch = 16, col = Transparent(ever_col, 0.3)
    )
    fit_ever <- lm(f ~ m, data = midgup_ever)
    abline(fit_ever, col = ever_col, lwd = 2, xpd = FALSE)
    preds_ever <- predict(fit_ever, newdata = data.frame(m = newx), 
        interval = "confidence")
    polygon(c(rev(newx), newx), c(rev(preds_ever[, 3]), preds_ever[, 2]), 
        col = Transparent(ever_col, 0.3), border = NA)

    legend("topleft", bty = "n", xpd = TRUE, 
        legend = c("Evergreen", "Deciduous"), 
        pch = 16, col = c(ever_col, deci_col)
    )

    # ~ Right panel
    # ~~~~~~~~~~~~~~~~~
    # deciduous midgdown
    plot(midgdown_deci[, .(m, f)],
        cex = 0.5, lab = c(5, 5, 0),
        xlab = NA, ylab = NA, pch = 16, col = Transparent(deci_col, 0.3),
        xlim = range(range(midgdown_deci$m), range(midgdown_deci$f)),
        ylim = range(range(midgdown_deci$m), range(midgdown_deci$f))
    )
    mtext("MODIS (DOY)", side = 1, line = 1.5, cex = 1)
    mtext("Flux (DOY)", side = 2, line = 1.5, cex = 1)
    abline(a = 0, b = 1, lty = 2, xpd = FALSE) # 1:1 line
    fit_deci <- lm(f ~ m, data = midgdown_deci)
    abline(fit_deci, col = deci_col, lwd = 2)
    newx <- seq(grconvertX(0, from = "nfc", to = "user"), 
        grconvertX(1, from = "nfc", to = "user"), length.out = 100)
    preds_deci <- predict(fit_deci, newdata = data.frame(m = newx), 
        interval = "confidence")
    polygon(c(rev(newx), newx), c(rev(preds_deci[, 3]), preds_deci[, 2]), 
        col = Transparent(deci_col, 0.3), border = NA)

    # evergreen midgdown
    points(midgdown_ever[, .(m, f)],
        cex = 0.5,
        xlab = NA, ylab = NA, pch = 16, col = Transparent(ever_col, 0.3)
    )
    fit_ever <- lm(f ~ m, data = midgdown_ever)
    abline(fit_ever, col = ever_col, lwd = 2, xpd = FALSE)
    preds_ever <- predict(fit_ever, newdata = data.frame(m = newx),  
        interval = "confidence")
    polygon(c(rev(newx), newx), c(rev(preds_ever[, 3]), preds_ever[, 2]),
        col = Transparent(ever_col, 0.3), border = NA)

    legend("topleft", bty = "n", xpd = TRUE, 
        legend = c("Evergreen", "Deciduous"),
        pch = 16, col = c(ever_col, deci_col)
    )


    # ~ boxplot subpanel
    # ~~~~~~~~~~~~~~~~~
    par(fig = c(0.35, 0.48, 0.2, 0.5), new = TRUE)
    par(mar = c(0, 0, 0, 0), mgp = c(1.5, 0.5, 0))
    boxplot(c(midgup_ever[, .(f, m)], midgup_deci[, .(f, m)]), xaxt = "n", 
        yaxt = "n", col = rep(c(ever_col, deci_col), each = 2), 
        outline = FALSE, boxwex = 0.5
    )
    axis(side = 1, at = 1:4, padj = -1, labels = c("F", "M", "F", "M"), 
        cex.axis = 0.8, tck = -0.03, gap.axis = 0.05)
    axis(side = 2, at = seq(50, 250, length.out = 5), labels = seq(50, 250, 
        length.out = 5), cex.axis = 0.8, tck = -0.03, gap.axis = 0.05, 
        las = 2, hadj = 0.7
    )

    par(fig = c(0.85, 0.98, 0.2, 0.5), new = TRUE)
    par(mar = c(0, 0, 0, 0), mgp = c(1.5, 0.5, 0))
    boxplot(c(midgdown_ever[, .(f, m)], midgdown_deci[, .(f, m)]), xaxt = "n", 
        yaxt = "n", col = rep(c(ever_col, deci_col), each = 2), 
        outline = FALSE, boxwex = 0.5
    )
    axis(side = 1, at = 1:4, padj = -1, labels = c("F", "M", "F", "M"), 
        cex.axis = 0.8, tck = -0.03, gap.axis = 0.05)
    axis(side = 2, at = seq(150, 350, length.out = 5), labels = seq(150, 350, 
        length.out = 5), cex.axis = 0.8, tck = -0.03, gap.axis = 0.05, 
        las = 2, hadj = 0.7
    )


    par(fig = c(0, 1, 0, 1), new = TRUE, xpd = NA)
    text(grconvertX(0, "nfc", "user"), grconvertY(0.95, "nfc", "user"), 
        labels = "a", pos = 4, font = 2)
    text(grconvertX(0.5, "nfc", "user"), grconvertY(0.95, "nfc", "user"), 
        labels = "b", pos = 4, font = 2)

    dev.off()
}