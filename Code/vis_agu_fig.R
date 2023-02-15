#************************************************************************
# Description: Some figures for the AGU 2021 fall meeting 
# Author: Xiaojie(J) Gao
# Date: 2022-01-01
#************************************************************************

source("Code/base.R")
source("Code/mod_base.R")
library(ggplot2)
library(patchwork)
library(egg)
# library(ggpubr)
library(gridExtra)

focus_biomes <- c("CRO", "DBF", "ENF", "GRA", "MF", "WET")

color_pal <- brewer.pal(8, "Dark2")

modis_color <- color_pal[4]
flux_color <- color_pal[3]


#region [Fig 1 baseplot version] ####
#************************************************************
#! Data comes from "bayesian_model.R"
load(file.path("Pipeline/best_model_fit.RData"))
load(file.path("Pipeline/model_3_fit.RData"))


igbp_colors <- c("#51B6F5", "#218A21", "#31CD31", "#9ACD31", "#97FA97", "#8FBB8F", "#BB8F8F", 
    "#F5DEB3", "#DBEB9D", "#FFD600", "#EFB766", "#4682B2", "#FAED73", "#FF0000", "#999355", "#F5F5DC", "#BDBDBD", "#000000")
# igbp_names <- c("water", "enf", "ebf", "dnf", "dbf", "mixed", "closed shrubs", "open shrubs", "woody savannas", "savannas", 
    # "grasslands", "perm wetlands", "croplands", "urban", "crop/natural mosaic", "snow and ice", "barren/sparse veg", "unclassified")
igbp_names <- c("water", "ENF", "EBF", "DNF", "DBF", "MF", "CSH", "OSH", "WSA", "SAV", 
    "GRA", "WET", "CRO", "urban", "crop/natural mosaic", "snow and ice", "barren/sparse veg", "unclassified")

IGBP <- unique(as.character(f_igbp_factor))
# cols <- RColorBrewer::brewer.pal(11, "Paired") # Set these to be static
colsIGBP <- merge(data.frame(IGBP = IGBP), data.frame(IGBP = igbp_names, cols = igbp_colors))
colsIGBP$IGBP <- as.character(colsIGBP$IGBP)
colsIGBP$cols <- as.character(colsIGBP$cols)

f_dt_model_3 <- data.table(fit = f_model_3$fit * 100, agpp = Y * 100, IGBP = f_igbp_factor)
f_dt_model_3 <- merge(f_dt_model_3, colsIGBP)
f_dt_model_3[, ":="(IGBP = as.character(IGBP), cols = as.character(cols))]

m_dt_model_3 <- data.table(fit = m_model_3$fit * 100, agpp = Y * 100, IGBP = f_igbp_factor)
m_dt_model_3 <- merge(m_dt_model_3, colsIGBP)
m_dt_model_3[, ":="(IGBP = as.character(IGBP), cols = as.character(cols))]

f_dt <- data.table(fit = f_model$fit * 100, agpp = Y * 100, IGBP = f_igbp_factor)
f_dt <- merge(f_dt, colsIGBP)
f_dt[, ":="(IGBP = as.character(IGBP), cols = as.character(cols))]

m_dt <- data.table(fit = m_model$fit * 100, agpp = Y * 100, IGBP = f_igbp_factor)
m_dt <- merge(m_dt, colsIGBP)
m_dt[, ":="(IGBP = as.character(IGBP), cols = as.character(cols))]


#Fig: Fig 1
{
    png("Output/fig1.png", width = 920, height = 800, res = 150)
    layout(matrix(c(1, 2, 3, 4), nrow = 2))
    par(
        mgp = c(2, 0.5, 0), cex.lab = 1.3, cex.axis = 1.2, oma = c(5, 0, 0, 0),
        bg = NA, fg = "white", col.axis = "white"
    )

    # ~ 1a
    # ~~~~~~~~~~~~~~~~~
    par(mar = c(0, 4, 2.5, 0))
    # GPP model 3 fit
    plot(f_dt_model_3[, .(fit.fit_median, agpp)],
        bg = f_dt_model_3$cols, pch = 21, lwd = 0,
        xlim = range(f_dt$fit.fit_median, f_dt$agpp),
        ylim = range(f_dt$fit.fit_median, f_dt$agpp),
        xlab = "", ylab = "",
        xaxt = "n", yaxt = "n"
    )
    # axis(side = 1, at = c(1000, 2000, 3000), labels = c(1000, 2000, 3000))
    axis(side = 2, at = c(1000, 2000, 3000), labels = c(1000, 2000, 3000))
    # mtext(side = 2, expression("True ΣGPP" ~ (g ~ C ~ m^"-2" ~ yr^"-1")), line = 1.5)
    text(1800, 900, cex = 1.3, labels = c(
        paste0("R2:", round(cor(f_dt_model_3$fit.fit_median, f_dt_model_3$agpp)^2, 2)),
        paste0("\n\nRMSE:", round(sqrt(mean((f_dt_model_3$fit.fit_median - f_dt_model_3$agpp)^2)), 2))
    ), adj = 0)
    abline(0, 1, lty = 2)

    # ~ 1b
    # ~~~~~~~~~~~~~~~~~
    par(mar = c(3, 4, 0, 0))
    # GPP best model fit
    plot(f_dt[, .(fit.fit_median, agpp)],
        bg = f_dt$cols, pch = 21, lwd = 0,
        xlim = range(f_dt$fit.fit_median, f_dt$agpp),
        ylim = range(f_dt$fit.fit_median, f_dt$agpp),
        xlab = "", ylab = "",
        xaxt = "n", yaxt = "n"
    )
    axis(side = 1, at = c(1000, 2000, 3000), labels = c(1000, 2000, 3000))
    axis(side = 2, at = c(1000, 2000, 3000), labels = c(1000, 2000, 3000))
    # mtext(side = 1, expression("Fitted ΣGPP" ~ (g ~ C ~ m^"-2"~ yr^"-1")), line = 2.5)
    # mtext(side = 2, expression("True ΣGPP" ~ (g ~ C ~ m^"-2"~ yr^"-1")), line = 1.5)
    text(1800, 900, cex = 1.3, labels = c(
        paste0("R2:", round(cor(f_dt$fit.fit_median, f_dt$agpp)^2, 2)),
        paste0("\n\nRMSE:", round(sqrt(mean((f_dt$fit.fit_median - f_dt$agpp)^2)), 2))
    ), adj = 0)
    abline(0, 1, lty = 2)


    # ~ 1c
    # ~~~~~~~~~~~~~~~~~
    par(mar = c(0, 0, 2.5, 3.5))
    # EVI2 model 3 fit
    plot(m_dt_model_3[, .(fit.fit_median, agpp)],
        bg = m_dt_model_3$cols, pch = 21, lwd = 0,
        xlim = range(m_dt_model_3$fit.fit_median, m_dt_model_3$agpp),
        ylim = range(m_dt_model_3$fit.fit_median, m_dt_model_3$agpp),
        xlab = "", ylab = "",
        xaxt = "n", yaxt = "n"
    )
    # axis(side = 1, at = c(1000, 2000, 3000), labels = c(1000, 2000, 3000))
    # axis(side = 2, at = c(1000, 2000, 3000), labels = c(1000, 2000, 3000))
    text(1800, 900, cex = 1.3, labels = c(
        paste0("R2:", round(cor(m_dt_model_3$fit.fit_median, m_dt_model_3$agpp)^2, 2)),
        paste0("\n\nRMSE:", round(sqrt(mean((m_dt_model_3$fit.fit_median - m_dt_model_3$agpp)^2)), 2))
    ), adj = 0)
    abline(0, 1, lty = 2)

    # ~ 1d
    # ~~~~~~~~~~~~~~~~~
    par(mar = c(3, 0, 0, 3.5))
    # EVI2 best model fit
    plot(m_dt[, .(fit.fit_median, agpp)],
        bg = m_dt$cols, pch = 21, lwd = 0,
        xlim = range(m_dt$fit.fit_median, m_dt$agpp),
        ylim = range(m_dt$fit.fit_median, m_dt$agpp),
        xlab = "", ylab = "",
        xaxt = "n", yaxt = "n"
    )
    axis(side = 1, at = c(1000, 2000, 3000), labels = c(1000, 2000, 3000))
    # mtext(side = 1, expression("Fitted ΣGPP" ~ (g ~ C ~ m^"-2"~ yr^"-1")), line = 2.5)
    # axis(side = 2, at = c(1000, 2000, 3000), labels = c(1000, 2000, 3000))
    text(1800, 900, cex = 1.3, labels = c(
        paste0("R2:", round(cor(m_dt$fit.fit_median, m_dt$agpp)^2, 2)),
        paste0("\n\nRMSE:", round(sqrt(mean((m_dt$fit.fit_median - m_dt$agpp)^2)), 2))
    ), adj = 0)
    abline(0, 1, lty = 2)

    text(grconvertX(0.04, "ndc", "user"), grconvertY(0.55, "ndc", "user"),
        expression("True annual GPP" ~ (g ~ C ~ m^"-2" ~ yr^"-1")),
        xpd = NA, srt = 90, cex = 1.3
    )
    text(grconvertX(0.5, "ndc", "user"), grconvertY(0.17, "ndc", "user"),
        expression("Fitted annual GPP" ~ (g ~ C ~ m^"-2" ~ yr^"-1")),
        xpd = NA, cex = 1.3
    )

    # ~ Legend
    # ~~~~~~~~~~~~~~~~~
    legend(grconvertX(0.5, "ndc", "user"), grconvertY(0.13, "ndc", "user"),
        pch = 21, pt.lwd = 0.3, pt.bg = as.character(colsIGBP$cols),
        legend = colsIGBP$IGBP,
        xpd = NA, ncol = 6, bty = "n", cex = 1.3, xjust = 0.5
    )

    dev.off()
}
#endregion [Fig 1 baseplot version]


#region [Normalized model] ####
#************************************************************
#! Data comes from "bayesian_model.R"
load(file.path("Pipeline/best_model_fit.RData"))
load(file.path("Pipeline/model_3_fit.RData"))


igbp_colors <- c(
    "#51B6F5", "#218A21", "#31CD31", "#9ACD31", "#97FA97", "#8FBB8F", "#BB8F8F",
    "#F5DEB3", "#DBEB9D", "#FFD600", "#EFB766", "#4682B2", "#FAED73", "#FF0000", "#999355", "#F5F5DC", "#BDBDBD", "#000000"
)
# igbp_names <- c("water", "enf", "ebf", "dnf", "dbf", "mixed", "closed shrubs", "open shrubs", "woody savannas", "savannas",
# "grasslands", "perm wetlands", "croplands", "urban", "crop/natural mosaic", "snow and ice", "barren/sparse veg", "unclassified")
igbp_names <- c(
    "water", "ENF", "EBF", "DNF", "DBF", "MF", "CSH", "OSH", "WSA", "SAV",
    "GRA", "WET", "CRO", "urban", "crop/natural mosaic", "snow and ice", "barren/sparse veg", "unclassified"
)

IGBP <- unique(as.character(f_igbp_factor))
# cols <- RColorBrewer::brewer.pal(11, "Paired") # Set these to be static
colsIGBP <- merge(data.frame(IGBP = IGBP), data.frame(IGBP = igbp_names, cols = igbp_colors))
colsIGBP$IGBP <- as.character(colsIGBP$IGBP)
colsIGBP$cols <- as.character(colsIGBP$cols)

f_dt_model_3 <- data.table(fit = f_model_3$fit * 100, agpp = Y * 100, IGBP = f_igbp_factor)
f_dt_model_3 <- merge(f_dt_model_3, colsIGBP)
f_dt_model_3[, ":="(IGBP = as.character(IGBP), cols = as.character(cols))]

m_dt_model_3 <- data.table(fit = m_model_3$fit * 100, agpp = Y * 100, IGBP = f_igbp_factor)
m_dt_model_3 <- merge(m_dt_model_3, colsIGBP)
m_dt_model_3[, ":="(IGBP = as.character(IGBP), cols = as.character(cols))]

f_dt <- data.table(fit = f_model$fit * 100, agpp = Y * 100, IGBP = f_igbp_factor)
f_dt <- merge(f_dt, colsIGBP)
f_dt[, ":="(IGBP = as.character(IGBP), cols = as.character(cols))]

m_dt <- data.table(fit = m_model$fit * 100, agpp = Y * 100, IGBP = f_igbp_factor)
m_dt <- merge(m_dt, colsIGBP)
m_dt[, ":="(IGBP = as.character(IGBP), cols = as.character(cols))]


{
    png("Output/normal_effect.png", width = 800, height = 1100, res = 150)
    opar <- par()
    par(mar = c(4, 7, 3, 1), bg = NA, fg = "white", col.axis = "white")
    # normalized effect
    focus_biomes <- c("CRO", "DBF", "ENF", "GRA", "MF", "WET")

    plot(NA, xlim = c(-0.4, 1.4), ylim = c(0.5, 6.7), xlab = "", ylab = "", yaxt = "n", mgp = c(1.5, 0.5, 0))
    axis(side = 2, at = 1:6, labels = c("EVI2 min", "GPP min", "EVI2 max", "GPP max", "EVI2 GSL", "GPP GSL"), las = 2, cex.axis = 1.3)
    mtext(side = 1, "Normalized Effect", line = 2, cex = 1.3)
    abline(v = 0, lty = 2, col = "grey80")

    EtaPointRange <- function(eta_df, colsIGBP_df, focus_biomes, ycoord) {
        tmp_df <- merge(eta_df, colsIGBP_df)
        tmp_df$cols <- as.character(tmp_df$cols)
        points(tmp_df[IGBP %in% focus_biomes, med / 100], ycoord, pch = 21, lwd = 0.3, bg = tmp_df[IGBP %in% focus_biomes, cols], cex = 1.3)
        segments(
            x0 = tmp_df[IGBP %in% focus_biomes, lwr / 100], y0 = ycoord,
            x1 = tmp_df[IGBP %in% focus_biomes, upr / 100], y1 = ycoord, col = tmp_df[IGBP %in% focus_biomes, cols],
            lwd = 1.3
        )
    }

    EtaPointRange(f_model_norm$eta2, colsIGBP, focus_biomes, seq(6.1, 6.4, length = 6) - 0.2)
    EtaPointRange(f_model_norm$eta3, colsIGBP, focus_biomes, seq(4.1, 4.4, length = 6) - 0.2)
    EtaPointRange(f_model_norm$eta4, colsIGBP, focus_biomes, seq(2.1, 2.4, length = 6) - 0.2)
    ycoord <- c(2, 4, 6) - 0.25
    segments(x0 = f_model_norm$mu[4:2, lwr / 100], y0 = ycoord, x1 = f_model_norm$mu[4:2, upr / 100], y1 = ycoord, "white", lwd = 3)
    points(f_model_norm$mu[4:2, med / 100], ycoord, pch = 21, bg = flux_color, col = 1, cex = 3, lwd = 1.5)

    EtaPointRange(m_model_norm$eta2, colsIGBP, focus_biomes, seq(5.1, 5.4, length = 6) - 0.2)
    EtaPointRange(m_model_norm$eta3, colsIGBP, focus_biomes, seq(3.1, 3.4, length = 6) - 0.2)
    EtaPointRange(m_model_norm$eta4, colsIGBP, focus_biomes, seq(1.1, 1.4, length = 6) - 0.2)
    ycoord <- c(1, 3, 5) + 0.3
    segments(x0 = m_model_norm$mu[4:2, lwr / 100], y0 = ycoord, x1 = m_model_norm$mu[4:2, upr / 100], y1 = ycoord, col = "white", lwd = 3)
    points(m_model_norm$mu[4:2, med / 100], ycoord, pch = 21, bg = modis_color, col = 1, cex = 3, lwd = 1.5)
    legend("bottomright",
        legend = c("MODIS EVI2", "Flux GPP"),
        pch = 21, pt.bg = c(modis_color, flux_color), pt.cex = 2,
        pt.lwd = 1.5, col = 1, bty = "n", cex = 1
    )

    # legend(grconvertX(0.5, "ndc", "user"), grconvertY(0.95, "ndc", "user"),
    #     pch = 21, pt.lwd = 0.3, pt.bg = as.character(colsIGBP$cols),
    #     legend = colsIGBP$IGBP,
    #     xpd = NA, ncol = 6, bty = "n", cex = 1.3, xjust = 0.5)
    legend(grconvertX(0.56, "ndc", "user"), grconvertY(1, "ndc", "user"),
        pch = 21, pt.lwd = 0.3, pt.bg = as.character(colsIGBP$cols),
        legend = colsIGBP$IGBP,
        xpd = NA, ncol = 6, bty = "n", cex = 1, xjust = 0.5
    )

    dev.off()
}
#endregion [Normalized model]



#region [Pheno compare] ####
#************************************************************
# The comparion result
comp_result_all <- fread("Data/comparison_new.csv")
comp_result_all$V1 <- NULL
# remove site year that flux are 0
comp_result <- comp_result_all[f_gup_1 != 0 | f_gup_2 != 0, ]

# I want to try melt and dcast here
comp_result_melt <- melt(comp_result[, c(1:30, 43), with = FALSE], id = c("site", "year", "IGBP"))
comp_result_melt[, c("source", "pheno", "cycle") := tstrsplit(variable, "_")]
comp_result_dcast <- dcast(comp_result_melt, site + year + pheno + cycle + IGBP ~ source, value.var = "value")

# convert to DOY
comp_result_doy <- na.omit(comp_result_dcast)[, .(IGBP, site, pheno, cycle, year,
    f = f - as.integer(as.Date(paste0(year, "-1-1")) - as.Date("1970-1-1")),
    m = m - as.integer(as.Date(paste0(year, "-1-1")) - as.Date("1970-1-1"))
)]

midgup_ever <- comp_result_doy[pheno == "midgup" & IGBP %in% c("ENF", "EBF"), ]
midgup_deci <- comp_result_doy[pheno == "midgup" & !(IGBP %in% c("ENF", "EBF")), ]

midgdown_ever <- comp_result_doy[pheno == "midgdown" & IGBP %in% c("ENF", "EBF"), ]
midgdown_deci <- comp_result_doy[pheno == "midgdown" & !(IGBP %in% c("ENF", "EBF")), ]


ever_col <- color_pal[1]
deci_col <- color_pal[2]

{
    opar <- par()
    # png("0_Manuscript/Fig_3.png", width = 2100, height = 1200, res = 300)
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
    newx <- seq(grconvertX(0, from = "nfc", to = "user"), grconvertX(1, from = "nfc", to = "user"), length.out = 100)
    preds_deci <- predict(fit_deci, newdata = data.frame(m = newx), interval = "confidence")
    polygon(c(rev(newx), newx), c(rev(preds_deci[, 3]), preds_deci[, 2]), col = Transparent(deci_col, 0.3), border = NA)

    # evergreen midgup
    points(midgup_ever[, .(m, f)],
        cex = 0.5,
        xlab = NA, ylab = NA, pch = 16, col = Transparent(ever_col, 0.3)
    )
    fit_ever <- lm(f ~ m, data = midgup_ever)
    abline(fit_ever, col = ever_col, lwd = 2, xpd = FALSE)
    preds_ever <- predict(fit_ever, newdata = data.frame(m = newx), interval = "confidence")
    polygon(c(rev(newx), newx), c(rev(preds_ever[, 3]), preds_ever[, 2]), col = Transparent(ever_col, 0.3), border = NA)

    legend("topleft", bty = "n", xpd = TRUE, legend = c("Evergreen", "Deciduous"), pch = 16, col = c(ever_col, deci_col))

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
    newx <- seq(grconvertX(0, from = "nfc", to = "user"), grconvertX(1, from = "nfc", to = "user"), length.out = 100)
    preds_deci <- predict(fit_deci, newdata = data.frame(m = newx), interval = "confidence")
    polygon(c(rev(newx), newx), c(rev(preds_deci[, 3]), preds_deci[, 2]), col = Transparent(deci_col, 0.3), border = NA)

    # evergreen midgdown
    points(midgdown_ever[, .(m, f)],
        cex = 0.5,
        xlab = NA, ylab = NA, pch = 16, col = Transparent(ever_col, 0.3)
    )
    fit_ever <- lm(f ~ m, data = midgdown_ever)
    abline(fit_ever, col = ever_col, lwd = 2, xpd = FALSE)
    preds_ever <- predict(fit_ever, newdata = data.frame(m = newx), interval = "confidence")
    polygon(c(rev(newx), newx), c(rev(preds_ever[, 3]), preds_ever[, 2]), col = Transparent(ever_col, 0.3), border = NA)

    legend("topleft", bty = "n", xpd = TRUE, legend = c("Evergreen", "Deciduous"), pch = 16, col = c(ever_col, deci_col))


    # ~ boxplot subpanel
    # ~~~~~~~~~~~~~~~~~
    par(fig = c(0.35, 0.48, 0.2, 0.5), new = TRUE)
    par(mar = c(0, 0, 0, 0), mgp = c(1.5, 0.5, 0))
    boxplot(c(midgup_ever[, .(f, m)], midgup_deci[, .(f, m)]), xaxt = "n", yaxt = "n", 
        whisklty = c(2, 1), col = rep(c(ever_col, deci_col), each = 2), outline = FALSE, boxwex = 0.5)
    axis(side = 1, at = 1:4, padj = -1, labels = c("F", "M", "F", "M"), cex.axis = 0.8, tck = -0.03, gap.axis = 0.05)
    axis(side = 2, at = seq(50, 250, length.out = 5), labels = seq(50, 250, length.out = 5), cex.axis = 0.8, tck = -0.03, gap.axis = 0.05, las = 2, hadj = 0.7)

    par(fig = c(0.85, 0.98, 0.2, 0.5), new = TRUE)
    par(mar = c(0, 0, 0, 0), mgp = c(1.5, 0.5, 0))
    boxplot(c(midgdown_ever[, .(f, m)], midgdown_deci[, .(f, m)]), xaxt = "n", yaxt = "n", 
        whisklty = c(2, 1), col = rep(c(ever_col, deci_col), each = 2), outline = FALSE, boxwex = 0.5)
    axis(side = 1, at = 1:4, padj = -1, labels = c("F", "M", "F", "M"), cex.axis = 0.8, tck = -0.03, gap.axis = 0.05)
    axis(side = 2, at = seq(150, 350, length.out = 5), labels = seq(150, 350, length.out = 5), cex.axis = 0.8, tck = -0.03, gap.axis = 0.05, las = 2, hadj = 0.7)


    par(fig = c(0, 1, 0, 1), new = TRUE, xpd = NA)
    # dev.off()
    par(opar)
}


#endregion [Pheno compare]
