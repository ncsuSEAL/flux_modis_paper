rm(list=ls())
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

ever_col <- color_pal[1]
deci_col <- color_pal[2]

if (dir.exists("SI") == FALSE) {
    dir.create("SI")
}


#region [Fig 1 baseplot version] ####

#! Data comes from "bayesian_model.R"
load(file.path("Pipeline/best_model_fit_VUT.RData"))
load(file.path("Pipeline/model_3_fit_VUT.RData"))


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
png("Output/Fig_1_xx.png", width = 2500, height = 1700, res = 300)

layout(matrix(c(1, 2, 3, 4, 5, 5), nrow = 2))
par(mgp = c(2, 0.5, 0), cex.lab = 1.3, cex.axis = 1.2, oma = c(3, 0, 0, 0))

# ~ 1a
# ~~~~~~~~~~~~~~~~~
par(mar = c(0, 4, 2.5, 0))
# GPP model 3 fit
plot(f_dt_model_3[, .(fit.fit_median, agpp)],
    bg = f_dt_model_3$cols, pch = 21, lwd = 0.3,
    xlim = range(f_dt$fit.fit_median, f_dt$agpp),
    ylim = range(f_dt$fit.fit_median, f_dt$agpp),
    xlab = "", ylab = "",
    xaxt = "n", yaxt = "n"
)
# axis(side = 1, at = c(1000, 2000, 3000), labels = c(1000, 2000, 3000))
axis(side = 2, at = c(1000, 2000, 3000), labels = c(1000, 2000, 3000))
mtext(side = 2, expression("True ΣGPP" ~ (g ~ C ~ m^"-2" ~ yr^"-1")), line = 1.5)
text(2000, 900, cex = 1.3, labels = c(
    paste0("R2:", round(cor(f_dt_model_3$fit.fit_median, f_dt_model_3$agpp)^2, 2)),
    paste0("\n\nRMSE:", round(sqrt(mean((f_dt_model_3$fit.fit_median - f_dt_model_3$agpp)^2)), 2))
), adj = 0)
abline(0, 1, lty = 2)
text(grconvertX(0, "npc", "user"), grconvertY(0.93, "npc", "user"), labels = "a", pos = 4, font = 2, cex = 1.5, xpd = NA)

# ~ 1b
# ~~~~~~~~~~~~~~~~~
par(mar = c(3, 4, 0, 0))
# GPP best model fit
plot(f_dt[, .(fit.fit_median, agpp)],
    bg = f_dt$cols, pch = 21, lwd = 0.3,
    xlim = range(f_dt$fit.fit_median, f_dt$agpp),
    ylim = range(f_dt$fit.fit_median, f_dt$agpp),
    xlab = "", ylab = "",
    xaxt = "n", yaxt = "n"
)
axis(side = 1, at = c(1000, 2000, 3000), labels = c(1000, 2000, 3000))
axis(side = 2, at = c(1000, 2000, 3000), labels = c(1000, 2000, 3000))
mtext(side = 1, expression("Fitted ΣGPP" ~ (g ~ C ~ m^"-2"~ yr^"-1")), line = 2.5)
mtext(side = 2, expression("True ΣGPP" ~ (g ~ C ~ m^"-2"~ yr^"-1")), line = 1.5)
text(2000, 900, cex = 1.3, labels = c(
    paste0("R2:", round(cor(f_dt$fit.fit_median, f_dt$agpp)^2, 2)),
    paste0("\n\nRMSE:", round(sqrt(mean((f_dt$fit.fit_median - f_dt$agpp)^2)), 2))
), adj = 0)
abline(0, 1, lty = 2)
text(grconvertX(0, "npc", "user"), grconvertY(0.93, "npc", "user"), labels = "b", pos = 4, font = 2, cex = 1.5, xpd = NA)

# ~ 1c
# ~~~~~~~~~~~~~~~~~
par(mar = c(0, 0, 2.5, 3.5))
# EVI2 model 3 fit
plot(m_dt_model_3[, .(fit.fit_median, agpp)],
    bg = m_dt_model_3$cols, pch = 21, lwd = 0.3,
    xlim = range(m_dt_model_3$fit.fit_median, m_dt_model_3$agpp),
    ylim = range(m_dt_model_3$fit.fit_median, m_dt_model_3$agpp),
    xlab = "", ylab = "",
    xaxt = "n", yaxt = "n"
)
# axis(side = 1, at = c(1000, 2000, 3000), labels = c(1000, 2000, 3000))
# axis(side = 2, at = c(1000, 2000, 3000), labels = c(1000, 2000, 3000))
text(2000, 900, cex = 1.3, labels = c(
    paste0("R2:", round(cor(m_dt_model_3$fit.fit_median, m_dt_model_3$agpp)^2, 2)),
    paste0("\n\nRMSE:", round(sqrt(mean((m_dt_model_3$fit.fit_median - m_dt_model_3$agpp)^2)), 2))
), adj = 0)
abline(0, 1, lty = 2)
text(grconvertX(0, "npc", "user"), grconvertY(0.93, "npc", "user"), labels = "c", pos = 4, font = 2, cex = 1.5, xpd = NA)

# ~ 1d
# ~~~~~~~~~~~~~~~~~
par(mar = c(3, 0, 0, 3.5))
# EVI2 best model fit
plot(m_dt[, .(fit.fit_median, agpp)],
    bg = m_dt$cols, pch = 21, lwd = 0.3,
    xlim = range(m_dt$fit.fit_median, m_dt$agpp),
    ylim = range(m_dt$fit.fit_median, m_dt$agpp),
    xlab = "", ylab = "",
    xaxt = "n", yaxt = "n"
)
axis(side = 1, at = c(1000, 2000, 3000), labels = c(1000, 2000, 3000))
mtext(side = 1, expression("Fitted ΣGPP" ~ (g ~ C ~ m^"-2"~ yr^"-1")), line = 2.5)
# axis(side = 2, at = c(1000, 2000, 3000), labels = c(1000, 2000, 3000))
text(2000, 900, cex = 1.3, labels = c(
    paste0("R2:", round(cor(m_dt$fit.fit_median, m_dt$agpp)^2, 2)),
    paste0("\n\nRMSE:", round(sqrt(mean((m_dt$fit.fit_median - m_dt$agpp)^2)), 2))
), adj = 0)
abline(0, 1, lty = 2)
text(grconvertX(0, "npc", "user"), grconvertY(0.93, "npc", "user"), labels = "d", pos = 4, font = 2, cex = 1.5, xpd = NA)

# ~ 1e
# ~~~~~~~~~~~~~~~~~
par(mar = c(3, 3, 2.5, 1))
# normalized effect
focus_biomes <- c("CRO", "DBF", "ENF", "GRA", "MF", "WET")

plot(NA, xlim = c(-0.4, 1.4), ylim = c(0.5, 6.7), xlab = "", ylab = "", yaxt = "n", mgp = c(1.5, 0.5, 0))
axis(side = 2, at = 1:6, labels = c("EVI2 min", "GPP min", "EVI2 max", "GPP max", "EVI2 GSL", "GPP GSL"), las = 2, cex.axis = 1.3)
mtext(side = 1, "Normalized Effect", line = 2)
abline(v = 0, lty = 2, col = "grey50")

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
segments(x0 = f_model_norm$mu[4:2, lwr / 100], y0 = ycoord, x1 = f_model_norm$mu[4:2, upr / 100], y1 = ycoord, col = 1, lwd = 2)
points(f_model_norm$mu[4:2, med / 100], ycoord, pch = 21, bg = flux_color, col = 1, cex = 2, lwd = 1.5)

EtaPointRange(m_model_norm$eta2, colsIGBP, focus_biomes, seq(5.1, 5.4, length = 6) - 0.2)
EtaPointRange(m_model_norm$eta3, colsIGBP, focus_biomes, seq(3.1, 3.4, length = 6) - 0.2)
EtaPointRange(m_model_norm$eta4, colsIGBP, focus_biomes, seq(1.1, 1.4, length = 6) - 0.2)
ycoord <- c(1, 3, 5) + 0.3
segments(x0 = m_model_norm$mu[4:2, lwr / 100], y0 = ycoord, x1 = m_model_norm$mu[4:2, upr / 100], y1 = ycoord, col = 1, lwd = 2)
points(m_model_norm$mu[4:2, med / 100], ycoord, pch = 21, bg = modis_color, col = 1, cex = 2, lwd = 1.5)
legend("bottomright", legend = c("MODIS EVI2", "Flux GPP"), pch = 21, pt.bg = c(modis_color, flux_color), pt.lwd = 1.5, col = 1, bty = "n", cex = 1.3)
text(grconvertX(0, "npc", "user"), grconvertY(0.966, "npc", "user"), labels = "e", pos = 4, font = 2, cex = 1.5, xpd = NA)

# ~ Annoatations
# ~~~~~~~~~~~~~~~~~
# text(grconvertX(0, "ndc", "user"), grconvertY(0.95, "ndc", "user"), labels = "A", pos = 4, font = 2, cex = 1.5, xpd = NA)
# text(grconvertX(0, "ndc", "user"), grconvertY(0.48, "ndc", "user"), labels = "B", pos = 4, font = 2, cex = 1.5, xpd = NA)
# text(grconvertX(0.35, "ndc", "user"), grconvertY(0.95, "ndc", "user"), labels = "C", pos = 4, font = 2, cex = 1.5, xpd = NA)
# text(grconvertX(0.35, "ndc", "user"), grconvertY(0.48, "ndc", "user"), labels = "D", pos = 4, font = 2, cex = 1.5, xpd = NA)
# text(grconvertX(0.68, "ndc", "user"), grconvertY(0.95, "ndc", "user"), labels = "E", pos = 4, font = 2, cex = 1.5, xpd = NA)

# ~ Legend
# ~~~~~~~~~~~~~~~~~
legend(grconvertX(0.5, "ndc", "user"), grconvertY(0.07, "ndc", "user"), 
    pch = 21, pt.lwd = 0.3, pt.bg = as.character(colsIGBP$cols), 
    legend = colsIGBP$IGBP, 
    xpd = NA, ncol = 11, bty = "n", cex = 1.3, xjust = 0.5)

dev.off()
}
#endregion [Fig 1 baseplot version]



#region [Fig 2] ####

#! Data comes from "bayesian_model.R"
load(file.path("Pipeline/best_model_fit_VUT.RData"))

gsl_coef_fig <- ggplot() +
    # add modis
    geom_pointrange(aes(x = IGBP, ymin = lwr, y = med, ymax = upr, color = modis_color), shape = 16, alpha = 0.5, 
        size = 0.3, data = m_model$beta2[IGBP %in% focus_biomes], position = position_dodge2(width = 0.5)) +
    geom_pointrange(aes(x = IGBP, ymin = lwr, y = med, ymax = upr, fill = modis_color), shape = 21, stroke = 1, 
        size = 0.8, data = m_model$eta2[IGBP %in% focus_biomes],  position = position_nudge(x = -0.35)) +
    geom_pointrange(aes(x = IGBP, ymin = lwr, y = med, ymax = upr, color = flux_color), shape = 16, alpha = 0.5, 
        size = 0.3, data = f_model$beta2[IGBP %in% focus_biomes], position = position_dodge2(width = 0.5)) +
    geom_pointrange(aes(x = IGBP, ymin = lwr, y = med, ymax = upr, fill = flux_color), shape = 21, size = 0.8, 
        data = f_model$eta2[IGBP %in% focus_biomes],  position = position_nudge(x = -0.45)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
    xlab("Biome type") +
    ylab(expression("GSL Coefficient" ~ (g ~ C ~ m^"-2" ~ yr^"-1"))) +
    geom_text(aes(x = IGBP, y = rep(-20, 6), label = N), 
        data = north_sites_dt[, .N, by = "IGBP"][IGBP %in% focus_biomes]) +
    geom_text(aes(x = IGBP, y = rep(-15, 6), label = N), 
        data = unique(north_sites_dt[, .(site, IGBP)])[IGBP %in% focus_biomes, .N, by = c("IGBP")]) +
    theme_article() +
    scale_y_continuous(breaks = seq(-40, 40, by = 10)) +
    scale_color_identity() +
    scale_fill_manual(name = "", values = c("#E7298A" = "#E7298A", "#7570B3" = "#7570B3"), labels = c("MODIS EVI2", "FLUX GPP"), guide = "legend") +
    theme(legend.position = c(0.9, 0.93)) +
    guides(fill = guide_legend(
        keywidth = 0.1,
        keyheight = 0.3,
        default.unit = "inch"
    ))

gsl_coef_fig

ggsave(filename = "Output/Fig_2_xx.png", plot = gsl_coef_fig, width = 7, height = 4, bg = "white")

#endregion [Fig 2]



#region [Fig 3] ####

# ! Data comes from "flux_modis_summary.Rmd".

# The comparion result
comp_result_all <- fread(file.path("Data/comparison_new.csv"))
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




{
png("Output/Fig_3.png", width = 2100, height = 1200, res = 300)
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
boxplot(c(midgup_ever[, .(f, m)], midgup_deci[, .(f, m)]), xaxt = "n", yaxt = "n", col = rep(c(ever_col, deci_col), each = 2), outline = FALSE, boxwex = 0.5)
axis(side = 1, at = 1:4, padj = -1, labels = c("F", "M", "F", "M"), cex.axis = 0.8, tck = -0.03, gap.axis = 0.05)
axis(side = 2, at = seq(50, 250, length.out = 5), labels = seq(50, 250, length.out = 5), cex.axis = 0.8, tck = -0.03, gap.axis = 0.05, las = 2, hadj = 0.7)

par(fig = c(0.85, 0.98, 0.2, 0.5), new = TRUE)
par(mar = c(0, 0, 0, 0), mgp = c(1.5, 0.5, 0))
boxplot(c(midgdown_ever[, .(f, m)], midgdown_deci[, .(f, m)]), xaxt = "n", yaxt = "n", col = rep(c(ever_col, deci_col), each = 2), outline = FALSE, boxwex = 0.5)
axis(side = 1, at = 1:4, padj = -1, labels = c("F", "M", "F", "M"), cex.axis = 0.8, tck = -0.03, gap.axis = 0.05)
axis(side = 2, at = seq(150, 350, length.out = 5), labels = seq(150, 350, length.out = 5), cex.axis = 0.8, tck = -0.03, gap.axis = 0.05, las = 2, hadj = 0.7)


par(fig = c(0, 1, 0, 1), new = TRUE, xpd = NA)
text(grconvertX(0, "nfc", "user"), grconvertY(0.95, "nfc", "user"), labels = "a", pos = 4, font = 2)
text(grconvertX(0.5, "nfc", "user"), grconvertY(0.95, "nfc", "user"), labels = "b", pos = 4, font = 2)

dev.off()
}
#endregion [Fig 3]



#region [Fig 4] ####

plotEviGPPForSite <- function(cur_site, start_date, end_date) {
    start_date <- as.Date(start_date)
    end_date <- as.Date(end_date)

    par(mfrow = c(2, 1), mgp = c(1.5, 0.5, 0))

    # ~ EVI2 time series
    # ~~~~~~~~~~~~~~~~~
    par(mar = c(0, 3, 3, 1))

    evi_ts <- modis_evi2[site == cur_site, .(date, evi2)]
    plot(NA, xlim = c(start_date, end_date), ylim = c(0, max(evi_ts$evi2)), xaxt = "n", xlab = "", ylab = "EVI2", cex.axis = 0.8)
    rect(
        xleft = seq(start_date, end_date, by = "year")[seq(1, 15, by = 2)], ybottom = rep(grconvertY(0.005, "npc", "user"), 8),
        xright = seq(start_date, end_date + 1, by = "year")[seq(2, 16, by = 2)], ytop = rep(grconvertY(0.995, "npc", "user"), 8),
        col = "grey90", border = NA
    )
    # title(paste0(cur_site, " (", site_meta[siteID == cur_site, IGBP], ")"), line = 0.5, cex.main = 0.8)

    points(evi_ts, pch = 16, cex = 0.5, col = Transparent("grey50", 0.3))
    # first cycle
    apply(modis_pheno[site_name == cur_site & between(year, year(start_date), year(end_date)), 
        .(m_midgup_1, m_midgdown_1)], 1, function(x) {
        abline(v = x, col = rev(viridis(8)[c(2, 6)]))
    })
    # second cycle
    apply(modis_pheno[site_name == cur_site & between(year, year(start_date), year(end_date)), .(
        m_midgup_2, m_midgdown_2
    )], 1, function(x) {
        abline(v = x, col = rev(viridis(8)[c(2, 6)]), lty = 2)
    })
    legend(grconvertX(0.5, "nfc", "user"), grconvertY(0.85, "nfc", "user"),
        xpd = TRUE, xjust = 0.5, bty = "n",
        legend = c("1st cycle", "2nd cycle", "MidGreenup", "MidGreendown"),
        lty = c(1, 2, rep(1, 2)), col = c("black", "black", rev(viridis(8)[c(2, 6)])), ncol = 4, bg = "white", cex = 0.8
    )


    # ~ GPP time series
    # ~~~~~~~~~~~~~~~~~
    par(mar = c(3, 3, 0, 1))
    
    gpp_ts <- flux_gpp_ts[site == cur_site, .(date, gpp)]
    # gpp_ts <- flux_gpp_ts[site == cur_site, .(date, gpp_raw)]
    plot(NA, xlim = c(start_date, end_date), ylim = c(0, max(gpp_ts$gpp)), xaxt = "n", xlab = "", ylab = "GPP", cex.axis = 0.8)
    rect(
        xleft = seq(start_date, end_date, by = "year")[seq(1, 15, by = 2)], ybottom = rep(grconvertY(0.005, "npc", "user"), 8),
        xright = seq(start_date, end_date + 1, by = "year")[seq(2, 16, by = 2)], ytop = rep(grconvertY(0.995, "npc", "user"), 8),
        col = "grey90", border = NA
    )
    points(flux_gpp_ts[site == cur_site, .(date, gpp)], pch = 16, cex = 0.5, col = Transparent("grey50", 0.3))
    # points(flux_gpp_ts[site == cur_site, .(date, gpp_raw)], pch = 16, cex = 0.5, col = Transparent("grey50", 0.3))
    # lines(flux_gpp_ts[site == cur_site, .(date, abs(gpp_pred))], pch = 16, cex = 0.5, col = Transparent("black", 0.5), lwd = 2)
    # first cyclflux_gpp_ts_2e
    apply(flux_pheno[site == cur_site & between(year, year(start_date), year(end_date)), .(
        f_midgup_1, f_midgdown_1
    )], 1, function(x) {
        abline(v = x, col = rev(viridis(8)[c(2, 6)]))
    })
    # second cycle
    apply(flux_pheno[site == cur_site & between(year, year(start_date), year(end_date)), .(
        f_midgup_2, f_midgdown_2
    )], 1, function(x) {
        abline(v = x, col = rev(viridis(8)[1:7]), lty = 2)
    })

    # x-axis
    axis(side = 1, at = seq(start_date, end_date + 1, by = "year"), labels = year(start_date):(year(end_date) + 1))
}

modis_evi2 <- LoadModisEvi2TS()
modis_pheno <- LoadModisPheno()
flux_gpp_ts <- LoadFluxGppTs(v2 = FALSE)
flux_pheno <- LoadFluxPheno()

pdf("Output/Fig_4_components_xx.pdf", width = 8, height = 4)
plotEviGPPForSite(cur_site = "US-Ha1", start_date = "2004-01-01", end_date = "2010-12-31")
plotEviGPPForSite(cur_site = "IT-BCi", start_date = "2009-01-01", end_date = "2012-12-31")
plotEviGPPForSite(cur_site = "CZ-BK1", start_date = "2005-01-01", end_date = "2008-12-31")
plotEviGPPForSite(cur_site = "CH-Fru", start_date = "2008-01-01", end_date = "2014-12-31")
# plotEviGPPForSite(cur_site = "US-Los", start_date = "2001-01-01", end_date = "2007-12-31")

dev.off()


#endregion [Fig 4]







#region [SI Fig GPPmax vs. EVI2max] ####

# ~ Log transformation plot on the original scale
# ~~~~~~~~~~~~~~~~~
# compute site means across years
north_sites_dt[, ":=" (m_EVImax_1_mean = mean(m_EVImax_1, na.rm = TRUE), 
    f_gppmax_1_mean = mean(f_gppmax_1, na.rm = TRUE)), by = "site"]

{
    png("SI/SI_Fig_6.png", width = 800, height = 800, res = 170)

    colorsTable <- brewer.pal(5, "Dark2")
    plot(north_sites_dt[cat == "ever", .(m_EVImax_1_mean, f_gppmax_1_mean)],
        mgp = c(1.5, 0.5, 0), xlab = "EVI2max", ylab = "GPPmax", pch = 16,
        col = ever_col, xlim = c(0.1, 0.9), ylim = c(2, 30)
    )
    fit_ever <- lm(log(f_gppmax_1_mean) ~ m_EVImax_1_mean, data = north_sites_dt[cat == "ever"])
    # abline(fit_ever, col = ever_col, lwd = 2)
    newx <- seq(grconvertX(0, from = "nfc", to = "user"), grconvertX(1, from = "nfc", to = "user"), length.out = 100)
    preds_ever <- predict(fit_ever, newdata = data.frame(m_EVImax_1_mean = newx), interval = "confidence")
    lines(newx, exp(preds_ever[, 1]), lwd = 2, col = ever_col)
    polygon(c(rev(newx), newx), c(rev(exp(preds_ever[, 3])), exp(preds_ever[, 2])), col = Transparent(ever_col, 0.3), border = NA)
    text(0.1, 20, col = ever_col, labels = c(
        paste0("y = ", "exp(", round(coef(fit_ever)[2], 2), "x +", round(coef(fit_ever)[1], 2), ")"),
        paste0("\n\nR2:", round(summary(fit_ever)$r.squared, 2))
    ), adj = 0)

    points(north_sites_dt[cat == "deci", .(m_EVImax_1_mean, f_gppmax_1_mean)], col = deci_col, pch = 16)
    fit_deci <- lm(log(f_gppmax_1_mean) ~ m_EVImax_1_mean, data = north_sites_dt[cat == "deci"])
    # abline(fit_deci, col = deci_col, lwd = 2)
    newx <- seq(grconvertX(0, from = "nfc", to = "user"), grconvertX(1, from = "nfc", to = "user"), length.out = 100)
    preds_deci <- predict(fit_deci, newdata = data.frame(m_EVImax_1_mean = newx), interval = "confidence")
    lines(newx, exp(preds_deci[, 1]), lwd = 2, col = deci_col)
    polygon(c(rev(newx), newx), c(rev(exp(preds_deci[, 3])), exp(preds_deci[, 2])), col = Transparent(deci_col, 0.3), border = NA)
    text(0.1, 15, col = deci_col, labels = c(
        paste0("y = ", "exp(", round(coef(fit_deci)[2], 2), "x +", round(coef(fit_deci)[1], 2), ")"),
        paste0("\n\nR2:", round(summary(fit_deci)$r.squared, 2))
    ), adj = 0)

    legend("topleft", bty = "n", col = colorsTable[1:2], pch = 16, legend = c("Evergreen", "Deciduous"))
    
    dev.off()
}

#endregion [SI Fig GPPmax vs. EVI2max]



#region [SI fig. variability for all biome types] ####

load(file.path("Pipeline/best_model_fit_VUT.RData"))

gsl_coef_fig <- ggplot() +
    # add modis
    geom_pointrange(aes(x = IGBP, ymin = lwr, y = med, ymax = upr, color = modis_color), shape = 16, alpha = 0.3, 
        size = 0.3, data = m_model$beta2, position = position_dodge2(width = 0.5)) +
    geom_pointrange(aes(x = IGBP, ymin = lwr, y = med, ymax = upr, fill = modis_color), shape = 21, stroke = 1, 
        size = 0.8, data = m_model$eta2,  position = position_nudge(x = -0.35)) +
    geom_pointrange(aes(x = IGBP, ymin = lwr, y = med, ymax = upr, color = flux_color), shape = 16, alpha = 0.3, 
        size = 0.3, data = f_model$beta2, position = position_dodge2(width = 0.5)) +
    geom_pointrange(aes(x = IGBP, ymin = lwr, y = med, ymax = upr, fill = flux_color), shape = 21, size = 0.8, 
        data = f_model$eta2,  position = position_nudge(x = -0.45)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
    ylab(expression("GSL Coefficient" ~ (g ~ C ~ m^"-2"))) +
    geom_text(aes(x = IGBP, y = rep(-20, 11), label = N), 
        data = north_sites_dt[, .N, by = "IGBP"]) +
    geom_text(aes(x = IGBP, y = rep(-15, 11), label = N), 
        data = unique(north_sites_dt[, .(site, IGBP)])[, .N, by = c("IGBP")]) +
    theme_article() +
    scale_y_continuous(breaks = seq(-40, 40, by = 10)) +
    scale_color_identity() +
    scale_fill_manual(name = "", values = c("#E7298A" = "#E7298A", "#7570B3" = "#7570B3"), labels = c("MODIS EVI2", "FLUX GPP"), guide = "legend") +
    theme(legend.position = c(0.9, 0.93)) +
    guides(fill = guide_legend(
        keywidth = 0.1,
        keyheight = 0.3,
        default.unit = "inch"
    ))

gsl_coef_fig

ggsave(filename = "SI/SI_Fig_7.png", plot = gsl_coef_fig, width = 10, height = 4, bg = "white")

#endregion [SI fig. variability for all biome types]



#region [SI fig. comparsion on the original scale] ####

load(file.path("Pipeline/best_model_fit_VUT.RData"))

ggplot() +
    geom_boxplot(aes(x = "GPPmax", y = med), data = f_model$eta2, width = 0.3) +
    geom_boxplot(aes(x = "EVI2max", y = med), data = m_model$eta2, width = 0.3) +
    xlab("EVI2/GPP max effect") +
    ylab("Biome-level effect variability") +
    theme_article()

#endregion [SI fig. comparsion on the original scale]



#region [SI fig. comparison between models with GSL effect only] ####

load(file.path("Pipeline", "gsl_only_model.RData"))
load(file.path("Pipeline", "model_3_fit.RData"))

f_model_gsl_only$mu$med[2] - mean(unlist(f_model_3$mu_beta[2])) * 100
m_model_gsl_only$mu$med[2] - mean(unlist(m_model_3$mu_beta[2])) * 100

f_mu_beta2 <- data.frame(t(quantile(unlist(f_model_3$mu_beta[2]) * 100, c(0.025, 0.5, 0.975))))
colnames(f_mu_beta2) <- c("lwr", "med", "upr")

m_mu_beta2 <- data.frame(t(quantile(unlist(m_model_3$mu_beta[2]) * 100, c(0.025, 0.5, 0.975))))
colnames(m_mu_beta2) <- c("lwr", "med", "upr")

ggplot() +
    geom_pointrange(aes(x = "GPP GSL only", ymin = lwr, y = med, ymax = upr), data = f_model_gsl_only$mu[2]) +
    geom_pointrange(aes(x = "GPP with max and min", ymin = lwr, y = med, ymax = upr), data = f_mu_beta2) +
    geom_pointrange(aes(x = "EVI2 GSL only", ymin = lwr, y = med, ymax = upr), data = m_model_gsl_only$mu[2]) +
    geom_pointrange(aes(x = "EVI2 with max and min", ymin = lwr, y = med, ymax = upr), data = m_mu_beta2) +
    xlab("Models") +
    ylab("GSL effect") +
    theme_article()


quantile(m_model_3$beta2[,5], c(0.025, 0.5, 0.975))

#endregion [SI fig. comparison between models with GSL effect only]



#region [GSL among biomes] ####

gsl_plot <- ggplot(north_sites_dt) +
    geom_boxplot(aes(x = IGBP, y = m_gsl, fill = modis_color), position = position_nudge(x = -0.2), width = 0.3) +
    geom_boxplot(aes(x = IGBP, y = f_gsl, fill = flux_color), position = position_nudge(x = 0.2), width = 0.3) +
    scale_fill_manual(name = "", values = c("#E7298A" = "#E7298A", "#7570B3" = "#7570B3"), labels = c("MODIS EVI2", "FLUX GPP"), guide = "legend") +
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
    annotate("text", x = 0.5, y = 100, label = "Difference = FLUX - MODIS", hjust = 0, vjust = 0) +
    theme_article() +
    scale_fill_identity()

compose <- gsl_plot / gsl_diff_plot + plot_annotation(tag_levels = "a")

ggsave(file = "SI/SI_Fig_9.png", plot = compose, width = 8, height = 6, bg = "white")

#endregion [GSL among biomes]



#region [Some reports] ####

load(file.path("Pipeline/best_model_fit_VUT.RData"))

f_model$beta2[IGBP == "CRO"]
m_model$beta2[IGBP == "CRO"]

t.test(f_model$beta2[IGBP == "CRO", med], m_model$beta2[IGBP == "CRO", med])
t.test(f_model$beta2[IGBP == "DBF", med], m_model$beta2[IGBP == "DBF", med])
t.test(f_model$beta2[IGBP == "MF", med], m_model$beta2[IGBP == "MF", med])
t.test(f_model$beta2[IGBP == "ENF", med], m_model$beta2[IGBP == "ENF", med])
t.test(f_model$beta2[IGBP == "GRA", med], m_model$beta2[IGBP == "GRA", med])
t.test(f_model$beta2[IGBP == "WET", med], m_model$beta2[IGBP == "WET", med])



# ~ For all biome types
# ~~~~~~~~~~~~~~~~~
# overall GSL-GPP effect estimated by flux and LSP data
f_model$mu[2]
m_model$mu[2]

# how much did LSP underestimate GSL-GPP effect
1 - mean(m_model$eta2$med) / mean(f_model$eta2$med)

# compare biome-level difference
var(f_model$eta2$med) # biome-level variability from flux data
var(m_model$eta2$med) # biome-level variability from LSP data

var(m_model$eta2$med) / var(f_model$eta2$med)

# compare site-level difference
var(f_model$beta2$med) # site-level variability from flux data
var(m_model$beta2$med) # site-level variability from LSP data

var(m_model$beta2$med) / var(f_model$beta2$med)


# ~ For focused biome types
# ~~~~~~~~~~~~~~~~~
# overall GSL-GPP effect estimated by flux and LSP data
mean(f_model$eta2[IGBP %in% focus_biomes, med])
mean(m_model$eta2[IGBP %in% focus_biomes, med])


# how much did LSP underestimate GSL-GPP effect
1 - mean(m_model$eta2[IGBP %in% focus_biomes, med]) / mean(f_model$eta2[IGBP %in% focus_biomes, med])

# compare biome-level difference
var(f_model$eta2[IGBP %in% focus_biomes, med]) # biome-level variability from flux data
var(m_model$eta2[IGBP %in% focus_biomes, med]) # biome-level variability from LSP data

var(m_model$eta2[IGBP %in% focus_biomes, med]) / var(f_model$eta2[IGBP %in% focus_biomes, med])

# compare site-level difference
var(f_model$beta2[IGBP %in% focus_biomes, med]) # site-level variability from flux data
var(m_model$beta2[IGBP %in% focus_biomes, med]) # site-level variability from LSP data

var(m_model$beta2[IGBP %in% focus_biomes, med]) / var(f_model$beta2[IGBP %in% focus_biomes, med])

#endregion [Some reports]



#region [SI Fig cross validation result] ####

load(file.path("Pipeline", "cv.RData"))

uni_yrs <- sort(unique(north_sites_dt$year))

{
    png("SI/SI_Fig_5.png", width = 1000, height = 1000, res = 150)
    par(mar = c(3, 4, 2, 2))
    plot(uni_yrs, sqrt(f_cv_mse_1), type = "b", ylim = c(0, 500), col = flux_color, pch = 16, 
        xlab = "Year", ylab = expression("RMSE" ~ (g ~ C ~ m^"-2")), mgp = c(1.5, 0.5, 0), lwd = 1)
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
        ncol = 3, legend = c("Model 1", "Model 2", "Model 3", "Model 4", "Model 5", "GPP based", "EVI2 based"),
        pch = c(16, 22, 2, 4, 7, NA, NA), lwd = 1, col = c(rep("black", 5), flux_color, modis_color), bty = "n"
    )
    dev.off()
}

#endregion [SI Fig S4 cross validation result]



#region [Sankey diagram] ####

library(networkD3)

load(file.path("Pipeline", "match_processed.RData"))

pheno_cols_1 <- paste("f", pheno_names, 1, sep = "_")
pheno_cols_2 <- paste("f", pheno_names, 2, sep = "_")

match_dt[, matchstate_1 := "unkown"]
match_dt[, matchstate_2 := "unkown"]

for (i in 1:nrow(match_dt)) {
    non_na_sum_1 <- 0
    for (j in pheno_cols_1) {
        non_na_sum_1 <- non_na_sum_1 + ifelse(!is.na(match_dt[i, get(j)]), 1, 0)
    }
    if (non_na_sum_1 == 7) {
        match_dt[i, matchstate_1 := "full"]
    } else if (non_na_sum_1 < 7 & non_na_sum_1 > 0) {
        match_dt[i, matchstate_1 := "part"]
    } else if (non_na_sum_1 == 0) {
        match_dt[i, matchstate_1 := "none"]
    }

    non_na_sum_2 <- 0
    for (j in pheno_cols_2) {
        non_na_sum_2 <- non_na_sum_2 + ifelse(!is.na(match_dt[i, get(j)]), 1, 0)
    }
    if (non_na_sum_2 == 7) {
        match_dt[i, matchstate_2 := "full"]
    } else if (non_na_sum_2 < 7 & non_na_sum_2 > 0) {
        match_dt[i, matchstate_2 := "part"]
    } else if (non_na_sum_2 == 0) {
        match_dt[i, matchstate_2 := "none"]
    }
}

full_match <- match_dt[matchstate_1 == "full" | matchstate_2 == "full", ]
part_match <- match_dt[matchstate_1 == "part" | matchstate_2 == "part", ]

links <- data.frame(source=NA, target=NA, value=NA)

# first/second link
links <- rbind(links, c("Fluxnet2015 (GPP)", "No GPP metrics", nrow(flux[is.na(f_gup_1),])))
links <- rbind(links, c("Fluxnet2015 (GPP)", "GPP metrics", nrow(flux[!is.na(f_gup_1),])))
links <- rbind(links, c("MCD12Q2 (EVI2)", "No EVI2 phenometrics", nrow(gpp_evi_low) + nrow(no_match_modis_na)))
links <- rbind(links, c("MCD12Q2 (EVI2)", "EVI2 phenometrics", nrow(flux) - nrow(gpp_evi_low) - nrow(no_match_modis_na)))

links <- rbind(links, c("No GPP metrics", "MODIS non-NA, but FLux NA", nrow(flux_gpp_but_NA) - nrow(gpp_evi_low)))
links <- rbind(links, c("No GPP metrics", "GPP and EVI2 amplitude both too low", nrow(gpp_evi_low)))
links <- rbind(links, c("GPP metrics", "Full cycle match", nrow(full_match)))
links <- rbind(links, c("GPP metrics", "Part cycle match", nrow(part_match)))
links <- rbind(links, c("GPP metrics", "Both have data but don't match", nrow(no_match_real)))
links <- rbind(links, c("GPP metrics", "Flux non-NA, but MODIS NA", nrow(no_match_modis_na)))

links <- rbind(links, c("No EVI2 phenometrics", "GPP and EVI2 amplitude both too low", nrow(gpp_evi_low)))
links <- rbind(links, c("No EVI2 phenometrics", "Flux non-NA, but MODIS NA", nrow(no_match_modis_na)))
links <- rbind(links, c("EVI2 phenometrics", "Full cycle match", nrow(full_match)))
links <- rbind(links, c("EVI2 phenometrics", "Part cycle match", nrow(part_match)))
links <- rbind(links, c("EVI2 phenometrics", "Both have data but don't match", nrow(no_match_real)))
links <- rbind(links, c("EVI2 phenometrics", "MODIS non-NA, but FLux NA", nrow(no_match_flux_na)))

links <- links[-1,]

nodes <- data.frame(
  name = unique(c(as.character(links$source), as.character(links$target)))
)

nodes$group <- as.factor(c("flux", "modis", "gpp_low", "gpp_metrics", "evi_low", 
    "evi_phenometrics", "flux_na", "both_low", "match", "match", "no_match", "modis_na"))

links$IDsource <- match(links$source, nodes$name) - 1
links$IDtarget <- match(links$target, nodes$name) - 1

links$group <- as.factor(c("gpp_low", "flux", "evi_low", "modis", "gpp_low", 
    "both_low", "match", "match", "no_match","modis_na", "both_low", "evi_low", "match", "match", 
    "no_match", "flux_na"))


my_color <- "d3.scaleOrdinal()
  .domain(['flux', 'modis', 'gpp_low', 'evi_low', 'both_low', 'match', 'no_match', 'modis_na', 'flux_na'])
  .range(['#7570B3', '#E7298A', 'grey', 'grey', '#1B9E77', '#1B9E77', 'red', 'grey', 'grey'])"

sankeyNetwork(Links = links, Nodes = nodes,
              Source = "IDsource", Target = "IDtarget",
              Value = "value", NodeID = "name",
              sinksRight=FALSE, fontSize = 14, colourScale = my_color,
              LinkGroup = "group", NodeGroup = "group")

#endregion [Sankey diagram]




#region [MODIS GSL vs. FLUX GSL] ####

# compute site means across years
north_sites_dt[, ":="(m_gsl_mean = mean(m_gsl, na.rm = TRUE),
    f_gsl_mean = mean(f_gsl, na.rm = TRUE)), by = "site"]

{
    png("SI/SI_Fig_gsl_evi_gpp.png", width = 1000, height = 1000, res = 170)

    colorsTable <- brewer.pal(5, "Dark2")
    plot(north_sites_dt[cat == "ever", .(m_gsl_mean, f_gsl_mean)],
        mgp = c(1.5, 0.5, 0), xlab = "EVI2 GSL", ylab = "GPP GSL", pch = 16,
        col = ever_col, xlim = c(50, 250), ylim = c(50, 250)
    )
    abline(0, 1, lty = 2)
    fit_ever <- lm(f_gsl_mean ~ m_gsl_mean, data = north_sites_dt[cat == "ever"])
    # abline(fit_ever, col = ever_col, lwd = 2)
    newx <- seq(grconvertX(0, from = "nfc", to = "user"), grconvertX(1, from = "nfc", to = "user"), length.out = 100)
    preds_ever <- predict(fit_ever, newdata = data.frame(m_gsl_mean = newx), interval = "confidence")
    lines(newx, preds_ever[, 1], lwd = 2, col = ever_col)
    polygon(c(rev(newx), newx), c(rev(preds_ever[, 3]), preds_ever[, 2]), col = Transparent(ever_col, 0.3), border = NA)
    text(50, 220, col = ever_col, labels = c(
        paste0("y = ", round(coef(fit_ever)[2], 2), "x +", round(coef(fit_ever)[1], 2)),
        paste0("\n\nR2:", round(summary(fit_ever)$r.squared, 2))
    ), adj = 0)

    points(north_sites_dt[cat == "deci", .(m_gsl_mean, f_gsl_mean)], col = deci_col, pch = 16)
    fit_deci <- lm(f_gsl_mean ~ m_gsl_mean, data = north_sites_dt[cat == "deci"])
    # abline(fit_deci, col = deci_col, lwd = 2)
    newx <- seq(grconvertX(0, from = "nfc", to = "user"), grconvertX(1, from = "nfc", to = "user"), length.out = 100)
    preds_deci <- predict(fit_deci, newdata = data.frame(m_gsl_mean = newx), interval = "confidence")
    lines(newx, preds_deci[, 1], lwd = 2, col = deci_col)
    polygon(c(rev(newx), newx), c(rev(preds_deci[, 3]), preds_deci[, 2]), col = Transparent(deci_col, 0.3), border = NA)
    text(50, 200, col = deci_col, labels = c(
        paste0("y = ", round(coef(fit_deci)[2], 2), "x +", round(coef(fit_deci)[1], 2)),
        paste0("\n\nR2:", round(summary(fit_deci)$r.squared, 2))
    ), adj = 0)

    legend("topleft", bty = "n", col = colorsTable[1:2], pch = 16, legend = c("Evergreen", "Deciduous"))

    dev.off()
}


#endregion [MODIS GSL vs. FLUX GSL]


#region [Annual GPP vs. MODIS GSL and FLUX GSL] ####
north_sites_dt <- LoadModelARD(north_only = TRUE)
dt <- unique(
    north_sites_dt[, .(
        m_gsl_mean = mean(m_gsl, na.rm = TRUE),
        f_gsl_mean = mean(f_gsl, na.rm = TRUE),
        m_EVImax_1_mean = mean(m_EVImax_1, na.rm = TRUE),
        f_gppmax_1_mean = mean(f_gppmax_1, na.rm = TRUE),
        agpp_mean = mean(annual_gpp, na.rm = TRUE),
        cat,
        IGBP
    ),
    by = "site"
    ]
)


{
    png("SI/SI_fig_xx.png", width = 1000, height = 1000, res = 120)
    par(mfrow = c(2, 2), mgp = c(1.5, 0.5, 0), mar = c(3, 3, 2, 2))

    # ~ EVI2 GSL vs annual GPP
    # ~~~~~~~~~~~~~~~~~
    xrange <- range(dt$m_gsl_mean)
    yrange <- range(dt$agpp_mean)
    # plot evergreen points
    plot(dt[cat == "ever", .(m_gsl_mean, agpp_mean)], pch = 16, col = ever_col,
        xlim = xrange, ylim = yrange, xlab = "EVI2 GSL", ylab = "Annual GPP")
    # fit lm model
    fit_ever <- lm(agpp_mean ~ m_gsl_mean, data = dt[cat == "ever"])
    abline(fit_ever, col = ever_col, lwd = 2)
    newx <- seq(grconvertX(0, from = "nfc", to = "user"), grconvertX(1, from = "nfc", to = "user"), length.out = 100)
    preds_ever <- predict(fit_ever, newdata = data.frame(m_gsl_mean = newx), interval = "confidence")
    polygon(c(rev(newx), newx), c(rev(preds_ever[, 3]), preds_ever[, 2]), col = Transparent(ever_col, 0.3), border = NA)
    text(150, 800, col = ever_col, labels = c(
        paste0("y = ", round(coef(fit_ever)[2], 2), "x +", round(coef(fit_ever)[1], 2)),
        paste0("\n\nR2:", round(summary(fit_ever)$r.squared, 2))
    ), adj = 0)

    # plot deciduous points
    points(dt[cat == "deci", .(m_gsl_mean, agpp_mean)], pch = 16, col = deci_col)
    # fit lm model
    fit_deci <- lm(agpp_mean ~ m_gsl_mean, data = dt[cat == "deci"])
    abline(fit_deci, col = deci_col, lwd = 2)
    newx <- seq(grconvertX(0, from = "nfc", to = "user"), grconvertX(1, from = "nfc", to = "user"), length.out = 100)
    preds_deci <- predict(fit_deci, newdata = data.frame(m_gsl_mean = newx), interval = "confidence")
    polygon(c(rev(newx), newx), c(rev(preds_deci[, 3]), preds_deci[, 2]), col = Transparent(deci_col, 0.3), border = NA)
    text(150, 500, col = deci_col, labels = c(
        paste0("y = ", round(coef(fit_deci)[2], 2), "x +", round(coef(fit_deci)[1], 2)),
        paste0("\n\nR2:", round(summary(fit_deci)$r.squared, 2))
    ), adj = 0)

    legend("topleft", bty = "n", col = colorsTable[1:2], pch = 16, legend = c("Evergreen", "Deciduous"))


    # ~ GPP GSL vs annual GPP
    # ~~~~~~~~~~~~~~~~~
    xrange <- range(dt$f_gsl_mean)
    yrange <- range(dt$agpp_mean)
    # plot evergreen points
    plot(dt[cat == "ever", .(f_gsl_mean, agpp_mean)],
        pch = 16, col = ever_col,
        xlim = c(50, 250), ylim = c(350, 2700), xlab = "GPP GSL", ylab = "Annual GPP"
    )
    # fit lm model
    fit_ever <- lm(agpp_mean ~ f_gsl_mean, data = dt[cat == "ever"])
    abline(fit_ever, col = ever_col, lwd = 2)
    newx <- seq(grconvertX(0, from = "nfc", to = "user"), grconvertX(1, from = "nfc", to = "user"), length.out = 100)
    preds_ever <- predict(fit_ever, newdata = data.frame(f_gsl_mean = newx), interval = "confidence")
    polygon(c(rev(newx), newx), c(rev(preds_ever[, 3]), preds_ever[, 2]), col = Transparent(ever_col, 0.3), border = NA)
    text(150, 800, col = ever_col, labels = c(
        paste0("y = ", round(coef(fit_ever)[2], 2), "x +", round(coef(fit_ever)[1], 2)),
        paste0("\n\nR2:", round(summary(fit_ever)$r.squared, 2))
    ), adj = 0)

    # plot deciduous points
    points(dt[cat == "deci", .(f_gsl_mean, agpp_mean)], pch = 16, col = deci_col)
    # fit lm model
    fit_deci <- lm(agpp_mean ~ f_gsl_mean, data = dt[cat == "deci"])
    abline(fit_deci, col = deci_col, lwd = 2)
    newx <- seq(grconvertX(0, from = "nfc", to = "user"), grconvertX(1, from = "nfc", to = "user"), length.out = 100)
    preds_deci <- predict(fit_deci, newdata = data.frame(f_gsl_mean = newx), interval = "confidence")
    polygon(c(rev(newx), newx), c(rev(preds_deci[, 3]), preds_deci[, 2]), col = Transparent(deci_col, 0.3), border = NA)
    text(150, 500, col = deci_col, labels = c(
        paste0("y = ", round(coef(fit_deci)[2], 2), "x +", round(coef(fit_deci)[1], 2)),
        paste0("\n\nR2:", round(summary(fit_deci)$r.squared, 2))
    ), adj = 0)
    

    # ~ EVI2 max vs annual GPP
    # ~~~~~~~~~~~~~~~~~
    xrange <- range(dt$m_EVImax_1_mean)
    yrange <- range(dt$agpp_mean)
    # plot evergreen points
    plot(dt[cat == "ever", .(m_EVImax_1_mean, agpp_mean)],
        pch = 16, col = ever_col,
        xlim = c(0.3, 0.7), ylim = c(350, 2700), xlab = "EVI2max", ylab = "Annual GPP"
    )
    # fit lm model
    fit_ever <- lm(agpp_mean ~ m_EVImax_1_mean, data = dt[cat == "ever"])
    abline(fit_ever, col = ever_col, lwd = 2)
    newx <- seq(grconvertX(0, from = "nfc", to = "user"), grconvertX(1, from = "nfc", to = "user"), length.out = 100)
    preds_ever <- predict(fit_ever, newdata = data.frame(m_EVImax_1_mean = newx), interval = "confidence")
    polygon(c(rev(newx), newx), c(rev(preds_ever[, 3]), preds_ever[, 2]), col = Transparent(ever_col, 0.3), border = NA)
    text(0.3, 2600, col = ever_col, labels = c(
        paste0("y = ", round(coef(fit_ever)[2], 2), "x +", round(coef(fit_ever)[1], 2)),
        paste0("\n\nR2:", round(summary(fit_ever)$r.squared, 2))
    ), adj = 0)

    # plot deciduous points
    points(dt[cat == "deci", .(m_EVImax_1_mean, agpp_mean)], pch = 16, col = deci_col)
    # fit lm model
    fit_deci <- lm(agpp_mean ~ m_EVImax_1_mean, data = dt[cat == "deci"])
    abline(fit_deci, col = deci_col, lwd = 2)
    newx <- seq(grconvertX(0, from = "nfc", to = "user"), grconvertX(1, from = "nfc", to = "user"), length.out = 100)
    preds_deci <- predict(fit_deci, newdata = data.frame(m_EVImax_1_mean = newx), interval = "confidence")
    polygon(c(rev(newx), newx), c(rev(preds_deci[, 3]), preds_deci[, 2]), col = Transparent(deci_col, 0.3), border = NA)
    text(0.3, 2300, col = deci_col, labels = c(
        paste0("y = ", round(coef(fit_deci)[2], 2), "x +", round(coef(fit_deci)[1], 2)),
        paste0("\n\nR2:", round(summary(fit_deci)$r.squared, 2))
    ), adj = 0)


    # ~ GPP max vs annual GPP
    # ~~~~~~~~~~~~~~~~~
    xrange <- range(dt$f_gppmax_1_mean)
    yrange <- range(dt$agpp_mean)
    # plot evergreen points
    plot(dt[cat == "ever", .(f_gppmax_1_mean, agpp_mean)],
        pch = 16, col = ever_col,
        xlim = xrange, ylim = yrange, xlab = "GPPmax", ylab = "Annual GPP"
    )
    # fit lm model
    fit_ever <- lm(agpp_mean ~ f_gppmax_1_mean, data = dt[cat == "ever"])
    abline(fit_ever, col = ever_col, lwd = 2)
    newx <- seq(grconvertX(0, from = "nfc", to = "user"), grconvertX(1, from = "nfc", to = "user"), length.out = 100)
    preds_ever <- predict(fit_ever, newdata = data.frame(f_gppmax_1_mean = newx), interval = "confidence")
    polygon(c(rev(newx), newx), c(rev(preds_ever[, 3]), preds_ever[, 2]), col = Transparent(ever_col, 0.3), border = NA)
    text(11, 1000, col = ever_col, labels = c(
        paste0("y = ", round(coef(fit_ever)[2], 2), "x +", round(coef(fit_ever)[1], 2)),
        paste0("\n\nR2:", round(summary(fit_ever)$r.squared, 2))
    ), adj = 0)

    # plot deciduous points
    points(dt[cat == "deci", .(f_gppmax_1_mean, agpp_mean)], pch = 16, col = deci_col)
    # fit lm model
    fit_deci <- lm(agpp_mean ~ f_gppmax_1_mean, data = dt[cat == "deci"])
    abline(fit_deci, col = deci_col, lwd = 2)
    newx <- seq(grconvertX(0, from = "nfc", to = "user"), grconvertX(1, from = "nfc", to = "user"), length.out = 100)
    preds_deci <- predict(fit_deci, newdata = data.frame(f_gppmax_1_mean = newx), interval = "confidence")
    polygon(c(rev(newx), newx), c(rev(preds_deci[, 3]), preds_deci[, 2]), col = Transparent(deci_col, 0.3), border = NA)
    text(11, 600, col = deci_col, labels = c(
        paste0("y = ", round(coef(fit_deci)[2], 2), "x +", round(coef(fit_deci)[1], 2)),
        paste0("\n\nR2:", round(summary(fit_deci)$r.squared, 2))
    ), adj = 0)
    
    dev.off()
}

#endregion [Annual GPP vs. MODIS GSL and FLUX GSL]
