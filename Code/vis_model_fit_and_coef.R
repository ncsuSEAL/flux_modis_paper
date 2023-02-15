#************************************************************************************
# Description: GPP and EVI2 model fit performance and the standard coefficients.
# Author: Xiaojie(J) Gao
# Date: 2022-03-03
#************************************************************************************
rm(list=ls())

source("Code/base.R")
source("Code/mod_base.R")



focus_biomes <- c("CRO", "DBF", "ENF", "GRA", "MF", "WET")

# Make color consistent
color_pal <- brewer.pal(8, "Dark2")

modis_color <- color_pal[4]
flux_color <- color_pal[3]


# Data comes from "bayesian_model.R"
load(file.path("Pipeline/best_model_fit_VUT.RData"))
load(file.path("Pipeline/model_3_fit_VUT.RData"))


igbp_colors <- c(
    "#51B6F5", "#218A21", "#31CD31", "#9ACD31", "#97FA97", "#8FBB8F", 
    "#BB8F8F", "#F5DEB3", "#DBEB9D", "#FFD600", "#EFB766", "#4682B2", 
    "#FAED73", "#FF0000", "#999355", "#F5F5DC", "#BDBDBD", "#000000"
)
# igbp_names <- c("water", "enf", "ebf", "dnf", "dbf", "mixed", "closed shrubs", 
#    "open shrubs", "woody savannas", "savannas",
# "grasslands", "perm wetlands", "croplands", "urban", "crop/natural mosaic", 
#    "snow and ice", "barren/sparse veg", "unclassified")
igbp_names <- c(
    "water", "ENF", "EBF", "DNF", "DBF", "MF", "CSH", "OSH", "WSA", "SAV",
    "GRA", "WET", "CRO", "urban", "crop/natural mosaic", "snow and ice", 
    "barren/sparse veg", "unclassified"
)

IGBP <- unique(as.character(f_igbp_factor))
# cols <- RColorBrewer::brewer.pal(11, "Paired") # Set these to be static
colsIGBP <- merge(data.frame(IGBP = IGBP), 
    data.frame(IGBP = igbp_names, cols = igbp_colors))
colsIGBP$IGBP <- as.character(colsIGBP$IGBP)
colsIGBP$cols <- as.character(colsIGBP$cols)

f_dt_model_3 <- data.table(fit = f_model_3$fit * 100, agpp = Y * 100, 
    IGBP = f_igbp_factor)
f_dt_model_3 <- merge(f_dt_model_3, colsIGBP)
f_dt_model_3[, ":="(IGBP = as.character(IGBP), cols = as.character(cols))]

m_dt_model_3 <- data.table(fit = m_model_3$fit * 100, agpp = Y * 100, 
    IGBP = f_igbp_factor)
m_dt_model_3 <- merge(m_dt_model_3, colsIGBP)
m_dt_model_3[, ":="(IGBP = as.character(IGBP), cols = as.character(cols))]

f_dt <- data.table(fit = f_model$fit * 100, agpp = Y * 100, IGBP = f_igbp_factor)
f_dt <- merge(f_dt, colsIGBP)
f_dt[, ":="(IGBP = as.character(IGBP), cols = as.character(cols))]

m_dt <- data.table(fit = m_model$fit * 100, agpp = Y * 100, IGBP = f_igbp_factor)
m_dt <- merge(m_dt, colsIGBP)
m_dt[, ":="(IGBP = as.character(IGBP), cols = as.character(cols))]


# Fig: Fig 1
png("Output/model_fit_and_coef.png", width = 2500, height = 1700, res = 300)
# png("0_Manuscript/Fig_1_xx.png", width = 2500, height = 1700, res = 300)

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
mtext(side = 2, expression(EC~~Sigma~GPP), line = 1.5)
text(2000, 900, cex = 1.3, labels = c(
    paste0("R2:", round(cor(f_dt_model_3$fit.fit_median, f_dt_model_3$agpp)^2, 2)),
    paste0("\n\nRMSE:", round(sqrt(mean((f_dt_model_3$fit.fit_median - 
        f_dt_model_3$agpp)^2)), 2))
), adj = 0)
abline(0, 1, lty = 2)
text(grconvertX(0, "npc", "user"), grconvertY(0.93, "npc", "user"), labels = "a", 
    pos = 4, font = 2, cex = 1.5, xpd = NA)

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
mtext(side = 1, expression(Fitted~EC~~Sigma~GPP), line = 2.5)
mtext(side = 2, expression(EC~~Sigma~GPP), line = 1.5)
text(2000, 900, cex = 1.3, labels = c(
    paste0("R2:", round(cor(f_dt$fit.fit_median, f_dt$agpp)^2, 2)),
    paste0("\n\nRMSE:", round(sqrt(mean((f_dt$fit.fit_median - f_dt$agpp)^2)), 2))
), adj = 0)
abline(0, 1, lty = 2)
text(grconvertX(0, "npc", "user"), grconvertY(0.93, "npc", "user"), labels = "b", 
    pos = 4, font = 2, cex = 1.5, xpd = NA)

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
    paste0("\n\nRMSE:", round(sqrt(mean((m_dt_model_3$fit.fit_median - 
        m_dt_model_3$agpp)^2)), 2))
), adj = 0)
abline(0, 1, lty = 2)
text(grconvertX(0, "npc", "user"), grconvertY(0.93, "npc", "user"), labels = "c", 
    pos = 4, font = 2, cex = 1.5, xpd = NA)

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
mtext(side = 1, expression(Fitted~MODIS~~Sigma~GPP), line = 2.5)
# axis(side = 2, at = c(1000, 2000, 3000), labels = c(1000, 2000, 3000))
text(2000, 900, cex = 1.3, labels = c(
    paste0("R2:", round(cor(m_dt$fit.fit_median, m_dt$agpp)^2, 2)),
    paste0("\n\nRMSE:", round(sqrt(mean((m_dt$fit.fit_median - m_dt$agpp)^2)), 2))
), adj = 0)
abline(0, 1, lty = 2)
text(grconvertX(0, "npc", "user"), grconvertY(0.93, "npc", "user"), labels = "d", 
    pos = 4, font = 2, cex = 1.5, xpd = NA)

# ~ 1e
# ~~~~~~~~~~~~~~~~~
par(mar = c(3, 3, 2.5, 1))
# normalized effect
focus_biomes <- c("CRO", "DBF", "ENF", "GRA", "MF", "WET")

plot(NA, xlim = c(-0.4, 1.4), ylim = c(0.5, 6.7), xlab = "", ylab = "", yaxt = "n", 
    mgp = c(1.5, 0.5, 0))
axis(side = 2, at = 1:6, labels = c("EVI2 min", "GPP min", "EVI2 max", "GPP max", 
    "EVI2 GSL", "GPP GSL"), las = 2, cex.axis = 1.3)
mtext(side = 1, "Normalized Effect", line = 2)
abline(v = 0, lty = 2, col = "grey50")

EtaPointRange <- function(eta_df, colsIGBP_df, focus_biomes, ycoord) {
    tmp_df <- merge(eta_df, colsIGBP_df)
    tmp_df$cols <- as.character(tmp_df$cols)
    points(tmp_df[IGBP %in% focus_biomes, med / 100], ycoord, pch = 21, lwd = 0.3, 
        bg = tmp_df[IGBP %in% focus_biomes, cols], cex = 1.3)
    segments(
        x0 = tmp_df[IGBP %in% focus_biomes, lwr / 100], y0 = ycoord,
        x1 = tmp_df[IGBP %in% focus_biomes, upr / 100], y1 = ycoord, 
        col = tmp_df[IGBP %in% focus_biomes, cols], lwd = 1.3
    )
}

EtaPointRange(f_model_norm$eta2, colsIGBP, focus_biomes, 
    seq(6.1, 6.4, length = 6) - 0.2)
EtaPointRange(f_model_norm$eta3, colsIGBP, focus_biomes, 
    seq(4.1, 4.4, length = 6) - 0.2)
EtaPointRange(f_model_norm$eta4, colsIGBP, focus_biomes, 
    seq(2.1, 2.4, length = 6) - 0.2)
ycoord <- c(2, 4, 6) - 0.25
segments(x0 = f_model_norm$mu[4:2, lwr / 100], y0 = ycoord, 
    x1 = f_model_norm$mu[4:2, upr / 100], y1 = ycoord, col = 1, lwd = 2)
points(f_model_norm$mu[4:2, med / 100], ycoord, pch = 21, bg = flux_color, 
    col = 1, cex = 2, lwd = 1.5)

EtaPointRange(m_model_norm$eta2, colsIGBP, focus_biomes, 
    seq(5.1, 5.4, length = 6) - 0.2)
EtaPointRange(m_model_norm$eta3, colsIGBP, focus_biomes, 
    seq(3.1, 3.4, length = 6) - 0.2)
EtaPointRange(m_model_norm$eta4, colsIGBP, focus_biomes, 
    seq(1.1, 1.4, length = 6) - 0.2)
ycoord <- c(1, 3, 5) + 0.3
segments(x0 = m_model_norm$mu[4:2, lwr / 100], y0 = ycoord, 
    x1 = m_model_norm$mu[4:2, upr / 100], y1 = ycoord, col = 1, lwd = 2)
points(m_model_norm$mu[4:2, med / 100], ycoord, pch = 21, bg = modis_color, 
    col = 1, cex = 2, lwd = 1.5)
legend("bottomright", legend = c("MODIS EVI2", "Flux GPP"), pch = 21, 
    pt.bg = c(modis_color, flux_color), pt.lwd = 1.5, col = 1, bty = "n", cex = 1.3)
text(grconvertX(0, "npc", "user"), grconvertY(0.966, "npc", "user"), 
    labels = "e", pos = 4, font = 2, cex = 1.5, xpd = NA)

# ~ Annoatations
# ~~~~~~~~~~~~~~~~~
# text(grconvertX(0, "ndc", "user"), grconvertY(0.95, "ndc", "user"), 
#    labels = "A", pos = 4, font = 2, cex = 1.5, xpd = NA)
# text(grconvertX(0, "ndc", "user"), grconvertY(0.48, "ndc", "user"), 
#    labels = "B", pos = 4, font = 2, cex = 1.5, xpd = NA)
# text(grconvertX(0.35, "ndc", "user"), grconvertY(0.95, "ndc", "user"), 
#    labels = "C", pos = 4, font = 2, cex = 1.5, xpd = NA)
# text(grconvertX(0.35, "ndc", "user"), grconvertY(0.48, "ndc", "user"), 
#    labels = "D", pos = 4, font = 2, cex = 1.5, xpd = NA)
# text(grconvertX(0.68, "ndc", "user"), grconvertY(0.95, "ndc", "user"), 
#    labels = "E", pos = 4, font = 2, cex = 1.5, xpd = NA)

# ~ Legend
# ~~~~~~~~~~~~~~~~~
legend(grconvertX(0.5, "ndc", "user"), grconvertY(0.07, "ndc", "user"),
    pch = 21, pt.lwd = 0.3, pt.bg = as.character(colsIGBP$cols),
    legend = colsIGBP$IGBP,
    xpd = NA, ncol = 11, bty = "n", cex = 1.3, xjust = 0.5
)

dev.off()
