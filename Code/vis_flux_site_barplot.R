#************************************************************************************
# Description: Flux site and site years barplots
# Author: Xiaojie(J) Gao
# Date: 2022-03-03
#************************************************************************************
rm(list=ls())

source("Code/base.R")
library(data.table)


# Extracted flux pheno data
flux <- fread(file.path(gdir, "flux_q2algo_pheno.csv"))

# GPP missing value
gpp_miss <- fread(file.path(gdir, "gpp_lookup.csv"))

flux <- merge(flux, gpp_miss,
    by.x = c("site_name", "year"),
    by.y = c("siteID", "year"), all = TRUE)
# Get site years that have GPP measurements, no matter if phenometrics extracted
flux <- flux[miss == FALSE,]


# Site metadata
site_meta <- fread(file.path(gdir, "sites_fluxnet2015.csv"))
site_meta <- site_meta[data_tier2015 == 1, ]

# Merge flux data with site metadata
flux <- merge(flux[, .(site_name, year)], site_meta[, .(siteID, IGBP)],
    by.x = "site_name", by.y = "siteID")

# Number of sites per IGBP
lc_sites <- unique(flux[, .(site_name, IGBP)])[, .N, by = "IGBP"]
# sum(lc_sites$N)

# Number of site years per IGBP
lc_site_years <- flux[, .N, by = c("IGBP")]
# sum(lc_site_years$N)



# Color palette
colorsTable <- brewer.pal(5, "Dark2")


# fig: Number of flux site and site years analyzed
png(filename = "Output/flux_site_barplot.png", bg = NA,
    width = 1000, height = 600, res = 150)

par(mgp = c(1.5, 0.5, 0), mar = c(3, 4, 1, 2))
layout(matrix(c(1:2), nrow = 1))
lc_sites$color <- colorsTable[2]
lc_sites[IGBP %in% c("ENF", "EBF")]$color <- colorsTable[1]

bar <- barplot(lc_sites$N,
    las = 2, names.arg = lc_sites$IGBP, horiz = TRUE,
    ylab = "", xlab = "Number of sites", xlim = c(0, 40),
    col = lc_sites$color, las = 1
)
mtext(side = 2, "IGBP type", line = 2.5)
# title("flux sites per LC")
text(bar, x = lc_sites$N + 2.5, labels = lc_sites$N, xpd = TRUE)
# legend("topright", fill = c(colorsTable[1], colorsTable[2]), legend = c("Evergreen", "Deciduous"), bty = "n")
text(grconvertX(0.03, "nfc", "user"), grconvertY(0.92, "nfc", "user"), "a", xpd = TRUE)


# setorder(lc_site_years, cols = -N)
lc_site_years$color <- colorsTable[2]
lc_site_years[IGBP %in% c("ENF", "EBF")]$color <- colorsTable[1]
bar <- barplot(lc_site_years$N,
    las = 2, names.arg = lc_site_years$IGBP, horiz = TRUE,
    ylab = "", xlab = "Number of site years", xlim = c(0, 300),
    col = lc_site_years$color, las = 1
)
# title("flux site-years per LC")
mtext(side = 2, "IGBP type", line = 2.5)
text(bar, x = lc_site_years$N + 27, labels = lc_site_years$N, xpd = TRUE)
text(grconvertX(0.03, "nfc", "user"), grconvertY(0.92, "nfc", "user"), "b", xpd = TRUE)
legend("topright",
    fill = c(colorsTable[1], colorsTable[2]),
    legend = c("Evergreen", "Deciduous"), bty = "n"
)

dev.off()