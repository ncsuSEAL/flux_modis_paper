#************************************************************************
# Description: Plot the EVI2 and GPP time series for all site years
# Author: Xiaojie Gao
# Date: 2021-03-20
#************************************************************************
source("Code/base.R")


start_date <- as.Date("2000-01-01")
end_date <- as.Date("2014-12-31")


plotEviGPPForSite <- function(cur_site) {
    # cur_site <- "AR-SLu"
    par(mfrow = c(2, 1), mgp = c(1.5, 0.5, 0))

    # ~ EVI2 time series
    # ~~~~~~~~~~~~~~~~~
    par(mar = c(0, 3, 2, 1))
    plot(NA, xlim = c(start_date, end_date), ylim = c(0, 1), xaxt = "n", xlab = "", ylab = "EVI2")
    rect(
        xleft = seq(start_date, end_date, by = "year")[seq(1, 15, by = 2)], ybottom = rep(grconvertY(0.005, "npc", "user"), 8),
        xright = seq(start_date, end_date + 1, by = "year")[seq(2, 16, by = 2)], ytop = rep(grconvertY(0.995, "npc", "user"), 8),
        col = "grey90", border = NA
    )
    title(paste0(cur_site, " (", site_meta[siteID == cur_site, IGBP], ")"), line = 0.5, cex.main = 0.8)

    points(modis_evi2[site == cur_site, .(date, evi2)], pch = 16, cex = 0.5, col = Transparent("grey50", 0.3))
    # first cycle
    apply(modis_pheno[site_name == cur_site, .(
        m_gup_1, m_midgup_1, m_maturity_1, m_peak_1, m_sene_1, m_midgdown_1, m_dormancy_1
    )], 1, function(x) {
        abline(v = x, col = rev(viridis(8)[1:7]))
    })
    # second cycle
    apply(modis_pheno[site_name == cur_site, .(
        m_gup_2, m_midgup_2, m_maturity_2, m_peak_2, m_sene_2, m_midgdown_2, m_dormancy_2
    )], 1, function(x) {
        abline(v = x, col = rev(viridis(8)[1:7]), lty = 2)
    })
    legend("top",
        xpd = TRUE, xjust = 0.5,
        legend = c("1st cycle", "2nd cycle", "Greenup", "MidGreenup", "Maturity", "Peak", "Senescence", "MidGreendown", "Dormancy"),
        lty = c(1, 2, rep(1, 7)), col = c("black", "black", rev(viridis(8)[1:7])), ncol = 9, bg = "white", cex = 0.8
    )


    # ~ GPP time series
    # ~~~~~~~~~~~~~~~~~
    par(mar = c(3, 3, 0, 1))
    plot(NA, xlim = c(start_date, end_date), ylim = c(0, 32), xaxt = "n", xlab = "", ylab = "GPP")
    rect(
        xleft = seq(start_date, end_date, by = "year")[seq(1, 15, by = 2)], ybottom = rep(grconvertY(0.005, "npc", "user"), 8),
        xright = seq(start_date, end_date + 1, by = "year")[seq(2, 16, by = 2)], ytop = rep(grconvertY(0.995, "npc", "user"), 8),
        col = "grey90", border = NA
    )
    points(flux_gpp_ts_2[site == cur_site, .(date, gpp_raw)], pch = 16, cex = 0.5, col = Transparent("grey50", 0.3))
    lines(flux_gpp_ts_2[site == cur_site, .(date, abs(gpp_pred))], pch = 16, cex = 0.5, col = Transparent("black", 0.5), lwd = 2)
    # first cyclflux_gpp_ts_2e
    apply(flux_pheno[site == cur_site, .(
        f_gup_1, f_midgup_1, f_maturity_1, f_peak_1, f_sene_1, f_midgdown_1, f_dormancy_1
    )], 1, function(x) {
        abline(v = x, col = rev(viridis(8)[1:7]))
    })
    # second cycle
    apply(flux_pheno[site == cur_site, .(
        f_gup_2, f_midgup_2, f_maturity_2, f_peak_2, f_sene_2, f_midgdown_2, f_dormancy_2
    )], 1, function(x) {
        abline(v = x, col = rev(viridis(8)[1:7]), lty = 2)
    })

    # x-axis
    axis(side = 1, at = seq(start_date, end_date + 1, by = "year"), labels = year(start_date):(year(end_date) + 1))
}

modis_pheno <- LoadModisPheno()
modis_evi2 <- LoadModisEvi2TS()
flux_pheno <- LoadFluxPheno()
flux_gpp_ts_2 <- LoadFluxGppTs()

site_meta <- LoadFluxSiteMeta()

pdf("evi2_with_gpp_for_all_sites.pdf", width = 16, height = 6)
res <- sapply(unique(modis_pheno$site), function(x) { plotEviGPPForSite(x)})
dev.off()

# unique(modis_pheno$site)[!(unique(modis_pheno$site) %in% unique(modis_evi2$site))]
# modis_pheno[site_name == "AU-DaS"]
# plotEviGPPForSite("US-Ha1")
