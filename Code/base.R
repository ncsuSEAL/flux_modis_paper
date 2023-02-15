#*******************************************************************************
# Description: Basic functions and variables needed for other scripts
# 
# Author: Xiaojie Gao
# Date: 2021-03-02
#*******************************************************************************
library(data.table)
library(viridis)
library(RColorBrewer)
library(rgdal)
library(raster)
library(ggplot2)
library(patchwork)

# Since this repo is used to reproduce the paper, data should be stored in the
#  repo folder as well. So, the `gdir` and `hpc_dir` should be the same.

gdir <- "Data"
hpc_dir <- gdir

if (dir.exists("Pipeline") == FALSE) {
    dir.create("Pipeline")
}

# ~ Load data
# ~~~~~~~~~~~~~~~~~
# NOTE: site metadata
LoadFluxSiteMeta <- function() {
    site_meta <- fread(file.path(
        gdir,
        "Data/Flux/metadata_fluxnet2015/sites_fluxnet2015.csv"
    ))
    site_meta <- site_meta[data_tier2015 == 1, ] # only use Tier 1
    return(site_meta)
}

# NOTE: MODIS EVI2 time series
LoadModisEvi2TS <- function() {
    modis_evi2 <- readRDS(file.path(gdir, "modis_evi2.Rds"))
    return(modis_evi2)
}

# NOTE: MODIS EVI2-based phenometrics
LoadModisPheno <- function() {
    modis_pheno <- readRDS(file.path(gdir, "modis_pheno.Rds"))
    return(modis_pheno)
}

# NOTE: FLUX GPP-based phenometrics
LoadFluxPheno <- function() {
    flux_pheno <- readRDS(file.path(gdir, "flux_pheno.Rds"))
    return(flux_pheno)
}

# NOTE: FLUX GPP time series from the `gpp_timeseries` folder
LoadFluxGppTs <- function(v2 = TRUE) {
    # NOTE: FLUX GPP VUT time series
    if (v2 == TRUE) {
        flux_gpp_ts <- readRDS(file.path(gdir, "flux_gpp_ts_2.Rds"))
    } else {
        flux_gpp_ts <- readRDS(file.path(gdir, "flux_gpp.Rds"))
    }
    return(flux_gpp_ts)
}

LoadMcd12q2 <- function() {
    modis <- fread(file.path(
        gdir, "Data/fluxsites-MCD12Q2-006-singlepixel.csv"
    ))
    return(modis)
}

LoadMatchDt <- function(north_only = FALSE) {
    match_dt <- if (north_only == FALSE) {
        # NOTE: FLUX-MODIS match result
        readRDS(file.path(gdir, "Pipeline/match_result.Rds"))
    } else {
        # NOTE: FLUX-MODIS match result of north hemisphere sites)
        readRDS(file.path(gdir, "Pipeline/north_sites_dt.Rds"))
    }

    return(match_dt)
}

LoadModelARD <- function(north_only = FALSE) {
    if (north_only == FALSE) {
        # A giant table containing MODIS phenos and Flux phenos & GPP
        dt <- readRDS(file.path(gdir, "giant_dt.Rds"))
    } else {
        # The giant table for north sites only.
        #    It also contains UMD, LAI, and other classification for each site.
        dt <- readRDS(file.path(gdir, "north_st_2.Rds"))
        # we only focus on the first cycle
        dt <- dt[
            !is.na(f_gppmin_1) & !is.na(f_gpparea_1) &
                !is.na(f_gppamp_1) & !is.na(m_EVImin_1) &
                !is.na(m_EVIarea_1) & !is.na(m_EVIamp_1),
        ]
        dt[, f_gsl := f_midgdown_1 - f_midgup_1]
        dt[, f_gppmax_1 := f_gppamp_1 + f_gppmin_1]
        dt[, m_gsl := m_midgdown_1 - m_midgup_1]
        dt[, m_EVImax_1 := m_EVIamp_1 + m_EVImin_1]
        
        # annual GPP provided by Ian
        # provided_annualgpp <- fread("north_sites_dtCompare.csv") 
        # dt$annual_gpp <- provided_annualgpp$GPP_DT_VUT_REF
        # dt <- dt[!is.na(annual_gpp),]
    }

    return(dt)
}



fill_val <- -9999 # fill value
pheno_names <- c("gup", "midgup", "maturity", "peak", 
    "sene", "midgdown", "dormancy"
)

#' Make a standard color transparent.
#' This function is borrowed from 'yarrr' package, but I changed the trans.val 
#' to use alpha value directly.
#' @param orig.col: the original color, can be a color name, a hexadecimal code, 
#' or a rgb vector.
#' @param alpha: define the transparent level.
#' @param maxColorValue: used to convert the color to rgb format before making 
#' it transparent.
#' @example: color <- Transparent("red", 0.5)
Transparent <- function(orig.col, alpha = 1, maxColorValue = 255) {
    n.cols <- length(orig.col)
    orig.col <- col2rgb(orig.col)
    final.col <- rep(NA, n.cols)
    for (i in 1:n.cols) {
        final.col[i] <- rgb(orig.col[1, i], orig.col[2, i], orig.col[3, i],
        alpha = alpha * 255,
        maxColorValue = maxColorValue)
    }
    return(final.col)
}


# Plot EVI2 and GPP time series along with phenometrics for the 
# interested site and year
PlotEviGpp <- function(cur_site, cur_year) {
    yoi <- c(
        as.Date(paste0(as.integer(cur_year) - 1, "-06-30")), 
        as.Date(paste0(as.integer(cur_year) + 1, "-06-30"))
    )
        
    # flux pheno
    f_pheno <- flux_pheno[
        site_name == cur_site & year %in% c(
            as.integer(cur_year) - 1, 
            cur_year, 
            as.integer(cur_year) + 1
        ), 
        .(
            f_gup_1, f_midgup_1, f_maturity_1, f_peak_1, f_sene_1, 
            f_midgdown_1, f_dormancy_1,
            f_gup_2, f_midgup_2, f_maturity_2, f_peak_2, f_sene_2, 
            f_midgdown_2, f_dormancy_2
    )]
    # modis pheno
    m_pheno <- modis_pheno[
        site_name == cur_site & year %in% c(
            as.integer(cur_year) - 1, 
            cur_year, 
            as.integer(cur_year) + 1
        ), .(
            m_gup_1, m_midgup_1, m_maturity_1, m_peak_1, m_sene_1, 
            m_midgdown_1, m_dormancy_1,
            m_gup_2, m_midgup_2, m_maturity_2, m_peak_2, m_sene_2, 
            m_midgdown_2, m_dormancy_2
    )]
    m_pheno_qa <- modis_pheno[
        site_name == cur_site & year %in% c(
            as.integer(cur_year) - 1, 
            cur_year, 
            as.integer(cur_year) + 1
        ), .(
            gup_1 = QA_Detailed_0_Greenup_Description, 
            midgup_1 = QA_Detailed_0_MidGreenup_Description, 
            maturity_1 = QA_Detailed_0_Maturity_Description,
            peak_1 = QA_Detailed_0_Peak_Description, 
            sene_1 = QA_Detailed_0_Senescence_Description, 
            midgdown_1 = QA_Detailed_0_MidGreendown_Description, 
            dormancy_1 = QA_Detailed_0_Dormancy_Description, 
            gup_2 = QA_Detailed_1_Greenup_Description, 
            midgup_2 = QA_Detailed_1_MidGreenup_Description, 
            maturity_2 = QA_Detailed_1_Maturity_Description,
            peak_2 = QA_Detailed_1_Peak_Description, 
            sene_2 = QA_Detailed_1_Senescence_Description, 
            midgdown_2 = QA_Detailed_1_MidGreendown_Description, 
            dormancy_2 = QA_Detailed_1_Dormancy_Description
        )
    ]
    
    
    layout(matrix(1:2, nrow = 2))
    par(mar = c(0, 4, 4, 4))

    # ~ plot EVI2 time series
    # ~~~~~~~~~~~~~~~~~
    evi2_ts <- modis_evi2[site == cur_site & between(date, yoi[1], yoi[2]), ]
    evi2_ts[evi2 == 0, evi2 := NA]
    
    # if there's no EVI2 observations
    if (nrow(evi2_ts) == 0) {
        plot(NA, xlim = yoi, ylim = c(0, 1), 
            xaxt = "n", xlab = "", ylab = "", 
            cex.axis = 0.6
        )
        title(c(cur_site, cur_year))
        mtext("EVI2", side = 2, line = 2)
        # modis phenometric
        legend("topright", legend = "EVI2 missing", bg = "white", cex = 0.8)
    
        par(mar = c(4, 4, 0, 4))
        plot(NA, xlim = xlim, ylim = c(0, 1), xaxt = "n", xlab = "", ylab = "")
        # PlotYearBackground(cur$year)
        mtext("GPP", side = 2, line = 2)
        mtext("Date", side = 1, line = 2)
    }

    plot(evi2_ts[, .(date, evi2)], 
        ylim = c(0, 1), col = rgb(0.3, 0.3, 0.3, 0.3), pch = 16, xaxt = "n", 
        xlab = "", ylab = "", cex.axis = 0.6
    )
    text(grconvertX(0, "npc", "user"), grconvertY(1.3, "npc", "user"),
         labels = paste(
            cur_site, cur_year, "(", site_meta[siteID == cur_site, IGBP], ",", 
            ifelse(site_meta[siteID == cur_site, lat] < 0, "South", "North"), ")"
        ), 
        xpd = TRUE, adj = c(1, 0), cex = 1.2, pos = 4
    )
    # PlotYearBackground(cur$year)
    mtext("EVI2", side = 2, line = 2)
    # modis phenometric
    if (nrow(m_pheno) != 0) {
        abline(v = m_pheno[1, ], col = rep(rev(viridis(8)[1:7]), 2), 
            lty = rep(c(1, 2), each = 7)
        )
        abline(v = m_pheno[2, ], col = rep(rev(viridis(8)[1:7]), 2), 
            lty = rep(c(1, 2), each = 7)
        )
        abline(v = m_pheno[3, ], col = rep(rev(viridis(8)[1:7]), 2), 
            lty = rep(c(1, 2), each = 7)
        )

        text(m_pheno[1, ], 0, labels = m_pheno_qa[1, ], cex = 0.5)
        text(m_pheno[2, ], 0, labels = m_pheno_qa[2, ], cex = 0.5)
        text(m_pheno[3, ], 0, labels = m_pheno_qa[3, ], cex = 0.5)
    }
     
    legend(grconvertX(0.5, "npc", "user"), grconvertY(1.3, "npc", "user"),
        xpd = TRUE, xjust = 0.5, bty = "n",
        legend = c("1st cycle", "2nd cycle", "Greenup", "MidGreenup", 
            "Maturity", "Peak", "Senescence", "MidGreendown", "Dormancy"
        ),
        lty = c(1, 2, rep(1, 7)), 
        col = c("black", "black", rev(viridis(8)[1:7])), 
        ncol = 5, bg = "white", cex = 1
    )


    # plot GPP time series
    # ~~~~~~~~~~~~~~~~~
    gpp_ts <- flux_gpp_ts[site == cur_site & between(date, yoi[1], yoi[2])]
    
    par(mar = c(4, 4, 0, 4))
    plot(gpp_ts[, .(date, gpp)], pch = 16, col = rgb(0.3, 0.3, 0.3, 0.3), 
        xlab = "", ylab = "", cex.axis = 0.6
    )
    # PlotYearBackground(cur$year)
    mtext("GPP", side = 2, line = 2)
    mtext("Date", side = 1, line = 2)
    points(gpp_ts[qc != 0, .(date, gpp)], pch = 16, 
        col = rgb(0.8, 0.1, 0, 0.5)
    )
    # flux phenometrics
    if(nrow(f_pheno) != 0) {
        abline(v = f_pheno[1, ], col = rep(rev(viridis(8)[1:7]), 2), 
            lty = rep(c(1, 2), each = 7)
        )
        abline(v = f_pheno[2, ], col = rep(rev(viridis(8)[1:7]), 2), 
            lty = rep(c(1, 2), each = 7)
        )
        abline(v = f_pheno[3, ], col = rep(rev(viridis(8)[1:7]), 2), 
            lty = rep(c(1, 2), each = 7)
        )
    }
}


PlotStretchLegend <- function(r, colorRamp, digits = 0, labSize = 1, 
    col.axis = "black", ...
) {
    pal <- colorRampPalette(colorRamp)
    qs <- quantile(r, c(0, 0.02, 0.98, 1))
    r_breaks <- c(qs[1], seq(qs[2], qs[3], len = 255), qs[4])
    plot(r, col = pal(length(r_breaks) - 1), breaks = r_breaks, axes = F, 
        box = F, legend = F, ...
    )
    # add a reasonable legend
    legend_at <- round(
        seq(r_breaks[2], r_breaks[length(r_breaks) - 1], len = 7), 
        digits
    )
    legend_labels <- c(
        paste("<", legend_at[1]),
        as.character(legend_at[2:(length(legend_at) - 1)]),
        paste(">", legend_at[length(legend_at)])
    )
    
    plot(raster(matrix(legend_at)),
        legend.only = T, col = pal(length(r_breaks) - 1),
        axis.args = list(at = legend_at, 
        labels = legend_labels, cex.axis = labSize, col.axis = col.axis)
    )
}





