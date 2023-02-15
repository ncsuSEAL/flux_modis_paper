#************************************************************************
# Description: Compare GPP- and EVI2-based phenometrics
# Author: Xiaojie Gao
# Date: 2021-03-02
#************************************************************************
rm(list=ls())

source("Code/base.R")

# FLux data ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
load("Data/flux_q2algo_GPP_DT_VUT_REF.Rdata") # NOTE

flux <- setDT(fullDT)
flux <- flux[year <= 2014, lapply(.SD, unlist)]
flux$year <- as.character(flux$year)

# print(paste("How many flux sites:", length(unique(flux$site_to_do))))
# print(paste("How many flux site years in total: ", nrow(flux)))

flux_na_sy <- flux[is.na(flux_gppamp_cycle1), .(site_to_do, year)]
# print(paste("How many flux site years are NA: ", nrow(flux_na_sy)))

# rename columns
flux <- flux[, .(
    site_name = site_to_do, year = year, gpp_var = gpp_var,
    f_gpparea_1 = flux_gpparea_cycle1, f_gppamp_1 = flux_gppamp_cycle1, f_gppmin_1 = flux_gppmin_cycle1,
    f_gup_1 = flux_gup_cycle1, f_midgup_1 = flux_midgup_cycle1, f_maturity_1 = flux_maturity_cycle1,
    f_peak_1 = flux_peak_cycle1, f_sene_1 = flux_senescence_cycle1, f_midgdown_1 = flux_midgreendown_cycle1,
    f_dormancy_1 = flux_dormancy_cycle1, f_qa_1 = flux_qa_cycle1, f_qadetailed_1 = flux_qadetailed_cycle1,
    f_gpparea_2 = flux_gpparea_cycle2, f_gppamp_2 = flux_gppamp_cycle2, f_gppmin_2 = flux_gppmin_cycle2,
    f_gup_2 = flux_gup_cycle2, f_midgup_2 = flux_midgup_cycle2, f_maturity_2 = flux_maturity_cycle2,
    f_peak_2 = flux_peak_cycle2, f_sene_2 = flux_senescence_cycle2, f_midgdown_2 = flux_midgreendown_cycle2,
    f_dormancy_2 = flux_dormancy_cycle2, f_qa_2 = flux_qa_cycle2, f_qadetailed_2 = flux_qadetailed_cycle3
)]

na_flux <- flux[is.na(f_gup_1), .(NA_count = .N, year), by = c("site_name")]
flux_always_na <- unique(na_flux[NA_count == 14, .(site_name)])
# print("These sites are alwayis NA in FLUX: ")
# flux_always_na

# delete NA rows and rename columns
flux_noNA <- flux[!is.na(f_gup_1), ]
# print(paste("So, remove NA rows, How many site years left: ", nrow(flux_noNA)))

# how many site years have 2 cycles
flux_two_cycle <- flux_noNA[, .(site_name, year, 
                      twoCycle = ifelse(!is.na(f_gup_1) & !is.na(f_gup_2), TRUE, FALSE))]
# print("How many site years have 2 cycles:")
# flux_two_cycle[twoCycle == TRUE, .(site_name, year)]


flux_miss <- flux[, .(site = as.factor(site_name), year = as.factor(year), 
    missing = ifelse(is.na(f_gup_1), TRUE, FALSE)), by = c("site_name", "year")]


gpp_miss <- fread("Data/gpp_lookup.csv") # NOTE

gpp_miss$year <- as.character(gpp_miss$year)
flux_com <- merge(flux, gpp_miss,
    by.x = c("site_name", "year"),
    by.y = c("siteID", "year"), all = TRUE
)

flux_gpp_but_NA <- flux_com[miss == FALSE & is.na(f_gup_1), .(site_name, year, miss, miss45)]

# print("How many site-years have GPP but no metrics: ")
# nrow(flux_gpp_but_NA)

# print("How many site-years don't have GPP: ")
flux_gpp_missing <- flux_com[miss == TRUE, ]
# nrow(flux_gpp_missing)


site_meta <- LoadFluxSiteMeta()

# ~ MODIS data ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
modis <- LoadMcd12q2() # NOTE

# filter only tier1 sites
modis_tier1 <- modis[!is.na(match(Category, site_meta$siteID)), ]
modis_tier1_dup <- duplicated(modis_tier1[, .(Category, Date)])
modis <- modis_tier1[!modis_tier1_dup, ]
modis$Date <- substr(modis$Date, 5, 8)
# clean the table
modis <- modis[, c(
    "ID", "Latitude", "Longitude", "MODIS_Tile",
    "Line_Y_500m", "Sample_X_500m"
) := NULL]

# rename columns
colnames(modis)[1:23] <- c(
    "site_name", "year", "m_dormancy_1", "m_dormancy_2", "m_EVIamp_1", "m_EVIamp_2", "m_EVIarea_1", "m_EVIarea_2",
    "m_EVImin_1", "m_EVImin_2", "m_gup_1", "m_gup_2", "m_maturity_1", "m_maturity_2", "m_midgdown_1",
    "m_midgdown_2", "m_midgup_1", "m_midgup_2", "NumCycles", "m_peak_1", "m_peak_2", "m_sene_1", "m_sene_2"
)

modis <- modis[year <= 2015, ]

# print(paste("How many sites:", length(unique(modis$site_name))))
# print(paste("How many MODIS site years in total:", nrow(modis)))

modis_na_sy <- modis[is.na(NumCycles), .(site_name, year)]
# print(paste("How many site years are NA: ", nrow(modis_na_sy)))

na_modis <- modis[is.na(NumCycles), .(NA_count = .N, year), by = c("site_name")]
modis_always_na <- unique(na_modis[NA_count == 15, .(site_name)])
# print("These sites are alwayis NA in MODIS: ")
# modis_always_na

modis_noNA <- modis[!is.na(NumCycles), ]
# print(paste("So, remove NA rows, How many site years left: ", nrow(modis_noNA)))

# how many site years have 2 cycles
# modis_noNA[NumCycles == 2, .(site_name, year)]


modis_miss <- modis[, .(site = as.factor(site_name), year = as.factor(year), 
    missing = ifelse(is.na(NumCycles), TRUE, FALSE)),
    by = c("site_name", "year")]



# ~ Compare ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
flux <- flux_com[miss == FALSE, ]
# the final comparison table<U+3001>
rowLength <- nrow(flux)
comp_result_all <- data.table(
    site = flux$site_name,
    year = flux$year, # this year will be the flux year

    f_gup_1 = numeric(rowLength),
    m_gup_1 = numeric(rowLength),
    f_gup_2 = numeric(rowLength),
    m_gup_2 = numeric(rowLength),
    f_midgup_1 = numeric(rowLength),
    m_midgup_1 = numeric(rowLength),
    f_midgup_2 = numeric(rowLength),
    m_midgup_2 = numeric(rowLength),
    f_maturity_1 = numeric(rowLength),
    m_maturity_1 = numeric(rowLength),
    f_maturity_2 = numeric(rowLength),
    m_maturity_2 = numeric(rowLength),
    f_peak_1 = numeric(rowLength),
    m_peak_1 = numeric(rowLength),
    f_peak_2 = numeric(rowLength),
    m_peak_2 = numeric(rowLength),
    f_sene_1 = numeric(rowLength),
    m_sene_1 = numeric(rowLength),
    f_sene_2 = numeric(rowLength),
    m_sene_2 = numeric(rowLength),
    f_midgdown_1 = numeric(rowLength),
    m_midgdown_1 = numeric(rowLength),
    f_midgdown_2 = numeric(rowLength),
    m_midgdown_2 = numeric(rowLength),
    f_dormancy_1 = numeric(rowLength),
    m_dormancy_1 = numeric(rowLength),
    f_dormancy_2 = numeric(rowLength),
    m_dormancy_2 = numeric(rowLength),
    f_gppmin_1 = numeric(rowLength),
    m_EVImin_1 = numeric(rowLength),
    f_gppmin_2 = numeric(rowLength),
    m_EVImin_2 = numeric(rowLength),
    f_gppamp_1 = numeric(rowLength),
    m_EVIamp_1 = numeric(rowLength),
    f_gppamp_2 = numeric(rowLength),
    m_EVIamp_2 = numeric(rowLength),
    f_gpparea_1 = numeric(rowLength),
    m_EVIarea_1 = numeric(rowLength),
    f_gpparea_2 = numeric(rowLength),
    m_EVIarea_2 = numeric(rowLength)
)
comp_result_all[comp_result_all == 0] <- -9999

# this list will record site years that don't match
no_match <- data.table(site = character(), year = numeric(), metric = character())
no_match_modis_na <- data.table(site = character(), year = numeric(), metric = character())
no_match_flux_na <- data.table(site = character(), year = numeric())
no_match_both_na <- data.table(site = character(), year = numeric())
no_match_real <- data.table(site = character(), year = numeric(), metric = character()) # real NO MATCH

# Doing a foor loop for each flux site year, and try to find one that match best in MODIS
for (i in 1:nrow(flux)) { # for each row of flux table
    tryCatch(
        {
            cur <- flux[i, ]

            # check if it is due to modis NA
            cur_m <- modis[site_name == cur$site & year == cur$year, ]
            # check if flux is NA
            cur_f <- flux[site_name == cur$site & year == cur$year, ]

            ifModisNA <- ifelse(is.na(cur_m$NumCycles), TRUE, FALSE)
            ifFluxNA <- ifelse(is.na(cur_f$f_gup_1), TRUE, FALSE)


            if (ifFluxNA) {
                if (ifModisNA & ifFluxNA) { # both flux and modis NA
                    no_match_both_na <- rbind(no_match_both_na, list(cur$site, cur$year))
                } else if (!ifModisNA & ifFluxNA) { # due to flux NA, but modis is not NA
                    no_match_flux_na <- rbind(no_match_flux_na, list(cur$site, cur$year))
                    # here we can only check if flux is NA but MODIS not
                    # we can't be sure if MODIS is NA this year because MODIS store the year as peak year,
                    #   it could be that MODIS is NA this year, but the matched cycle is actually in
                    #   the previouse/next year.
                }
                next # if flux is NA, we can't do any comparison
            }


            f <- flux[i, ]

            phenometrics <- c("gup", "midgup", "maturity", "peak", "sene", "midgdown", "dormancy")

            for (j in phenometrics) { # for all phenometrics
                # construct column names
                m_1 <- paste0("m_", j, "_1")
                m_2 <- paste0("m_", j, "_2")
                f_1 <- paste0("f_", j, "_1")
                f_2 <- paste0("f_", j, "_2")

                f_var_1 <- f[[f_1]]
                f_var_2 <- f[[f_2]]

                # just get all site-years of MODIS, and say find the one with min distance
                m <- modis_noNA[site_name == f$site_name, c("site_name", "year", m_1, m_2),
                    with = FALSE
                ]
                if (nrow(m) == 0) { # modis doesn't have this site in any year
                    # no_match <- rbind(no_match, list(f$site_name, f$year, j))
                    no_match_modis_na <- rbind(no_match_modis_na, list(cur$site, cur$year, j))
                    next
                }

                m_var_1 <- m[[m_1]]
                m_var_2 <- m[[m_2]]

                if (is.na(f_var_2)) { # only the first cycle exist
                    # modis first cycle
                    dis1 <- abs(f_var_1 - m_var_1)
                    min_dis1 <- min(dis1, na.rm = TRUE)
                    # modis second cycle
                    dis2 <- abs(f_var_1 - m_var_2)
                    if (!all(is.na(dis2))) {
                        min_dis2 <- min(dis2, na.rm = TRUE)
                    } else {
                        min_dis2 <- 999999
                    }

                    min_dis <- min(min_dis1, min_dis2)

                    if (min_dis > 185) { # this means no match for both cycles
                        # no_match <- rbind(no_match, list(f$site_name, f$year, j))
                        if (ifModisNA) {
                            no_match_modis_na <- rbind(no_match_modis_na, list(cur$site, cur$year, m_1))
                        } else {
                            no_match_real <- rbind(no_match_real, list(cur$site, cur$year, j))
                        }

                        next # only 1 cycle, and no match with this cycle, so no match for this site year
                    }

                    # if found match
                    if (min_dis1 < min_dis2) { # in first cycle
                        comp_result_all[site == f$site_name & year == f$year][[f_1]] <- f_var_1
                        comp_result_all[site == f$site_name & year == f$year][[m_1]] <- m_var_1[which(min_dis1 == dis1)[1]]
                    } else { # in second cycle
                        comp_result_all[site == f$site_name & year == f$year][[f_1]] <- f_var_1
                        # comp_result_all[site == f$site_name & year == f$year][[m_2]] <- m_var_2[which(min_dis2 == dis2)[1]]
                        comp_result_all[site == f$site_name & year == f$year][[m_1]] <- m_var_2[which(min_dis2 == dis2)[1]]
                    }
                } else { # has two cycles
                    # flux first cycle and modis first cycle
                    dis11 <- abs(f_var_1 - m_var_1)
                    min_dis11 <- min(dis11, na.rm = TRUE)
                    # flux first cycle and modis second cycle
                    dis12 <- abs(f_var_1 - m_var_2)
                    if (!all(is.na(dis12))) {
                        min_dis12 <- min(dis12, na.rm = TRUE)
                    } else {
                        min_dis12 <- 999999
                    }

                    # flux second cycle and modis first cycle
                    dis21 <- abs(f_var_2 - m_var_1)
                    min_dis21 <- min(dis21, na.rm = TRUE)
                    # flux second cycle and modis second cycle
                    dis22 <- abs(f_var_2 - m_var_2)
                    if (!all(is.na(dis22))) {
                        min_dis22 <- min(dis22, na.rm = TRUE)
                    } else {
                        min_dis22 <- 999999
                    }

                    min_dis <- min(min_dis11, min_dis12, min_dis21, min_dis22)

                    if (min_dis > 185) { # this means no match for both cycles
                        # no_match <- rbind(no_match, list(f$site_name, f$year, j))
                        if (ifModisNA) {
                            no_match_modis_na <- rbind(no_match_modis_na, list(cur$site, cur$year, m_2))
                        } else {
                            no_match_real <- rbind(no_match_real, list(cur$site, cur$year, j))
                        }
                        next # only 1 cycle, and no match with this cycle, so no match for this site year
                    }

                    # if found match, decide where to put values
                    if (min_dis == min_dis11) { # flux 1st cycle and modis 1st cycle
                        comp_result_all[site == f$site_name & year == f$year][[f_1]] <- f_var_1
                        comp_result_all[site == f$site_name & year == f$year][[m_1]] <- m_var_1[which(min_dis == dis11)[1]]
                    } else if (min_dis == min_dis12) { # flux 1st cycle and modis 2nd cycle
                        comp_result_all[site == f$site_name & year == f$year][[f_1]] <- f_var_1
                        # comp_result_all[site == f$site_name & year == f$year][[m_2]] <- m_var_2[which(min_dis == dis12)[1]]
                        comp_result_all[site == f$site_name & year == f$year][[m_1]] <- m_var_2[which(min_dis == dis12)[1]]
                    } else if (min_dis == min_dis21) { # flux 2nd cycle and modis 1st cycle
                        comp_result_all[site == f$site_name & year == f$year][[f_2]] <- f_var_2
                        # comp_result_all[site == f$site_name & year == f$year][[m_1]] <- m_var_1[which(min_dis == dis21)[1]]
                        comp_result_all[site == f$site_name & year == f$year][[m_2]] <- m_var_1[which(min_dis == dis21)[1]]
                    } else if (min_dis == min_dis22) { # flux 2nd cycle and modis 2nd cycle
                        comp_result_all[site == f$site_name & year == f$year][[f_2]] <- f_var_2
                        comp_result_all[site == f$site_name & year == f$year][[m_2]] <- m_var_2[which(min_dis == dis22)[1]]
                    }
                }
            }
        },
        error = function(e) {
            stop("something is wrong!")
        }
    )
}
comp_result_all[comp_result_all == -9999] <- NA

comp_result_all <- comp_result_all[!(
    is.na(f_gup_1) & is.na(f_midgup_1) & is.na(f_maturity_1) & is.na(f_peak_1) & is.na(f_sene_1) &
        is.na(f_midgdown_1) & is.na(f_dormancy_1) & is.na(f_gup_2) & is.na(f_midgup_2) & is.na(f_maturity_2) & is.na(f_peak_2) & is.na(f_sene_2) &
        is.na(f_midgdown_2) & is.na(f_dormancy_2)
), ]

no_match_flux_na <- unique(no_match_flux_na)
no_match_modis_na <- no_match_modis_na[, .(.N), by = c("site", "year")][N >= 7, ]
no_match_real <- no_match_real[, .(.N), by = c("site", "year")][N >= 7, ]


# How many site years don't match?
# no_match_site_years <- unique(no_match[, .(site, year)])

# How many site years both gpp and evi amp too low?
gpp_evi_low <- merge(flux_gpp_but_NA, modis_na_sy)

find_match_sy <- nrow(flux) - nrow(no_match_flux_na) - nrow(no_match_modis_na) - nrow(no_match_both_na) - nrow(no_match_real)


getmode <- function(v) {
    uniqv <- unique(v)
    uniqv[which.max(tabulate(match(v, uniqv)))]
}
# fill in min, area, and amplitude
for (i in 1:nrow(comp_result_all)) {
    cur_site <- comp_result_all[i, site]
    cur_year <- comp_result_all[i, year]

    # first cycle
    major_year_1 <- year(as.Date(as.numeric(comp_result_all[i, .(m_gup_1, m_midgup_1, m_maturity_1, m_peak_1, m_sene_1, m_midgdown_1, m_dormancy_1)]), origin = "1970-01-01"))
    major_year_1 <- getmode(major_year_1)

    cur_evi2_vals <- unlist(modis[site_name == cur_site & year == major_year_1, .(m_EVImin_1, m_EVIarea_1, m_EVIamp_1)])
    cur_gpp_vals <- unlist(flux[site_name == cur_site & year == cur_year, .(f_gppmin_1, f_gpparea_1, f_gppamp_1)])
    comp_result_all[i, ":="(m_EVImin_1 = cur_evi2_vals[1], m_EVIarea_1 = cur_evi2_vals[2], m_EVIamp_1 = cur_evi2_vals[3],
        f_gppmin_1 = cur_gpp_vals[1], f_gpparea_1 = cur_gpp_vals[2], f_gppamp_1 = cur_gpp_vals[3])]

    # second cycle
    major_year_2 <- year(as.Date(as.numeric(comp_result_all[i, .(m_gup_2, m_midgup_2, m_maturity_2, m_peak_2, m_sene_2, m_midgdown_2, m_dormancy_2)]), origin = "1970-01-01"))
    major_year_2 <- getmode(major_year_2)

    if (!is.na(major_year_2)) {
        cur_evi2_vals <- unlist(modis[site_name == cur_site & year == major_year_2, .(m_EVImin_2, m_EVIarea_2, m_EVIamp_2)])
        cur_gpp_vals <- unlist(flux[site_name == cur_site & year == cur_year, .(f_gppmin_2, f_gpparea_2, f_gppamp_2)])

        comp_result_all[i, ":="(m_EVImin_2 = cur_evi2_vals[1], m_EVIarea_2 = cur_evi2_vals[2], m_EVIamp_2 = cur_evi2_vals[3],
            f_gppmin_2 = cur_gpp_vals[1], f_gpparea_2 = cur_gpp_vals[2], f_gppamp_2 = cur_gpp_vals[3])]
    }
}
# names(comp_result_all)

# out:`match_result.Rds` and `match_processed.RData`
saveRDS(comp_result_all, file = "Pipeline/match_result.Rds")
save.image(file = "Pipeline/match_processed.RData")