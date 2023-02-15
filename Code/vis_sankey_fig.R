#************************************************************************************
# Description: The flux MODIS comparison sankey diagram
# Author: Xiaojie(J) Gao
# Date: 2022-03-04
#************************************************************************************
rm(list=ls())

source("Code/base.R")
source("Code/mod_base.R")



library(networkD3)

load(file.path("Pipeline/match_processed.RData"))

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
links <- rbind(links, c("Fluxnet2015 (GPP)", "No GPP metrics", 
    nrow(flux[is.na(f_gup_1),])))
links <- rbind(links, c("Fluxnet2015 (GPP)", "GPP metrics", 
    nrow(flux[!is.na(f_gup_1),])))
links <- rbind(links, c("MCD12Q2 (EVI2)", "No EVI2 phenometrics", 
    nrow(gpp_evi_low) + nrow(no_match_modis_na)))
links <- rbind(links, c("MCD12Q2 (EVI2)", "EVI2 phenometrics", 
    nrow(flux) - nrow(gpp_evi_low) - nrow(no_match_modis_na)))

links <- rbind(links, c("No GPP metrics", "MODIS non-NA, but FLux NA", 
    nrow(flux_gpp_but_NA) - nrow(gpp_evi_low)))
links <- rbind(links, c("No GPP metrics", "GPP and EVI2 amplitude both too low", 
    nrow(gpp_evi_low)))
links <- rbind(links, c("GPP metrics", "Full cycle match", nrow(full_match)))
links <- rbind(links, c("GPP metrics", "Part cycle match", nrow(part_match)))
links <- rbind(links, c("GPP metrics", "Both have data but don't match", 
    nrow(no_match_real)))
links <- rbind(links, c("GPP metrics", "Flux non-NA, but MODIS NA", 
    nrow(no_match_modis_na)))

links <- rbind(links, c("No EVI2 phenometrics", "GPP and EVI2 amplitude both too low", 
    nrow(gpp_evi_low)))
links <- rbind(links, c("No EVI2 phenometrics", "Flux non-NA, but MODIS NA", 
    nrow(no_match_modis_na)))
links <- rbind(links, c("EVI2 phenometrics", "Full cycle match", 
    nrow(full_match)))
links <- rbind(links, c("EVI2 phenometrics", "Part cycle match", 
    nrow(part_match)))
links <- rbind(links, c("EVI2 phenometrics", "Both have data but don't match", 
    nrow(no_match_real)))
links <- rbind(links, c("EVI2 phenometrics", "MODIS non-NA, but FLux NA", 
    nrow(no_match_flux_na)))

links <- links[-1,]

nodes <- data.frame(
  name = unique(c(as.character(links$source), as.character(links$target)))
)

nodes$group <- as.factor(c("flux", "modis", "gpp_low", "gpp_metrics", "evi_low", 
    "evi_phenometrics", "flux_na", "both_low", "match", "match", "no_match", 
    "modis_na"))

links$IDsource <- match(links$source, nodes$name) - 1
links$IDtarget <- match(links$target, nodes$name) - 1

links$group <- as.factor(c("gpp_low", "flux", "evi_low", "modis", "gpp_low", 
    "both_low", "match", "match", "no_match","modis_na", "both_low", "evi_low", 
    "match", "match", "no_match", "flux_na"))


my_color <- "d3.scaleOrdinal()
  .domain(['flux', 'modis', 'gpp_low', 'evi_low', 'both_low', 'match', 'no_match',
      'modis_na', 'flux_na'])
  .range(['#7570B3', '#E7298A', 'grey', 'grey', '#1B9E77', '#1B9E77', 'red',
      'grey', 'grey'])"

sankeyNetwork(Links = links, Nodes = nodes,
              Source = "IDsource", Target = "IDtarget",
              Value = "value", NodeID = "name",
              sinksRight=FALSE, fontSize = 14, colourScale = my_color,
              LinkGroup = "group", NodeGroup = "group")
