#************************************************************************
# Description: Model performace for all models selected for the manuscript
# Author: Xiaojie(J) Gao
# Date: 2022-01-06
#************************************************************************
rm(list=ls())
source("Code/mod_base.R")
source("Code/mod_manu_models_jags_str.R")


# ~ Model 1 ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~ GPP based
f_model_1 <- FitBayesian(
    model_string = model_1,
    params = model_1_params,
    data = list(
        Y = Y, n = n, gsl = f_gsl, emax = f_emax, emin = f_emin
    ), inits = list(
        beta1 = runif(1, 10, 30), beta2 = runif(1, 0, 20),
        beta3 = runif(1, 10, 30), beta4 = runif(1, 0, 30)
    )
)

# goodness of fit
f_model_1$dic
f_model_1$goodfit <- GetGoodnessOfFit(f_model_1$fit[,2] * 100, Y * 100)


#~ EVI2 based
m_model_1 <- FitBayesian(
    model_string = model_1,
    params = model_1_params,
    data = list(
        Y = Y, n = n, gsl = m_gsl, emax = m_emax, emin = m_emin
    ), inits = list(
        beta1 = runif(1, 10, 30), beta2 = runif(1, 0, 20),
        beta3 = runif(1, 10, 30), beta4 = runif(1, 0, 30)
    )
)
# goodness of fit
m_model_1$dic
m_model_1$goodfit <- GetGoodnessOfFit(m_model_1$fit[, 2] * 100, Y * 100)


#fig: Fit figure for model 1
f_model_1_fit_fig <- PlotFit(f_model_1$fit, north_sites_dt$IGBP, Y, 
    anno_xy = c(300, 3000)
)
m_model_1_fit_fig <- PlotFit(m_model_1$fit, north_sites_dt$IGBP, Y, 
    legendOn = 2, anno_xy = c(300, 3000)
)

f_model_1_beta2 <- quantile(f_model_1$beta2, c(0.025, 0.5, 0.975))
m_model_1_beta2 <- quantile(m_model_1$beta2, c(0.025, 0.5, 0.975))

gsl_model_1_coef <- data.frame(rbind(f_model_1_beta2, m_model_1_beta2))
gsl_model_1_coef$source <- c("Flux GPP-based", "MODIS EVI2-based")
colnames(gsl_model_1_coef) <- c("lwr", "median", "upper", "source")

model_1_beta2_fig <- ggplot(gsl_model_1_coef) +
    geom_pointrange(aes(x = 1:2, ymin = lwr, y = median, ymax = upper, 
        color = source), 
        position = position_dodge(width = 0.1)
    ) +
    xlab("") +
    ylab("GSL Coefficients") +
    theme_article() +
    geom_abline(intercept = 0, slope = 0, linetype = "dashed", 
        color = "grey"
    ) +
    theme(legend.position = c(0.8, 0.9), legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
    )

scatter_plts <- ggarrange(f_model_1_fit_fig, m_model_1_fit_fig, 
    common.legend = TRUE, legend.grob = GetLegend(m_model_1_fit_fig), 
    legend = "top", labels = c("A", "B"))
ggarrange(scatter_plts, model_1_beta2_fig, nrow = 2, labels = c("", "C"))



# ~ Model 2 ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~ GPP based
f_model_2 <- FitBayesian(
    model_string = model_2,
    params = model_2_params,
    data = list(
        Y = Y, n = n, gsl = f_gsl, emax = f_emax, emin = f_emin,
        igbp = f_igbp, nigbp = f_nigbp
    ), inits = list(
        beta1 = runif(f_nigbp, 10, 30), beta2 = runif(1, 0, 20),
        beta3 = runif(1, 10, 30), beta4 = runif(1, 0, 30)
    )
)
# goodness of fit
f_model_2$dic
f_model_2$goodfit <- GetGoodnessOfFit(f_model_2$fit[, 2] * 100, Y * 100)

# ~ EVI2 based
m_model_2 <- FitBayesian(
    model_string = model_2,
    params = model_2_params,
    data = list(
        Y = Y, n = n, gsl = m_gsl, emax = m_emax, emin = m_emin,
        igbp = f_igbp, nigbp = f_nigbp
    ), inits = list(
        beta1 = runif(f_nigbp, 10, 30), beta2 = runif(1, 0, 20),
        beta3 = runif(1, 10, 30), beta4 = runif(1, 0, 30)
    )
)
# goodness of fit
m_model_2$dic
m_model_2$goodfit <- GetGoodnessOfFit(m_model_2$fit[, 2] * 100, Y * 100)


# fig: Fit figure for model 2
f_model_2_fit_fig <- PlotFit(f_model_2$fit, north_sites_dt$IGBP, Y, 
    anno_xy = c(300, 3000))
m_model_2_fit_fig <- PlotFit(m_model_2$fit, north_sites_dt$IGBP, Y, 
    legendOn = 2, anno_xy = c(300, 3000))

f_model_2_beta2 <- quantile(f_model_2$beta2, c(0.025, 0.5, 0.975))
m_model_2_beta2 <- quantile(m_model_2$beta2, c(0.025, 0.5, 0.975))

gsl_model_2_coef <- data.frame(rbind(f_model_2_beta2, m_model_2_beta2))
gsl_model_2_coef$source <- c("Flux GPP-based", "MODIS EVI2-based")
colnames(gsl_model_2_coef) <- c("lwr", "median", "upper", "source")

model_2_beta2_fig <- ggplot(gsl_model_2_coef) +
    geom_pointrange(aes(x = 1:2, ymin = lwr, y = median, ymax = upper, 
        color = source), position = position_dodge(width = 0.1)) +
    xlab("") +
    ylab("GSL Coefficients") +
    theme_article() +
    geom_abline(intercept = 0, slope = 0, linetype = "dashed", 
        color = "grey") +
    theme(
        legend.position = c(0.8, 0.9), legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
    )

scatter_plts <- ggarrange(f_model_2_fit_fig, m_model_2_fit_fig,
    common.legend = TRUE, legend.grob = GetLegend(m_model_2_fit_fig),
    legend = "top", labels = c("A", "B")
)
ggarrange(scatter_plts, model_2_beta2_fig, nrow = 2, labels = c("", "C"))



# ~ Model 3 ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~ GPP based
f_model_3 <- FitBayesian(
    model_string = model_3,
    params = model_3_params,
    data = list(
        Y = Y, n = n, gsl = f_gsl, emax = f_emax, emin = f_emin,
        igbp = f_igbp, nigbp = f_nigbp
    ), inits = list(
        beta1 = runif(m_nigbp, 10, 30), beta2 = runif(m_nigbp, 0, 20),
        beta3 = runif(m_nigbp, 10, 30), beta4 = runif(m_nigbp, 0, 30)
    )
)
f_model_3$dic
f_model_3$goodfit <- GetGoodnessOfFit(f_model_3$fit[, 2] * 100, Y * 100)

# ~ EVI2 based
m_model_3 <- FitBayesian(
    model_string = model_3,
    params = model_3_params,
    data = list(
        Y = Y, n = n, gsl = m_gsl, emax = m_emax, emin = m_emin,
        igbp = f_igbp, nigbp = f_nigbp
    ), inits = list(
        beta1 = runif(f_nigbp, 10, 30), beta2 = runif(f_nigbp, 0, 20),
        beta3 = runif(f_nigbp, 10, 30), beta4 = runif(f_nigbp, 0, 30)
    )
)
m_model_3$dic
m_model_3$goodfit <- GetGoodnessOfFit(m_model_3$fit[, 2] * 100, Y * 100)

# out: biome level models -- model 3
save(f_model_3, m_model_3, file = "Pipeline/model_3_fit_VUT.RData")


# fig: Fit figure for model 3
f_model_3_fit_fig <- PlotFit(f_model_3$fit, north_sites_dt$IGBP, Y, 
    anno_xy = c(300, 3000))
m_model_3_fit_fig <- PlotFit(m_model_3$fit, north_sites_dt$IGBP, Y, 
    legendOn = 2, anno_xy = c(300, 3000))

f_model_3_beta2 <- as.data.frame(apply(
    f_model_3$beta2, 2,
    quantile, c(0.025, 0.5, 0.975)
))
gsl_coef_df <- data.frame(cbind(
    round(as.data.frame(t(f_model_3_beta2)), 2),
    igbp_names
))
gsl_coef_df$source <- "Flux GPP-based"
m_model_3_beta2 <- as.data.frame(apply(
    m_model_3$beta2, 2,
    quantile, c(0.025, 0.5, 0.975)
))
m_model_3_beta2 <- data.frame(cbind(
    round(as.data.frame(t(m_model_3_beta2)), 2),
    igbp_names
))
m_model_3_beta2$source <- "MODIS EVI2-based"
gsl_coef_df <- rbind(gsl_coef_df, m_model_3_beta2)

model_3_beta2_fig <- ggplot(gsl_coef_df) +
    geom_pointrange(aes(x = igbp_names, ymin = X2.5., y = X50., ymax = X97.5., 
        color = source), position = position_dodge(width = 0.5)) +
    xlab("IGBP type") +
    ylab("GSL Coefficients") +
    theme_article() +
    geom_abline(intercept = 0, slope = 0, linetype = "dashed", 
        color = "grey") +
    theme(legend.position = c(0.8, 0.3), legend.title = element_blank()) +
    geom_text(aes(x = IGBP, y = rep(-30, 11), label = N), 
        data = north_sites_dt[, .N, by = "IGBP"]) +
    geom_text(aes(x = IGBP, y = rep(-25, 11), label = N), 
        data = unique(north_sites_dt[, .(site, IGBP)])[, .N, by = c("IGBP")]) +
    scale_y_continuous(breaks = seq(-40, 40, by = 10))

scatter_plts <- ggarrange(f_model_3_fit_fig, m_model_3_fit_fig, 
    common.legend = TRUE, legend.grob = get_legend(m_model_3_fit_fig), 
    legend = "top", labels = c("A", "B")
)
ggarrange(scatter_plts, model_3_beta2_fig, nrow = 2, labels = c("", "C"))



# ~ Model 4 ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~ GPP based
f_model_4 <- FitBayesian(
    model_string = model_4,
    params = model_4_params,
    data = list(
        Y = Y, n = n, gsl = f_gsl, emax = f_emax, emin = f_emin,
        igbp = f_igbp, nigbp = f_nigbp, site = site, ns = ns
    ), inits = list(
        beta1 = runif(ns, 10, 30), beta2 = runif(f_nigbp, 0, 20),
        beta3 = runif(f_nigbp, 10, 30), beta4 = runif(f_nigbp, 0, 30)
    )
)
f_model_4$dic
f_model_4$goodfit <- GetGoodnessOfFit(f_model_4$fit[, 2] * 100, Y * 100)

# ~ EVI2 based
m_model_4 <- FitBayesian(
    model_string = model_4,
    params = model_4_params,
    data = list(
        Y = Y, n = n, gsl = m_gsl, emax = m_emax, emin = m_emin,
        igbp = f_igbp, nigbp = f_nigbp, site = site, ns = ns
    ), inits = list(
        beta1 = runif(ns, 10, 30), beta2 = runif(f_nigbp, 0, 20),
        beta3 = runif(f_nigbp, 10, 30), beta4 = runif(f_nigbp, 0, 30)
    )
)
m_model_4$dic
m_model_4$goodfit <- GetGoodnessOfFit(m_model_4$fit[, 2] * 100, Y * 100)


# fig: Fit figure for model 4
f_model_4_fit_fig <- PlotFit(f_model_4$fit, north_sites_dt$IGBP, Y, 
    anno_xy = c(300, 3000))
m_model_4_fit_fig <- PlotFit(m_model_4$fit, north_sites_dt$IGBP, Y, 
    legendOn = 2, anno_xy = c(300, 3000))

f_model_4_beta2 <- as.data.frame(apply(f_model_4$beta2, 2, quantile, 
    c(0.025, 0.5, 0.975)))
gsl_coef_df <- data.frame(cbind(round(as.data.frame(t(f_model_4_beta2)), 2), 
    igbp_names))
gsl_coef_df$source <- "Flux GPP-based"
m_model_4_beta2 <- as.data.frame(apply(m_model_4$beta2, 2, quantile, 
    c(0.025, 0.5, 0.975)))
m_model_4_beta2 <- data.frame(cbind(round(as.data.frame(t(m_model_4_beta2)), 2), 
    igbp_names))
m_model_4_beta2$source <- "MODIS EVI2-based"
gsl_coef_df <- rbind(gsl_coef_df, m_model_4_beta2)

model_4_beta2_fig <- ggplot(gsl_coef_df) +
    geom_pointrange(aes(x = igbp_names, ymin = X2.5., y = X50., ymax = X97.5., 
        color = source), position = position_dodge(width = 0.5)) +
    xlab("IGBP type") +
    ylab("GSL Coefficients") +
    theme_article() +
    geom_abline(intercept = 0, slope = 0, linetype = "dashed", 
        color = "grey") +
    theme(legend.position = c(0.8, 0.3), legend.title = element_blank()) +
    geom_text(aes(x = IGBP, y = rep(-30, 11), label = N), 
        data = north_sites_dt[, .N, by = "IGBP"]) +
    geom_text(aes(x = IGBP, y = rep(-25, 11), label = N), 
        data = unique(north_sites_dt[, .(site, IGBP)])[, .N, by = c("IGBP")]) +
    scale_y_continuous(breaks = seq(-40, 40, by = 10))

scatter_plts <- ggarrange(f_model_4_fit_fig, m_model_4_fit_fig, 
    common.legend = TRUE, legend.grob = get_legend(m_model_4_fit_fig), 
    legend = "top", labels = c("A", "B")
)
ggarrange(scatter_plts, model_4_beta2_fig, nrow = 2, labels = c("", "C"))




# ~ Model 5 ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~ Unnormalized models
# ~~~~~~~~~~~~~~~~~
f_model <- FitBayesian_model_5(
    data = list(
        Y = Y, site = site, n = n, ns = ns, nigbp = f_nigbp, 
        site_igbp = f_site_igbp,
        emin = f_emin,
        emax = f_emax,
        gsl = f_gsl
    ), inits = list(
        beta1 = runif(ns, 10, 30), beta2 = runif(ns, 0, 20),
        beta3 = runif(ns, 10, 30), beta4 = runif(ns, 0, 30)
    ), site_igbp_names = unique(north_sites_dt[, .(site, IGBP)]),
    igbp_levels = levels(f_igbp_factor)
)

m_model <- FitBayesian_model_5(
    data = list(
        Y = Y, site = site, n = n, ns = ns, nigbp = f_nigbp, 
        site_igbp = f_site_igbp,
        emin = m_emin,
        emax = m_emax,
        gsl = m_gsl
    ), inits = list(
        beta1 = runif(ns, 10, 30), beta2 = runif(ns, 0, 20),
        beta3 = runif(ns, 10, 30), beta4 = runif(ns, 0, 30)
    ), site_igbp_names = unique(north_sites_dt[, .(site, IGBP)]),
    igbp_levels = levels(f_igbp_factor)
)


# ~ Normalised models
# ~~~~~~~~~~~~~~~~~
f_model_norm <- FitBayesian_model_5(
    data = list(
        Y = as.vector(scale(Y)), site = site, n = n, ns = ns, nigbp = f_nigbp,
        site_igbp = f_site_igbp,
        emin = as.vector(scale(f_emin)),
        emax = as.vector(scale(f_emax)),
        gsl = as.vector(scale(f_gsl))
    ), inits = list(
        beta1 = runif(ns, 0, 1), beta2 = runif(ns, 0, 1),
        beta3 = runif(ns, 0, 1), beta4 = runif(ns, 0, 1)
    ), site_igbp_names = unique(north_sites_dt[, .(site, IGBP)]),
    igbp_levels = levels(f_igbp_factor)
)

m_model_norm <- FitBayesian_model_5(
    data = list(
        Y = as.vector(scale(Y)), site = site, n = n, ns = ns, nigbp = f_nigbp, 
        site_igbp = f_site_igbp,
        emin = as.vector(scale(m_emin)),
        emax = as.vector(scale(m_emax)),
        gsl = as.vector(scale(m_gsl))
    ), inits = list(
        beta1 = runif(ns, 0, 1), beta2 = runif(ns, 0, 1),
        beta3 = runif(ns, 0, 1), beta4 = runif(ns, 0, 1)
    ), site_igbp_names = unique(north_sites_dt[, .(site, IGBP)]),
    igbp_levels = levels(f_igbp_factor)
)

# out: best bayesian model fit
save(f_model, m_model, f_model_norm, m_model_norm, 
    file = "Pipeline/best_model_fit_VUT.RData"
)
