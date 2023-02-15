#************************************************************************************************************************
# Description: Report modeling results
#  So, the goal is to answer the question: can we use MODIS phenology to 
#  quanfity annual GPP? or, how much carbon dynamics can be explained by land surface phenology?
#  We started by comparing phenology derived from daily GPP observations and LSP from MODIS. We found that 
#  MODIS-LSP matched well with GPP-pheno at IGBP level, but has large mean absolute deviance from GPP-pheno
#  when comparing at site level. Then, inspired by Xia et al 2015, we wanted to see how good MODIS LSP
#  explain annual GPP. This script is to summarize the process of modeling and the main findings from the 
#  models we fit.
# Author: Xiaojie Gao
# Date: 2021-05-18
#************************************************************************************************************************
source("Code/base.R")
source("Code/mod_base.R")

library(ggplot2)
library(rjags)
library(patchwork)
library(egg)



FitBayesian <- function(model_string, params, data, inits) {
    require(rjags)

    model <- jags.model(textConnection(model_string), data = data, inits = inits, n.chains = 2, quiet = TRUE)
    update(model, 2000, progress.bar = "none")
    repeat({
        samp <- coda.samples(model,
            variable.names = params,
            n.iter = 10000,
            thin = 10,
            progress.bar = "none"
        )
        mpsrf <- gelman.diag(samp)$mpsrf

        if (abs(mpsrf) < 1.1 & max(abs(gelman.diag(samp)$psrf)) < 1.1) break
    })
    # max(abs(gelman.diag(samp)$psrf))
    # min(effectiveSize(samp))
    # gelman.diag(samp)


    # retrieve parameters
    samp_df <- as.data.frame(rbind(samp[[1]], samp[[2]]))
    samp_names <- colnames(samp_df)

    beta1 <- samp_df[, grep("beta1", samp_names)] * 100
    beta2 <- samp_df[, grep("beta2", samp_names)] * 100
    beta3 <- samp_df[, grep("beta3", samp_names)] * 100
    beta4 <- samp_df[, grep("beta4", samp_names)] * 100
    mu_beta <- samp_df[, grep("mu_beta", samp_names)]
    tau_beta <- samp_df[, grep("tau_beta", samp_names)]
    tau_y <- samp_df[, grep("tau_y", samp_names)]

    fit_samp <- coda.samples(model,
        variable.names = c("f"),
        n.iter = 5000,
        thin = 10,
        progress.bar = "none"
    )
    fit <- summary(fit_samp)$quantiles[, c(1, 3, 5)]
    colnames(fit) <- c("fit_lwr", "fit_median", "fit_upr")

    dic <- dic.samples(model, n.iter = 5000, progress.bar = "none")

    return(list(beta1 = beta1, beta2 = beta2, beta3 = beta3, beta4 = beta4, mu_beta = mu_beta, 
        tau_beta = tau_beta, tau_y = tau_y, fit = fit, dic = dic))
}

# this function needs global evn
PlotFit <- function(fit, IGBP, legendOn = TRUE, perGroup = FALSE) {
    df <- data.frame(fit = fit * 100, agpp = Y * 100, IGBP = as.factor(IGBP))
    fit_fig <- ggplot(df) +
        geom_point(aes(x = fit.fit_median, y = agpp, color = IGBP)) +
        xlim(range(df$fit.fit_median, df$agpp)) +
        ylim(range(df$fit.fit_median, df$agpp)) +
        xlab("Fitted") +
        ylab("True") +
        geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
        annotate("text", x = range(df$agpp)[1] + 50, y = range(df$fit.fit_median)[2], 
            label = paste0("R2:", round(cor(df$fit.fit_median, Y)^2, 2)), hjust = 0) +
        annotate("text", x = range(df$agpp)[1] + 50, y = range(df$fit.fit_median)[2] - 100, label = paste0("RMSE:", 
            round(sqrt(mean((df$fit.fit_median - df$agpp)^2, na.rm = TRUE)), 2)), hjust = 0)
    if (legendOn == FALSE) fit_fig <- fit_fig + theme(legend.position = "none")
    if (length(unique(IGBP)) < 11) fit_fig <- fit_fig + scale_color_brewer(palette = "Set3")
    if (perGroup == TRUE) fit_fig <- fit_fig + facet_wrap(~IGBP)
    return(fit_fig)
}


#* 1. The simple model means no random effects
#region [The simple model] ####
#************************************************************
model_1 <- "model {
    # likelihood
    for (i in 1:n) {
        Y[i] ~ dnorm(f[i], tau_y)
        f[i] <- beta1 + beta2 * gsl[i] + beta3 * emax[i] + beta4 * emin[i]
    }

    beta1 ~ dnorm(0, 0.0001) # intercept
    beta2 ~ dnorm(0, 0.001) # GSL effect
    beta3 ~ dnorm(0, 0.001) # Zmax effect
    beta4 ~ dnorm(0, 0.001) # Zmin effect
    
    tau_y ~ dgamma(0.1, 0.1)
}"
params <- c("beta1", "beta2", "beta3", "beta4", "tau_y")

# GPP based
f_model_1 <- FitBayesian(
    model_string = model_1,
    params = params,
    data = list(
        Y = Y, n = n, gsl = f_gsl, emax = f_emax, emin = f_emin
    ), inits = list(
        beta1 = runif(1, 10, 30), beta2 = runif(1, 0, 20),
        beta3 = runif(1, 10, 30), beta4 = runif(1, 0, 30)
    )
)
f_model_1$dic
f_model_1_fit_fig <- PlotFit(f_model_1$fit, north_sites_dt$IGBP, legendOn = FALSE)

# EVI2 based
m_model_1 <- FitBayesian(
    model_string = model_1,
    params = params,
    data = list(
        Y = Y, n = n, gsl = m_gsl, emax = m_emax, emin = m_emin
    ), inits = list(
        beta1 = runif(1, 10, 30), beta2 = runif(1, 0, 20),
        beta3 = runif(1, 10, 30), beta4 = runif(1, 0, 30)
    )
)
m_model_1$dic
m_model_1_fit_fig <- PlotFit(m_model_1$fit, north_sites_dt$IGBP, legendOn = FALSE)

f_model_1_fit_fig + m_model_1_fit_fig

#! The GPP-based model fits well, but EVI2-based model doesn't.

#~ What if replacing MODIS LSP with GPP-pheno?
m_f_model_1 <- FitBayesian(
    model_string = model_1,
    params = params,
    data = list(
        Y = Y, n = n, gsl = f_gsl, emax = m_emax, emin = m_emin
    ), inits = list(
        beta1 = runif(1, 10, 30), beta2 = runif(1, 0, 20),
        beta3 = runif(1, 10, 30), beta4 = runif(1, 0, 30)
    )
)
m_f_model_1$dic
m_f_model_1_fit_fig <- PlotFit(m_f_model_1$fit, north_sites_dt$IGBP, legendOn = FALSE)

#endregion [The simple model ]


#* 2. As we know, the GSL effect might vary among biome types, it's reasonable to allow
#* the intercept and the predictors' slopes to vary by biome types
#region [Add IGBP level random effects] ####
#************************************************************

# ~ First, add a varying intercept ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
model_2 <- "model {
    # likelihood
    for (i in 1:n) {
        Y[i] ~ dnorm(f[i], tau_y)
        f[i] <- beta1[igbp[i]] + beta2 * gsl[i] + beta3 * emax[i] + beta4 * emin[i]
    }

    for (i in 1:nigbp) {
        beta1[i] ~ dnorm(mu_beta[1], tau_beta[1])
    }
    beta2 ~ dnorm(0, 0.001) # GSL effect
    beta3 ~ dnorm(0, 0.001) # Zmax effect
    beta4 ~ dnorm(0, 0.001) # Zmin effect

    mu_beta[1] ~ dunif(0, 30)
    for (i in 1:1) { tau_beta[1] ~ dgamma(0.1, 0.1) }
    tau_y ~ dgamma(0.1, 0.1)
}"
params <- c("beta1", "beta2", "beta3", "beta4", "tau_y", "tau_beta", "mu_beta")

# GPP based
f_model_2 <- FitBayesian(
    model_string = model_2,
    params = params,
    data = list(
        Y = Y, n = n, gsl = f_gsl, emax = f_emax, emin = f_emin, 
        igbp = f_igbp, nigbp = f_nigbp
    ), inits = list(
        beta1 = runif(f_nigbp, 10, 30), beta2 = runif(1, 0, 20),
        beta3 = runif(1, 10, 30), beta4 = runif(1, 0, 30)
    )
)
f_model_2$dic
f_model_2_fit_fig <- PlotFit(f_model_2$fit, north_sites_dt$IGBP, legendOn = FALSE)
#! It seems doesn't improve the GPP-based model

#* what if using MODIS IGBP?
f_model_2_m_igbp <- FitBayesian(
    model_string = model_2,
    params = params,
    data = list(
        Y = Y, n = n, gsl = f_gsl, emax = f_emax, emin = f_emin,
        igbp = m_igbp, nigbp = m_nigbp
    ), inits = list(
        beta1 = runif(m_nigbp, 10, 30), beta2 = runif(1, 0, 20),
        beta3 = runif(1, 10, 30), beta4 = runif(1, 0, 30)
    )
)
f_model_2_m_igbp$dic
f_model_2_m_igbp_fit_fig <- PlotFit(f_model_2_m_igbp$fit, north_sites_dt$m_IGBP, legendOn = FALSE)
#! MODIS IGBP improved the model fit significantly by reducing mean and penalized deviance

# EVI2 based
m_model_2 <- FitBayesian(
    model_string = model_2,
    params = params,
    data = list(
        Y = Y, n = n, gsl = m_gsl, emax = m_emax, emin = m_emin,
        igbp = f_igbp, nigbp = f_nigbp
    ), inits = list(
        beta1 = runif(f_nigbp, 10, 30), beta2 = runif(1, 0, 20),
        beta3 = runif(1, 10, 30), beta4 = runif(1, 0, 30)
    )
)
m_model_2$dic
m_model_2_fit_fig <- PlotFit(m_model_2$fit, north_sites_dt$IGBP, legendOn = FALSE)

#* what if using MODIS IGBP?
m_model_2_m_igbp <- FitBayesian(
    model_string = model_2,
    params = params,
    data = list(
        Y = Y, n = n, gsl = m_gsl, emax = m_emax, emin = m_emin,
        igbp = m_igbp, nigbp = m_nigbp
    ), inits = list(
        beta1 = runif(m_nigbp, 10, 30), beta2 = runif(1, 0, 20),
        beta3 = runif(1, 10, 30), beta4 = runif(1, 0, 30)
    )
)
m_model_2_m_igbp$dic
m_model_2_m_igbp_fit_fig <- PlotFit(m_model_2_m_igbp$fit, north_sites_dt$m_IGBP, legendOn = FALSE)

# ! Both GPP-based and EVI2-based model improved by considering varying intercept by biome types.
# ! And, seems MODIS IGBP is better than flux site IGBP, although the improvement is not huge.

# ~ What if replacing MODIS LSP with GPP-pheno?
m_f_model_2 <- FitBayesian(
    model_string = model_2,
    params = params,
    data = list(
        Y = Y, n = n, gsl = f_gsl, emax = m_emax, emin = m_emin,
        igbp = m_igbp, nigbp = m_nigbp
    ), inits = list(
        beta1 = runif(m_nigbp, 10, 30), beta2 = runif(1, 0, 20),
        beta3 = runif(1, 10, 30), beta4 = runif(1, 0, 30)
    )
)
m_f_model_2$dic
m_f_model_2_fit_fig <- PlotFit(m_f_model_2$fit, north_sites_dt$IGBP, legendOn = FALSE)
#! Again, replacing MODIS LSP with GPP-pheno improved the fit significantly.




# ~ Second, add varying slopes ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
model_3 <- "model {
    # likelihood
    for (i in 1:n) {
        Y[i] ~ dnorm(f[i], tau_y)
        f[i] <- beta1[igbp[i]] + beta2[igbp[i]] * gsl[i] + beta3 * emax[i] + beta4 * emin[i]
    }

    for (i in 1:nigbp) {
        beta1[i] ~ eta1 * prcp + eta2 * vpd + eta3 * temp
        beta2[i] ~ dnorm(mu_beta[2], tau_beta[2])
    }
    beta3 ~ dnorm(mu_beta[3], tau_beta[3])
    beta4 ~ dnorm(mu_beta[4], tau_beta[4])

    mu_beta[1] ~ dunif(0, 30) # intercept
    mu_beta[2] ~ dunif(0, 20)   # GSL effect
    mu_beta[3] ~ dunif(0, 30) # EVI2max effect
    mu_beta[4] ~ dunif(0, 30) # EVI2min effect

    for (i in 1:4) { tau_beta[i] ~ dgamma(0.1, 0.1) }
    tau_y ~ dgamma(0.1, 0.1)
}"
params <- c("beta1", "beta2", "beta3", "beta4", "tau_y", "tau_beta", "mu_beta")

# GPP based
f_model_3 <- FitBayesian(
    model_string = model_3,
    params = params,
    data = list(
        Y = Y, n = n, gsl = f_gsl, emax = f_emax, emin = f_emin,
        igbp = f_igbp, nigbp = f_nigbp
    ), inits = list(
        beta1 = runif(f_nigbp, 10, 30), beta2 = runif(f_nigbp, 0, 20),
        beta3 = runif(1, 10, 30), beta4 = runif(1, 0, 30)
    )
)
f_model_3$dic
f_model_3_fit_fig <- PlotFit(f_model_3$fit, north_sites_dt$IGBP, legendOn = FALSE)

# ! Adding varying slope for GSL improved the GPP-based model

#* what if using MODIS IGBP?
f_model_3_m_igbp <- FitBayesian(
    model_string = model_3,
    params = params,
    data = list(
        Y = Y, n = n, gsl = f_gsl, emax = f_emax, emin = f_emin,
        igbp = m_igbp, nigbp = m_nigbp
    ), inits = list(
        beta1 = runif(m_nigbp, 10, 30), beta2 = runif(m_nigbp, 0, 20),
        beta3 = runif(1, 10, 30), beta4 = runif(1, 0, 30)
    )
)
f_model_3_m_igbp$dic
f_model_3_m_igbp_fit_fig <- PlotFit(f_model_3_m_igbp$fit, north_sites_dt$m_IGBP, legendOn = FALSE)
# ! Again, MODIS IGBP improved the model fit significantly by reducing mean and penalized deviance

#* what if using ECO_SYM?
f_model_3_eco <- FitBayesian(
    model_string = model_3,
    params = params,
    data = list(
        Y = Y, n = n, gsl = f_gsl, emax = f_emax, emin = f_emin,
        igbp = eco, nigbp = neco
    ), inits = list(
        beta1 = runif(neco, 10, 30), beta2 = runif(neco, 0, 20),
        beta3 = runif(1, 10, 30), beta4 = runif(1, 0, 30)
    )
)
f_model_3_eco$dic
f_model_3_eco_fit_fig <- PlotFit(f_model_3_eco$fit, north_sites_dt$ECO_SYM, legendOn = FALSE)
# ! Using ecoregion improved the model fit significantly by reducing mean and penalized deviance

# EVI2 based
m_model_3 <- FitBayesian(
    model_string = model_3,
    params = params,
    data = list(
        Y = Y, n = n, gsl = m_gsl, emax = m_emax, emin = m_emin,
        igbp = f_igbp, nigbp = f_nigbp
    ), inits = list(
        beta1 = runif(f_nigbp, 10, 30), beta2 = runif(f_nigbp, 0, 20),
        beta3 = runif(1, 10, 30), beta4 = runif(1, 0, 30)
    )
)
m_model_3$dic
m_model_3_fit_fig <- PlotFit(m_model_3$fit, north_sites_dt$IGBP, legendOn = FALSE)
#! Compare to EVI2-based model without varying slope for GSL, it improved, but not much

#* what if using MODIS IGBP?
m_model_3_m_igbp <- FitBayesian(
    model_string = model_3,
    params = params,
    data = list(
        Y = Y, n = n, gsl = m_gsl, emax = m_emax, emin = m_emin,
        igbp = m_igbp, nigbp = m_nigbp
    ), inits = list(
        beta1 = runif(m_nigbp, 10, 30), beta2 = runif(m_nigbp, 0, 20),
        beta3 = runif(1, 10, 30), beta4 = runif(1, 0, 30)
    )
)
m_model_3_m_igbp$dic
m_model_3_m_igbp_fit_fig <- PlotFit(m_model_3_m_igbp$fit, north_sites_dt$m_IGBP, legendOn = FALSE)
# ! Both GPP-based and EVI2-based model improved by considering varying intercept by biome types.
# ! And, seems MODIS IGBP is better than flux site IGBP, although the improvement is not huge.

#* what if using ECO_SYM?
m_model_3_eco <- FitBayesian(
    model_string = model_3,
    params = params,
    data = list(
        Y = Y, n = n, gsl = m_gsl, emax = m_emax, emin = m_emin,
        igbp = eco, nigbp = neco
    ), inits = list(
        beta1 = runif(neco, 10, 30), beta2 = runif(neco, 0, 20),
        beta3 = runif(1, 10, 30), beta4 = runif(1, 0, 30)
    )
)
m_model_3_eco$dic
m_model_3_eco_fit_fig <- PlotFit(m_model_3_eco$fit, north_sites_dt$ECO_SYM, legendOn = FALSE)

# ~ What if replacing MODIS LSP with GPP-pheno?
m_f_model_3 <- FitBayesian(
    model_string = model_3,
    params = params,
    data = list(
        Y = Y, n = n, gsl = f_gsl, emax = m_emax, emin = m_emin,
        igbp = eco, nigbp = neco
    ), inits = list(
        beta1 = runif(neco, 10, 30), beta2 = runif(neco, 0, 20),
        beta3 = runif(1, 10, 30), beta4 = runif(1, 0, 30)
    )
)
m_f_model_3$dic
m_f_model_3_fit_fig <- PlotFit(m_f_model_3$fit, north_sites_dt$IGBP, legendOn = FALSE)
# ! Again, replacing MODIS LSP with GPP-pheno improved the fit significantly.

# ~ What if use varying intercept among MODIS IGBP, but varying slope for GSL among ecoregion?
model_3_eco <- "model {
    # likelihood
    for (i in 1:n) {
        Y[i] ~ dnorm(f[i], tau_y)
        f[i] <- beta1[igbp[i]] + beta2[eco[i]] * gsl[i] + beta3 * emax[i] + beta4 * emin[i]
    }

    for (i in 1:nigbp) {
        beta1[i] ~ dnorm(mu_beta[1], tau_beta[1])
    }
    for (i in 1:neco) {
        beta2[i] ~ dnorm(mu_beta[2], tau_beta[2])
    }
    beta3 ~ dnorm(mu_beta[3], tau_beta[3])
    beta4 ~ dnorm(mu_beta[4], tau_beta[4])

    mu_beta[1] ~ dunif(0, 30) # intercept
    mu_beta[2] ~ dunif(0, 20)   # GSL effect
    mu_beta[3] ~ dunif(0, 30) # EVI2max effect
    mu_beta[4] ~ dunif(0, 30) # EVI2min effect

    for (i in 1:4) { tau_beta[i] ~ dgamma(0.1, 0.1) }
    tau_y ~ dgamma(0.1, 0.1)
}"
params <- c("beta1", "beta2", "beta3", "beta4", "tau_y", "tau_beta", "mu_beta")


# GPP based
f_model_3_eco4gsl <- FitBayesian(
    model_string = model_3_eco,
    params = params,
    data = list(
        Y = Y, n = n, gsl = f_gsl, emax = f_emax, emin = f_emin,
        igbp = m_igbp, nigbp = m_nigbp, eco = eco, neco = neco
    ), inits = list(
        beta1 = runif(m_nigbp, 10, 30), beta2 = runif(neco, 0, 20),
        beta3 = runif(1, 10, 30), beta4 = runif(1, 0, 30)
    )
)
f_model_3_eco4gsl$dic
# ! Ecoregion varying slope for GSL doesn't improve GPP-based model

# EVI2 based
m_model_3_eco4gsl <- FitBayesian(
    model_string = model_3_eco,
    params = params,
    data = list(
        Y = Y, n = n, gsl = m_gsl, emax = m_emax, emin = m_emin,
        igbp = m_igbp, nigbp = m_nigbp, eco = eco, neco = neco
    ), inits = list(
        beta1 = runif(m_nigbp, 10, 30), beta2 = runif(neco, 0, 20),
        beta3 = runif(1, 10, 30), beta4 = runif(1, 0, 30)
    )
)
m_model_3_eco4gsl$dic
# ! Same, ecoregion varying slope for GSL doesn't improve GPP-based model

# ! So, the section suggests that we should use MODIS IGBP instead of flux site IGBP
# ! Note that although ecoregion improves model fit when used as intercept category, we don't have training data for all
# ! ecoregion types, whcih means we still can't use it to map the world.



# ~ Add random effects for minimum and maximum ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
model_4 <- "model {
    # likelihood
    for (i in 1:n) {
        Y[i] ~ dnorm(f[i], tau_y)
        f[i] <- beta1[igbp[i]] + beta2[igbp[i]] * gsl[i] + beta3[igbp[i]] * emax[i] + beta4[igbp[i]] * emin[i]
    }

    for (i in 1:nigbp) {
        beta1[i] ~ dnorm(mu_beta[1], tau_beta[1])
        beta2[i] ~ dnorm(mu_beta[2], tau_beta[2])
        beta3[i] ~ dnorm(mu_beta[3], tau_beta[3])
        beta4[i] ~ dnorm(mu_beta[4], tau_beta[4])
    }

    mu_beta[1] ~ dunif(0, 30) # intercept
    mu_beta[2] ~ dunif(0, 20)   # GSL effect
    mu_beta[3] ~ dunif(0, 30) # EVI2max effect
    mu_beta[4] ~ dunif(0, 30) # EVI2min effect

    for (i in 1:4) { tau_beta[i] ~ dgamma(0.1, 0.1) }
    tau_y ~ dgamma(0.1, 0.1)
}"
params <- c("beta1", "beta2", "beta3", "beta4", "tau_y", "tau_beta", "mu_beta")



# GPP based
f_model_4 <- FitBayesian(
    model_string = model_4,
    params = params,
    data = list(
        Y = Y, n = n, gsl = f_gsl, emax = f_emax, emin = f_emin,
        igbp = f_igbp, nigbp = f_nigbp
    ), inits = list(
        beta1 = runif(m_nigbp, 10, 30), beta2 = runif(m_nigbp, 0, 20),
        beta3 = runif(m_nigbp, 10, 30), beta4 = runif(m_nigbp, 0, 30)
    )
)
f_model_4$dic
# f_model_4_fit_fig <- PlotFit(f_model_4$fit, north_sites_dt$m_IGBP, legendOn = FALSE)
f_df <- data.frame(fit = f_model_4$fit * 100, agpp = Y * 100, IGBP = as.factor(north_sites_dt$IGBP))
f_model_4_fit_fig <- ggplot(f_df) +
    geom_point(aes(x = fit.fit_median, y = agpp, color = IGBP)) +
    xlim(range(f_df$fit.fit_median, f_df$agpp)) +
    ylim(range(f_df$fit.fit_median, f_df$agpp)) +
    xlab("Fitted") +
    ylab("True") +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    annotate("text",
        x = 300, y = 3000,
        label = c(
            paste0("R2:", round(cor(f_df$fit.fit_median, Y)^2, 2)),
            paste0("\n\nRMSE:", round(sqrt(mean((f_df$fit.fit_median - f_df$agpp)^2)), 2))
        ), hjust = 0
    ) +
    theme_article() +
    theme(legend.position = "none") +
    scale_color_brewer(palette = "Set3")

# EVI2 based
m_model_4 <- FitBayesian(
    model_string = model_4,
    params = params,
    data = list(
        Y = Y, n = n, gsl = m_gsl, emax = m_emax, emin = m_emin,
        igbp = f_igbp, nigbp = f_nigbp
    ), inits = list(
        beta1 = runif(f_nigbp, 10, 30), beta2 = runif(f_nigbp, 0, 20),
        beta3 = runif(f_nigbp, 10, 30), beta4 = runif(f_nigbp, 0, 30)
    )
)
m_model_4$dic
# m_model_4_fit_fig <- PlotFit(m_model_4$fit, north_sites_dt$m_IGBP, legendOn = FALSE)
# m_model_4_fit_fig_perGroup <- PlotFit(m_model_4$fit, north_sites_dt$IGBP, legendOn = FALSE, perGroup = TRUE)
m_df <- data.frame(fit = m_model_4$fit * 100, agpp = Y * 100, IGBP = as.factor(north_sites_dt$IGBP))
m_model_4_fit_fig <- ggplot(m_df) +
    geom_point(aes(x = fit.fit_median, y = agpp, color = IGBP)) +
    xlim(range(m_df$fit.fit_median, m_df$agpp)) +
    ylim(range(m_df$fit.fit_median, m_df$agpp)) +
    xlab("Fitted") +
    ylab("True") +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    annotate("text",
        x = 300, y = 3000,
        label = c(
            paste0("R2:", round(cor(m_df$fit.fit_median, Y)^2, 2)),
            paste0("\n\nRMSE:", round(sqrt(mean((m_df$fit.fit_median - m_df$agpp)^2)), 2))
        ), hjust = 0
    ) +
    theme_article() +
    theme(legend.position = "top") +
    guides(color = guide_legend(ncol = 6)) +
    scale_color_brewer(palette = "Set3")

f_model_4_fit_fig + m_model_4_fit_fig

summary(f_model_4$mu_beta)
summary(m_model_4$mu_beta)

#! I think this is now the best model
#! But, note that cropland doesn't fit well.


# the coef boxplot
igbp_names <- levels(as.factor(north_sites_dt$IGBP))
# GPP based
f_model_4_beta2 <- as.data.frame(apply(f_model_4$beta2, 2, quantile, c(0.025, 0.5, 0.975)))

gsl_coef_df <- data.frame(cbind(round(as.data.frame(t(f_model_4_beta2)), 2), igbp_names))
gsl_coef_df$source <- "Flux GPP-based"

m_model_4_beta2 <- as.data.frame(apply(m_model_4$beta2, 2, quantile, c(0.025, 0.5, 0.975)))
m_model_4_beta2 <- data.frame(cbind(round(as.data.frame(t(m_model_4_beta2)), 2), igbp_names))
m_model_4_beta2$source <- "MODIS EVI2-based"
gsl_coef_df <- rbind(gsl_coef_df, m_model_4_beta2)

model_4_beta2_fig <- ggplot(gsl_coef_df) +
    geom_pointrange(aes(x = igbp_names, ymin = X2.5., y = X50., ymax = X97.5., color = source), position = position_dodge(width = 0.5)) +
    xlab("IGBP type") +
    ylab("GSL Coefficients") +
    theme_article() +
    geom_abline(intercept = 0, slope = 0, linetype = "dashed", color = "grey") +
    theme(legend.position = c(0.8, 0.3), legend.title = element_blank()) +
    geom_text(aes(x = IGBP, y = rep(-30, 11), label = N), data = north_sites_dt[, .N, by = "IGBP"]) +
    geom_text(aes(x = IGBP, y = rep(-25, 11), label = N), data = unique(north_sites_dt[, .(site, IGBP)])[, .N, by = c("IGBP")]) +
    scale_y_continuous(breaks = seq(-40, 40, by = 10))

scatter_plts <- ggarrange(f_model_4_fit_fig, m_model_4_fit_fig, common.legend = TRUE, legend.grob = get_legend(m_model_4_fit_fig), legend = "top", labels = c("A", "B"))
ggarrange(scatter_plts, model_4_beta2_fig, nrow = 2, labels = c("", "C"))




# ~ What if replace MODIS LSP with GPP-pheno?
# ~~~~~~~~~~~~~~~~~
m_f_model_4 <- FitBayesian(
    model_string = model_4,
    params = params,
    data = list(
        Y = Y, n = n, gsl = f_gsl, emax = m_emax, emin = m_emin,
        igbp = m_igbp, nigbp = m_nigbp
    ), inits = list(
        beta1 = runif(m_nigbp, 10, 30), beta2 = runif(m_nigbp, 0, 20),
        beta3 = runif(m_nigbp, 10, 30), beta4 = runif(m_nigbp, 0, 30)
    )
)
m_f_model_4$dic
m_f_model_4_fit_fig <- PlotFit(m_f_model_4$fit, north_sites_dt$m_IGBP, legendOn = FALSE)

m_f_model_4_fit_fig + m_f_model_3_fit_fig

f_model_4_fit_fig + m_model_4_fit_fig + m_f_model_4_fit_fig

m_f_model_4_fit_fig_perGroup <- PlotFit(m_f_model_4$fit, north_sites_dt$m_IGBP, legendOn = FALSE, perGroup = TRUE)

#endregion [Add IGBP level random effects]



#* Considering site level variablity would significantly increase the model
#region [site-level variability] ####
#************************************************************
model_5 <- "model {
    # likelihood
    for (i in 1:n) {
        Y[i] ~ dnorm(f[i], tau_y)
        f[i] <- beta1[site[i]] + beta2[igbp[i]] * gsl[i] + beta3[igbp[i]] * emax[i] + beta4[igbp[i]] * emin[i]
    }

    for (i in 1:nigbp) {
        beta2[i] ~ dnorm(mu_beta[2], tau_beta[2])
        beta3[i] ~ dnorm(mu_beta[3], tau_beta[3])
        beta4[i] ~ dnorm(mu_beta[4], tau_beta[4])
    }
    for (i in 1:ns) {
        beta1[i] ~ dnorm(mu_beta[1], tau_beta[1])
    }

    mu_beta[1] ~ dunif(0, 30) # intercept
    mu_beta[2] ~ dunif(0, 20)   # GSL effect
    mu_beta[3] ~ dunif(0, 30) # EVI2max effect
    mu_beta[4] ~ dunif(0, 30) # EVI2min effect

    for (i in 1:4) { tau_beta[i] ~ dgamma(0.1, 0.1) }
    tau_y ~ dgamma(0.1, 0.1)
}"
params <- c("beta1", "beta2", "beta3", "beta4", "tau_y", "tau_beta", "mu_beta")


# GPP based
f_model_5 <- FitBayesian(
    model_string = model_5,
    params = params,
    data = list(
        Y = Y, n = n, gsl = f_gsl, emax = f_emax, emin = f_emin,
        igbp = m_igbp, nigbp = m_nigbp, site = site, ns = ns
    ), inits = list(
        beta1 = runif(ns, 10, 30), beta2 = runif(m_nigbp, 0, 20),
        beta3 = runif(m_nigbp, 10, 30), beta4 = runif(m_nigbp, 0, 30)
    )
)
f_model_5$dic
f_model_5_fit_fig <- PlotFit(f_model_5$fit, north_sites_dt$m_IGBP, legendOn = FALSE)

# EVI2 based
m_model_5 <- FitBayesian(
    model_string = model_5,
    params = params,
    data = list(
        Y = Y, n = n, gsl = m_gsl, emax = m_emax, emin = m_emin,
        igbp = m_igbp, nigbp = m_nigbp, site = site, ns = ns
    ), inits = list(
        beta1 = runif(ns, 10, 30), beta2 = runif(m_nigbp, 0, 20),
        beta3 = runif(m_nigbp, 10, 30), beta4 = runif(m_nigbp, 0, 30)
    )
)
m_model_5$dic
m_model_5_fit_fig <- PlotFit(m_model_5$fit, north_sites_dt$m_IGBP, legendOn = FALSE)
m_model_5_fit_fig_perGroup <- PlotFit(m_model_5$fit, north_sites_dt$m_IGBP, legendOn = FALSE, perGroup = TRUE)

summary(m_model_5$mu_beta)
summary(f_model_5$mu_beta)



# ~ What if replace MODIS LSP with GPP-pheno?
# ~~~~~~~~~~~~~~~~~
m_f_model_5 <- FitBayesian(
    model_string = model_5,
    params = params,
    data = list(
        Y = Y, n = n, gsl = f_gsl, emax = m_emax, emin = m_emin,
        igbp = m_igbp, nigbp = m_nigbp, site = site, ns = ns
    ), inits = list(
         beta1 = runif(ns, 10, 30), beta2 = runif(m_nigbp, 0, 20),
        beta3 = runif(m_nigbp, 10, 30), beta4 = runif(m_nigbp, 0, 30)
    )
)
m_f_model_5$dic
m_f_model_5_fit_fig <- PlotFit(m_f_model_5$fit, north_sites_dt$m_IGBP, legendOn = FALSE)
m_f_model_5_fit_fig_perGroup <- PlotFit(m_f_model_5$fit, north_sites_dt$m_IGBP, legendOn = FALSE, perGroup = TRUE)

#endregion [site-level variability]




# * Now, after the above model analyses, I think the model we prefer to go is model 4 and maybe model 5.
# * Model 4 makes sense, as we have global IGBP map but we don't have flux site every where. However, model 4 has a 
# * problem that cropland doesn't fit well. From the fit of model 5, we know that considering site-level variability 
# * addressed the issue of cropland in model 4. But, now we face a new problem. How much does the GSL / Emin / Emax effects
# * vary from model 4 to model 5. If they don't vary much, we can pick model 5 and map out the GSL / Emin / Emax effects for 
# * the whole North Hemisphere, which is the best case. But, if their coefficients vary much, we'd have to go with model 4 and 
# * and find a way to address the cropland issue.
#region [Draw inference] ####
#************************************************************
#$ To run this section, go to above section and fit model 4 and model 5 for both GPP- and EVI2-based first!

# ~ Model 4
# ~~~~~~~~~~~~~~~~~
igbp_names <- levels(as.factor(north_sites_dt$m_IGBP))

# GPP based
f_model_4_beta1 <- as.data.frame(apply(f_model_4$beta1, 2, quantile, c(0.025, 0.5, 0.975)))
f_model_4_beta2 <- as.data.frame(apply(f_model_4$beta2, 2, quantile, c(0.025, 0.5, 0.975)))
f_model_4_beta3 <- as.data.frame(apply(f_model_4$beta3, 2, quantile, c(0.025, 0.5, 0.975)))
f_model_4_beta4 <- as.data.frame(apply(f_model_4$beta4, 2, quantile, c(0.025, 0.5, 0.975)))

f_model_4_beta1_fig <- ggplot(data.frame(cbind(round(as.data.frame(t(f_model_4_beta1)), 2), igbp_names))) +
    geom_pointrange(aes(x = igbp_names, ymin = X2.5., y = X50., ymax = X97.5.)) +
    xlab("") + ylab("Intercept") +
    ggtitle("GPP-based, varying intercept among IGBP")
f_model_4_beta2_fig <- ggplot(data.frame(cbind(round(as.data.frame(t(f_model_4_beta2)), 2), igbp_names))) +
    geom_pointrange(aes(x = igbp_names, ymin = X2.5., y = X50., ymax = X97.5.)) +
    xlab("") + ylab("GSL")
f_model_4_beta3_fig <- ggplot(data.frame(cbind(round(as.data.frame(t(f_model_4_beta3)), 2), igbp_names))) +
    geom_pointrange(aes(x = igbp_names, ymin = X2.5., y = X50., ymax = X97.5.)) +
    xlab("") + ylab("GPPmax")
f_model_4_beta4_fig <- ggplot(data.frame(cbind(round(as.data.frame(t(f_model_4_beta4)), 2), igbp_names))) +
    geom_pointrange(aes(x = igbp_names, ymin = X2.5., y = X50., ymax = X97.5.)) +
    xlab("MODIS IGBP") + ylab("GPPmin")

f_model_4_beta1_fig / f_model_4_beta2_fig / f_model_4_beta3_fig / f_model_4_beta4_fig


# EVI2 based
m_model_4_beta1 <- as.data.frame(apply(m_model_4$beta1, 2, quantile, c(0.025, 0.5, 0.975)))
m_model_4_beta2 <- as.data.frame(apply(m_model_4$beta2, 2, quantile, c(0.025, 0.5, 0.975)))
m_model_4_beta3 <- as.data.frame(apply(m_model_4$beta3, 2, quantile, c(0.025, 0.5, 0.975)))
m_model_4_beta4 <- as.data.frame(apply(m_model_4$beta4, 2, quantile, c(0.025, 0.5, 0.975)))

m_model_4_beta1_fig <- ggplot(data.frame(cbind(round(as.data.frame(t(m_model_4_beta1)), 2), igbp_names))) +
    geom_pointrange(aes(x = igbp_names, ymin = X2.5., y = X50., ymax = X97.5.)) +
    xlab("") + ylab("Intercept") +
    ggtitle("EVI2-based, varying intercept among IGBP")
m_model_4_beta2_fig <- ggplot(data.frame(cbind(round(as.data.frame(t(m_model_4_beta2)), 2), igbp_names))) +
    geom_pointrange(aes(x = igbp_names, ymin = X2.5., y = X50., ymax = X97.5.)) +
    xlab("") + ylab("GSL")
m_model_4_beta3_fig <- ggplot(data.frame(cbind(round(as.data.frame(t(m_model_4_beta3)), 2), igbp_names))) +
    geom_pointrange(aes(x = igbp_names, ymin = X2.5., y = X50., ymax = X97.5.)) +
    xlab("") + ylab("EVI2max")
m_model_4_beta4_fig <- ggplot(data.frame(cbind(round(as.data.frame(t(m_model_4_beta4)), 2), igbp_names))) +
    geom_pointrange(aes(x = igbp_names, ymin = X2.5., y = X50., ymax = X97.5.)) +
    xlab("MODIS IGBP") + ylab("EVI2min")

m_model_4_beta1_fig / m_model_4_beta2_fig / m_model_4_beta3_fig / m_model_4_beta4_fig


# ~ Model 5
# ~~~~~~~~~~~~~~~~~
# GPP based
f_model_5_beta1 <- as.data.frame(apply(f_model_5$beta1, 2, quantile, c(0.025, 0.5, 0.975)))
f_model_5_beta2 <- as.data.frame(apply(f_model_5$beta2, 2, quantile, c(0.025, 0.5, 0.975)))
f_model_5_beta3 <- as.data.frame(apply(f_model_5$beta3, 2, quantile, c(0.025, 0.5, 0.975)))
f_model_5_beta4 <- as.data.frame(apply(f_model_5$beta4, 2, quantile, c(0.025, 0.5, 0.975)))

f_model_5_beta1_fig <- ggplot(data.frame(cbind(round(as.data.frame(t(f_model_5_beta1)), 2), unique(north_sites_dt[, .(site, m_IGBP)])))) +
    geom_pointrange(aes(x = 1:ns, ymin = X2.5., y = X50., ymax = X97.5.)) +
    xlab("") + ylab("Intercept") +
    ggtitle("GPP-based, varying intercept among sites")
f_model_5_beta2_fig <- ggplot(data.frame(cbind(round(as.data.frame(t(f_model_5_beta2)), 2), igbp_names))) +
    geom_pointrange(aes(x = igbp_names, ymin = X2.5., y = X50., ymax = X97.5.)) +
    xlab("") + ylab("GSL")
f_model_5_beta3_fig <- ggplot(data.frame(cbind(round(as.data.frame(t(f_model_5_beta3)), 2), igbp_names))) +
    geom_pointrange(aes(x = igbp_names, ymin = X2.5., y = X50., ymax = X97.5.)) +
    xlab("") + ylab("GPPmax")
f_model_5_beta4_fig <- ggplot(data.frame(cbind(round(as.data.frame(t(f_model_5_beta4)), 2), igbp_names))) +
    geom_pointrange(aes(x = igbp_names, ymin = X2.5., y = X50., ymax = X97.5.)) +
    xlab("MODIS IGBP") + ylab("GPPmin")

f_model_5_beta1_fig / f_model_5_beta2_fig / f_model_5_beta3_fig / f_model_5_beta4_fig


# EVI2 based
m_model_5_beta1 <- as.data.frame(apply(m_model_5$beta1, 2, quantile, c(0.025, 0.5, 0.975)))
m_model_5_beta2 <- as.data.frame(apply(m_model_5$beta2, 2, quantile, c(0.025, 0.5, 0.975)))
m_model_5_beta3 <- as.data.frame(apply(m_model_5$beta3, 2, quantile, c(0.025, 0.5, 0.975)))
m_model_5_beta4 <- as.data.frame(apply(m_model_5$beta4, 2, quantile, c(0.025, 0.5, 0.975)))

m_model_5_beta1_fig <- ggplot(data.frame(cbind(round(as.data.frame(t(m_model_5_beta1)), 2), unique(north_sites_dt[, .(site, m_IGBP)])))) +
    geom_pointrange(aes(x = 1:ns, ymin = X2.5., y = X50., ymax = X97.5.)) +
    xlab("") + ylab("Intercept") +
    ggtitle("EVI2-based, varying intercept among sites")
m_model_5_beta2_fig <- ggplot(data.frame(cbind(round(as.data.frame(t(m_model_5_beta2)), 2), igbp_names))) +
    geom_pointrange(aes(x = igbp_names, ymin = X2.5., y = X50., ymax = X97.5.)) +
    xlab("") + ylab("GSL")
m_model_5_beta3_fig <- ggplot(data.frame(cbind(round(as.data.frame(t(m_model_5_beta3)), 2), igbp_names))) +
    geom_pointrange(aes(x = igbp_names, ymin = X2.5., y = X50., ymax = X97.5.)) +
    xlab("") + ylab("EVI2max")
m_model_5_beta4_fig <- ggplot(data.frame(cbind(round(as.data.frame(t(m_model_5_beta4)), 2), igbp_names))) +
    geom_pointrange(aes(x = igbp_names, ymin = X2.5., y = X50., ymax = X97.5.)) +
    xlab("MODIS IGBP") + ylab("EVI2min")

m_model_5_beta1_fig / m_model_5_beta2_fig / m_model_5_beta3_fig / m_model_5_beta4_fig

# colnames(m_model_4_beta2) <- igbp_names
# saveRDS(m_model_4_beta2, file = "m_model_4_beta2.Rds")

#endregion [Draw inference]

