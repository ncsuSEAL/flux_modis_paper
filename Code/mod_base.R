#*******************************************************************************
# Description: Store modeling needed data and functions
# 
# Author: Xiaojie Gao
# Date: 2021-06-01
#*******************************************************************************
source("Code/base.R")

library(rjags)
library(egg)
library(MASS)
# library(ggpubr)



#* Read the data in
north_sites_dt <- LoadModelARD(north_only = TRUE)

# ~ convert data to individual vectors
# ~~~~~~~~~~~~~~~~~
# annual GPP, scale it by 100 to make it compatible with predictors' scales
Y <- north_sites_dt$annual_gpp / 100
n <- length(Y) # number of obs.
# flux sites
site_factor <- as.factor(north_sites_dt$site)
site_factor_level <- levels(site_factor)
site <- as.numeric(site_factor)
ns <- length(unique(site)) # number of flux sites

# MODIS EVI2 based predictors and IGBP categories
m_igbp_factor <- as.factor(north_sites_dt$m_IGBP)
m_igbp <- as.numeric(m_igbp_factor)
m_nigbp <- length(unique(m_igbp)) # number of IGBP types
m_gsl <- as.vector(scale(north_sites_dt$m_gsl, scale = FALSE))
m_emax <- as.vector(scale(north_sites_dt$m_EVImax_1, scale = FALSE)) * 10
m_emin <- as.vector(scale(north_sites_dt$m_EVImin_1, scale = FALSE)) * 10

m_site_igbp <- as.numeric(
    as.factor(
        unique(north_sites_dt[, .(site, IGBP = m_IGBP)])$IGBP
    )
)

# Flux GPP based predictors and IGBP categories
f_igbp_factor <- as.factor(north_sites_dt$IGBP)
f_igbp <- as.numeric(f_igbp_factor)
f_nigbp <- length(unique(f_igbp)) # number of IGBP types
f_emin <- as.vector(scale(north_sites_dt$f_gppmin_1, scale = FALSE))
f_emax <- as.vector(scale(north_sites_dt$f_gppmax_1, scale = FALSE))
f_gsl <- as.vector(scale(north_sites_dt$f_gsl, scale = FALSE))

f_site_igbp <- as.numeric(
    as.factor(
        unique(north_sites_dt[, .(site, IGBP)])$IGBP
    )
)

# IGBP names
igbp_names <- levels(as.factor(north_sites_dt$IGBP))
# MODIS IGBP names
m_igbp_names <- levels(as.factor(north_sites_dt$m_IGBP))


# world climate ecoregion
eco <- as.numeric(north_sites_dt$ECO_SYM)
neco <- length(unique(eco)) # number of ecoregions





FitBayesian <- function(model_string, params, data, inits) {
    require(rjags)

    model <- jags.model(textConnection(model_string), 
        data = data, inits = inits, n.chains = 2, quiet = TRUE
    )
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

    return(list(
        beta1 = beta1, beta2 = beta2, beta3 = beta3, beta4 = beta4, 
        mu_beta = mu_beta, tau_beta = tau_beta, tau_y = tau_y, 
        fit = fit, dic = dic
    ))
}

# B/c model 5 will export more variables, I made a separate function for it
FitBayesian_model_5 <- function(data, inits, site_igbp_names, igbp_levels) {
    model_string <- "model {
        # likelihood
        for (i in 1:n) {
            Y[i] ~ dnorm(f[i], tau_y)
            f[i] <- beta1[site[i]] + beta2[site[i]] * gsl[i] + 
                beta3[site[i]] * emax[i] + beta4[site[i]] * emin[i]
        }
    
        # priors
        for (i in 1:ns) {
            beta1[i] ~ dnorm(eta1[site_igbp[i]], tau[1])
            beta2[i] ~ dnorm(eta2[site_igbp[i]], tau[2])
            beta3[i] ~ dnorm(eta3[site_igbp[i]], tau[3])
            beta4[i] ~ dnorm(eta4[site_igbp[i]], tau[4])
        }
        for (i in 1:nigbp) {
            eta1[i] ~ dnorm(mu[1], lambda[1])
            eta2[i] ~ dnorm(mu[2], lambda[2])
            eta3[i] ~ dnorm(mu[3], lambda[3])
            eta4[i] ~ dnorm(mu[4], lambda[4])
        }
        
        tau_y ~ dgamma(0.1, 0.1)
        
        for (i in 1:4) { 
            mu[i] ~ dnorm(0, 0.0001)
            tau[i] ~ dgamma(0.1, 0.1)
            lambda[i] ~ dgamma(0.1, 0.1)
        }

        sigma_y <- 1 / sqrt(tau_y)
        sigma_beta <- 1 / sqrt(tau)
        sigma_eta <- 1 / sqrt(lambda)
    }"

    model <- jags.model(textConnection(model_string), 
        data = data, inits = inits, n.chains = 2, quiet = TRUE
    )
    update(model, 20000, progress.bar = "none")
    iter <- 1
    repeat({
        samp <- coda.samples(model,
            variable.names = c(
                "beta1", "beta2", "beta3", "beta4",
                "eta1", "eta2", "eta3", "eta4",
                "sigma_y", "sigma_beta", "sigma_eta", "mu"
            ),
            n.iter = 50000, thin = 10,
            progress.bar = "none"
        )
        gelman_diag <- gelman.diag(samp)
        mpsrf <- gelman_diag$mpsrf
        iter <- iter + 1
        if (iter >= 6 | mpsrf < 1.1) {
            break
        } else {
            warning("Model might not converged!")
        }
    })

    # browser()

    samp_df <- as.data.frame(rbind(samp[[1]], samp[[2]]))
    samp_names <- colnames(samp_df)

    beta1 <- samp_df[, grep("^beta1", samp_names)] * 100
    beta2 <- samp_df[, grep("^beta2", samp_names)] * 100
    beta3 <- samp_df[, grep("^beta3", samp_names)] * 100
    beta4 <- samp_df[, grep("^beta4", samp_names)] * 100

    eta1 <- samp_df[, grep("^eta1", samp_names)] * 100
    eta2 <- samp_df[, grep("^eta2", samp_names)] * 100
    eta3 <- samp_df[, grep("^eta3", samp_names)] * 100
    eta4 <- samp_df[, grep("^eta4", samp_names)] * 100

    mu <- samp_df[, grep("^mu", samp_names)] * 100
    sigma_y <- samp_df[, grep("^sigma_y", samp_names)] * 100
    sigma_beta <- samp_df[, grep("^sigma_beta", samp_names)] * 100
    sigma_eta <- samp_df[, grep("^sigma_eta", samp_names)] * 100

    beta1_dt <- SummarySample(beta1)
    beta2_dt <- SummarySample(beta2)
    beta3_dt <- SummarySample(beta3)
    beta4_dt <- SummarySample(beta4)

    beta1_dt$site <- site_factor_level
    beta2_dt$site <- site_factor_level
    beta3_dt$site <- site_factor_level
    beta4_dt$site <- site_factor_level

    beta1_dt <- merge(beta1_dt, site_igbp_names, by = "site")
    beta2_dt <- merge(beta2_dt, site_igbp_names, by = "site")
    beta3_dt <- merge(beta3_dt, site_igbp_names, by = "site")
    beta4_dt <- merge(beta4_dt, site_igbp_names, by = "site")

    eta1_dt <- SummarySample(eta1)
    eta2_dt <- SummarySample(eta2)
    eta3_dt <- SummarySample(eta3)
    eta4_dt <- SummarySample(eta4)

    eta1_dt$IGBP <- igbp_levels
    eta2_dt$IGBP <- igbp_levels
    eta3_dt$IGBP <- igbp_levels
    eta4_dt$IGBP <- igbp_levels

    mu_dt <- SummarySample(mu)

    sigma_y_dt <- quantile(sigma_y, c(0.025, 0.5, 0.975))
    sigma_beta_dt <- SummarySample(sigma_beta)
    sigma_eta_dt <- SummarySample(sigma_eta)

    fit_samp <- coda.samples(model,
        variable.names = c("f"),
        n.iter = 50000,
        thin = 10,
        progress.bar = "none"
    )
    fit <- summary(fit_samp)$quantiles[, c(1, 3, 5)]
    colnames(fit) <- c("fit_lwr", "fit_median", "fit_upr")

    return(list(
        beta1 = beta1_dt, beta2 = beta2_dt, beta3 = beta3_dt, beta4 = beta4_dt,
        eta1 = eta1_dt, eta2 = eta2_dt, eta3 = eta3_dt, eta4 = eta4_dt,
        mu = mu_dt, fit = fit,
        sigma_y = sigma_y_dt, sigma_beta = sigma_beta_dt, 
        sigma_eta = sigma_eta_dt
    ))
}


# this function needs global env
PlotFit <- function(fit, IGBP, Y, legendOn = 0, perGroup = FALSE, 
    anno_xy = NULL
) {
    df <- data.frame(fit * 100, agpp = Y * 100, IGBP = as.factor(IGBP))
    
    fit_fig <- ggplot(df) +
        geom_point(aes(x = fit_median, y = agpp, color = IGBP)) +
        xlim(range(df$fit_median, df$agpp)) +
        ylim(range(df$fit_median, df$agpp)) +
        xlab("Fitted") +
        ylab("True") +
        geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
        theme_article()
    
    if (is.null(anno_xy)) { 
        anno_xy <- c(range(df$agpp)[1] + 50, range(df$fit_median)[2]) 
    }
    fit_fig <- fit_fig + annotate("text",
        x = anno_xy[1], y = anno_xy[2],
        label = c(
            paste0("R2:", round(cor(df$fit_median, Y)^2, 2)),
            paste0("\n\nRMSE:", round(
                sqrt(mean((df$fit_median - df$agpp)^2)), 
                2)
            )
        ), hjust = 0
    )

    if (legendOn == 0) {
        fit_fig <- fit_fig + theme(legend.position = "none")
    }
    if (legendOn == 1) {
        fit_fig <- fit_fig + theme(legend.position = "right")
    }
    if (legendOn == 2) {
        fit_fig <- fit_fig + theme(legend.position = "top") + 
            guides(color = guide_legend(ncol = 6))
    }
    if (length(unique(IGBP)) < 11) {
        fit_fig <- fit_fig + scale_color_brewer(palette = "Set3")
    }
    if (perGroup == TRUE) {
        fit_fig <- fit_fig + facet_wrap(~IGBP)
    }
    return(fit_fig)
}


GetGoodnessOfFit <- function(fit_median, true_val) {
    r2 <- round(cor(fit_median, true_val)^2, 2)
    rmse <- round(sqrt(mean((fit_median - true_val)^2)), 2)

    return(list(r2 = r2, rmse = rmse))
}


GetLegend <- function(myggplot) {
    tmp <- ggplot_gtable(ggplot_build(myggplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)
}


SummarySample <- function(samp) {
    samp_dt <- apply(samp, 2, quantile, c(0.025, 0.5, 0.975))
    samp_dt <- data.table(t(samp_dt))
    colnames(samp_dt) <- c("lwr", "med", "upr")
    return(samp_dt)
}
