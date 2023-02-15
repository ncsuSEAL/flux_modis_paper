library(rjags)
library(patchwork)

source("Code/base.R")
source("Code/mod_base.R")


# ~ Fit Bayesian heriarchical model ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
theBayesianModelGSL <- function(data, inits, site_igbp_names, igbp_levels) {
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

    model <- jags.model(textConnection(model_string), data = data, 
        inits = inits, n.chains = 2, quiet = TRUE
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



#region [Overall] ####
#************************************************************
# north_sites_dt[IGBP != "CSH", .N, by = "IGBP"]
# num_site_years <- north_sites_dt[IGBP != "CSH", .N, by = c("site", "IGBP")]
# num_site_years[N > 5, .N] / nrow(num_site_years)
# ! So, about a half of sites we can't make any inference from if doing classical regression


# ~ Unnormalized models
# ~~~~~~~~~~~~~~~~~
f_model <- theBayesianModelGSL(
    data = list(Y = Y, site = site, n = n, ns = ns, nigbp = f_nigbp, 
        site_igbp = f_site_igbp, emin = f_emin, emax = f_emax, gsl = f_gsl
    ), 
    inits = list(beta1 = runif(ns, 10, 30), beta2 = runif(ns, 0, 20), 
        beta3 = runif(ns, 10, 30), beta4 = runif(ns, 0, 30)
    ), 
    site_igbp_names = unique(north_sites_dt[, .(site, IGBP)]), 
    igbp_levels = levels(f_igbp_factor)
)

m_model <- theBayesianModelGSL(
    data = list(Y = Y, site = site, n = n, ns = ns, nigbp = f_nigbp, 
        site_igbp = f_site_igbp, emin = m_emin, emax = m_emax, gsl = m_gsl
    ), 
    inits = list(beta1 = runif(ns, 10, 30), beta2 = runif(ns, 0, 20), 
        beta3 = runif(ns, 10, 30), beta4 = runif(ns, 0, 30)
    ), 
    site_igbp_names = unique(north_sites_dt[, .(site, IGBP)]), 
    igbp_levels = levels(f_igbp_factor)
)


# ~ Normalised models
# ~~~~~~~~~~~~~~~~~
f_model_norm <- theBayesianModelGSL(
    data = list(Y = as.vector(scale(Y)), site = site, n = n, ns = ns, 
        nigbp = f_nigbp, site_igbp = f_site_igbp,
        emin = as.vector(scale(f_emin)),
        emax = as.vector(scale(f_emax)),
        gsl = as.vector(scale(f_gsl))
    ), 
    inits = list(beta1 = runif(ns, 0, 1), beta2 = runif(ns, 0, 1), 
        beta3 = runif(ns, 0, 1), beta4 = runif(ns, 0, 1)
    ), 
    site_igbp_names = unique(north_sites_dt[, .(site, IGBP)]), 
    igbp_levels = levels(f_igbp_factor)
)

m_model_norm <- theBayesianModelGSL(
    data = list(
        Y = as.vector(scale(Y)), site = site, n = n, ns = ns, nigbp = f_nigbp,
        site_igbp = f_site_igbp,
        emin = as.vector(scale(m_emin)),
        emax = as.vector(scale(m_emax)),
        gsl = as.vector(scale(m_gsl))
    ), 
    inits = list(
        beta1 = runif(ns, 0, 1), beta2 = runif(ns, 0, 1),
        beta3 = runif(ns, 0, 1), beta4 = runif(ns, 0, 1)
    ), 
    site_igbp_names = unique(north_sites_dt[, .(site, IGBP)]),
    igbp_levels = levels(f_igbp_factor)
)

# out: best bayesian model fit
save(f_model, m_model, f_model_norm, m_model_norm, 
    file = "best_model_fit.RData"
)


#fig: growing season length effect
ggplot() +
    geom_pointrange(aes(x = IGBP, ymin = lwr, y = med, ymax = upr),
        color = "brown2", size = 0.3, data = f_model$beta2,
        position = position_dodge2(width = 0.3)
    ) +
    geom_pointrange(aes(x = IGBP, ymin = lwr, y = med, ymax = upr),
        color = "blue4", size = 0.8, data = f_model$eta2
    ) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
    geom_pointrange(aes(x = IGBP, ymin = lwr, y = med, ymax = upr),
        color = "#1b9969", size = 0.3, data = f_model$beta3,
        position = position_dodge2(width = 0.3)
    ) +
    geom_pointrange(aes(x = IGBP, ymin = lwr, y = med, ymax = upr),
        color = "#0f746e", size = 0.8, data = f_model$eta3
    ) +
    geom_pointrange(aes(x = IGBP, ymin = lwr, y = med, ymax = upr),
        color = "#bfd460", size = 0.3, data = f_model$beta4,
        position = position_dodge2(width = 0.3)
    ) +
    geom_pointrange(aes(x = IGBP, ymin = lwr, y = med, ymax = upr), 
        color = "#a72fac", size = 0.8, data = f_model$eta4
    ) +
    ylab(expression("GSL Coefficient" ~ (g ~ Cm^"-2" ~ yr^"-2"))) +
    geom_text(aes(x = IGBP, y = rep(-18, 11), label = N), 
        data = north_sites_dt[, .N, by = "IGBP"]
    ) +
    geom_text(aes(x = IGBP, y = rep(-12, 11), label = N), 
        data = unique(north_sites_dt[, .(site, IGBP)])[, .N, by = c("IGBP")]
    ) +
    theme_article() +
    coord_flip()


gsl_coef_fig



ggplot() +
    geom_pointrange(aes(x = 1, ymin = lwr, y = med, ymax = upr),
        color = "brown2", size = 0.3, data = f_model$mu[2,],
        position = position_dodge2(width = 0.3))






gpp_based_model_m_igbp <- theBayesianModelGSL(
    data = list(
        Y = Y, site = site, n = n, ns = ns, nigbp = m_nigbp, 
        site_igbp = m_site_igbp,
        emin = f_emin,
        emax = f_emax,
        gsl = f_gsl
    ), 
    inits = list(
        beta1 = runif(ns, 10, 30), beta2 = runif(ns, 0, 20),
        beta3 = runif(ns, 10, 30), beta4 = runif(ns, 0, 30)
    )
)





sapply(igbp_names, function(x) {
    sd(g_beta2_1[IGBP == x, med])
})

g_beta2_2 <- data.table(gpp_based_model_m_igbp$beta2)
g_beta2_2$site <- site_factor_level
g_beta2_2 <- merge(g_beta2_2, unique(north_sites_dt[, .(site, IGBP = m_IGBP)]), 
    by = "site"
)

sapply(m_igbp_names, function(x) {
    sd(g_beta2_2[IGBP == x, med])
})




# ~ Check onvergence
# ~~~~~~~~~~~~~~~~~
plot(gpp_based_model$mu_samp)
autocorr(gpp_based_model$mu_samp[[1]], lag = 1)
effectiveSize(gpp_based_model$mu_samp)
gelman.diag(gpp_based_model$mu_samp) # greater than 1.1 indicates poor convergence
geweke.diag(gpp_based_model$mu_samp) # |z| greater than 2 indicates poor convergence

plot(gpp_based_model$eta_samp)
effectiveSize(gpp_based_model$eta_samp)
gelman.diag(gpp_based_model$eta_samp) # greater than 1.1 indicates poor convergence
geweke.diag(gpp_based_model$eta_samp) # |z| greater than 2 indicates poor convergence

g_beta_df <- data.table(y = gpp_based_model$beta_quan$beta3, unique(north_sites_dt[IGBP != "CSH", .(site, IGBP)]))
colnames(g_beta_df) <- c("gsl_min", "gsl_median", "gsl_max", "Site", "IGBP")
setorder(g_beta_df, cols = IGBP)
setkey(g_beta_df, "IGBP")

g_beta <- ggplot(g_beta_df) +
    geom_pointrange(aes(x = 1:ns, y = gsl_median, ymax = gsl_max, ymin = gsl_min, color = IGBP)) +
    xlab("Site") +
    ylab("GSL coefficient") +
    ggtitle("GPP-based model")
g_eta <- ggplot(data.frame(x = levels(as.factor(north_sites_dt[IGBP != "CSH"]$IGBP)), y = gpp_based_model$eta_quan$eta3)) +
    geom_pointrange(aes(x = x, y = y.50., ymax = y.97.5., ymin = y.2.5.)) +
    xlab("IGBP") +
    ylab("GSL coefficient") +
    ggtitle("GPP-based model")
gpp_based_model$mu_quan

g_fit <- summary(gpp_based_model$fit_samp)$quantiles[1:n, 3]
g_fit_fig <- ggplot(data.frame(fit = g_fit, agpp = Y, IGBP = north_sites_dt[IGBP != "CSH"]$IGBP)) +
    geom_point(aes(x = fit, y = agpp, color = IGBP)) +
    xlim(0, 3500) +
    ylim(0, 3500) +
    xlab("Fitted") +
    ylab("True") +
    scale_color_brewer(palette = "Set3") +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    ggtitle("GPP-based model")

# Plot fitted for each site
fitted_dt <- data.table(fit = g_fit, agpp = Y, site = north_sites_dt[IGBP != "CSH"]$site, igbp = north_sites_dt[IGBP != "CSH"]$IGBP)
st_plots <- list()
colors <- brewer.pal(10, "Paired")
k <- 1
for (i in unique(fitted_dt$igbp)) {
    point_color <- colors[k]
    k <- k + 1
    for (j in unique(fitted_dt[igbp == i]$site)) {
        fig <- ggplot(fitted_dt[igbp == i & site == j, ]) +
            geom_point(aes(x = fit, y = agpp), color = point_color) +
            xlim(0, 3500) +
            ylim(0, 3500) +
            xlab("Predicted") +
            ylab("True") +
            scale_color_brewer(palette = "Set3") +
            geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
            ggtitle(paste(j, "(", i, ")")) +
            theme(legend.position = "none")

        st_plots <- append(st_plots, list(fig))
    }
}
pdf("GPP_based_Bayesian_model_per_site.pdf", width = 40, height = 35)
wrap_plots(st_plots, nrow = 10)
dev.off()



# ~ EVI2-based
# ~~~~~~~~~~~~~~~~~
evi2_based_model <- theBayesianModelGSL(data = list(
    Y = Y, site = site, igbp = f_igbp, n = n, ns = ns, nigbp = f_nigbp, site_igbp = f_ site_igbp,
    emin = m_emin,
    emax = m_emax,
    gsl = m_gsl
),  inits = list(beta0 = runif(ns, 10, 30), beta1 = runif(ns, 0, 20), beta2 = runif(ns, 10, 30), beta3 = runif(ns, 0, 30)))

e_beta_df <- data.table(y = evi2_based_model$beta_quan$beta3, unique(north_sites_dt[IGBP != "CSH", .(site, IGBP)]))
colnames(e_beta_df) <- c("gsl_min", "gsl_median", "gsl_max", "Site", "IGBP")
setorder(e_beta_df, cols = IGBP)
setkey(e_beta_df, "IGBP")

e_beta <- ggplot(e_beta_df) +
    geom_pointrange(aes(x = 1:ns, y = gsl_median, ymax = gsl_max, ymin = gsl_min, color = IGBP)) +
    xlab("Site") +
    ylab("GSL coefficient") +
    ggtitle("EVI2-based model")
e_eta <- ggplot(data.frame(x = levels(as.factor(north_sites_dt[IGBP != "CSH"]$IGBP)), y = evi2_based_model$eta_quan$eta3)) +
    geom_pointrange(aes(x = x, y = y.50., ymax = y.97.5., ymin = y.2.5.)) +
    xlab("IGBP") +
    ylab("GSL coefficient") +
    ggtitle("EVI2-based model")
evi2_based_model$mu_quan

e_fit <- summary(evi2_based_model$fit_samp)$quantiles[1:n, 3]
e_fit_fig <- ggplot(data.frame(fit = e_fit, agpp = Y, IGBP = north_sites_dt[IGBP != "CSH"]$IGBP)) +
    geom_point(aes(x = fit, y = agpp, color = IGBP)) +
    xlim(0, 3500) +
    ylim(0, 3500) +
    xlab("Fitted") +
    ylab("True") +
    scale_color_brewer(palette = "Set3") +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    ggtitle("EVI2-based model")


# Plot fitted for each site
fitted_dt <- data.table(fit = e_fit, agpp = Y, site = north_sites_dt[IGBP != "CSH"]$site, igbp = north_sites_dt[IGBP != "CSH"]$IGBP)
st_plots <- list()
colors <- brewer.pal(10, "Paired")
k <- 1
for (i in unique(fitted_dt$igbp)) {
    point_color <- colors[k]
    k <- k + 1
    for (j in unique(fitted_dt[igbp == i]$site)) {
        fig <- ggplot(fitted_dt[igbp == i & site == j, ]) +
            geom_point(aes(x = fit, y = agpp), color = point_color) +
            xlim(0, 3500) +
            ylim(0, 3500) +
            xlab("Predicted") +
            ylab("True") +
            scale_color_brewer(palette = "Set3") +
            geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
            ggtitle(paste(j, "(", i, ")")) +
            theme(legend.position = "none")

        st_plots <- append(st_plots, list(fig))
    }
}
pdf("EVI2_based_Bayesian_model_per_site.pdf", width = 40, height = 35)
wrap_plots(st_plots, nrow = 10)
dev.off()



# ~ Combine figures
# ~~~~~~~~~~~~~~~~~
g_fit_fig + e_fit_fig
g_beta / e_beta
g_eta / e_eta


#endregion [Overall]




#region [Cross validation leave one year out] ####
#************************************************************
cv_figs <- list()
for (i in 2001:2014) {
    # left one year out
    ind <- grep(i, north_sites_dt[IGBP != "CSH", ]$year)

    # ~ Train data
    # ~~~~~~~~~~~~~~~~~
    Y_train <- north_sites_dt[IGBP != "CSH", ]$annual_gpp[-ind]
    site_train <- as.numeric(as.factor(north_sites_dt[IGBP != "CSH"]$site))[-ind]
    igbp_train <- as.numeric(as.factor(north_sites_dt[IGBP != "CSH"]$IGBP))[-ind]

    n_train <- length(Y_train) # number of samples

    # EVI2-based
    emin_train <- north_sites_dt[IGBP != "CSH"]$m_EVImin_1[-ind]
    emax_train <- north_sites_dt[IGBP != "CSH"]$m_EVImax_1[-ind]
    gsl_train <- north_sites_dt[IGBP != "CSH"]$m_gsl[-ind]


    # ~ Test data
    # ~~~~~~~~~~~~~~~~~
    Y_test <- north_sites_dt[IGBP != "CSH", ]$annual_gpp[ind]
    site_test <- as.numeric(as.factor(north_sites_dt[IGBP != "CSH"]$site))[ind]
    igbp_test <- as.numeric(as.factor(north_sites_dt[IGBP != "CSH"]$IGBP))[ind]

    n_test <- length(Y_test) # number of samples

    # EVI2-based
    emin_test <- north_sites_dt[IGBP != "CSH"]$m_EVImin_1[ind]
    emax_test <- north_sites_dt[IGBP != "CSH"]$m_EVImax_1[ind]
    gsl_test <- north_sites_dt[IGBP != "CSH"]$m_gsl[ind]


    # ~ Fit Bayesian heriarchical model ####
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    model_string <- "model {
        # likelihood
        for (i in 1:n) {
            Y[i] ~ dnorm(f[i], tau_y[igbp[i]])
            f[i] <- beta0[site[i]] + beta1[site[i]] * emax[i] + beta2[site[i]] * emin[i] + beta3[site[i]] * gsl[i]
        }
    
        # priors
        for (i in 1:ns) {
            # beta0[i] ~ dnorm(eta0[site_igbp[i]], tau0[site_igbp[i]])
            beta0[i] ~ dnorm(eta0[site_igbp[i]], tau0)
            beta1[i] ~ dnorm(eta1[site_igbp[i]], tau1)
            beta2[i] ~ dnorm(eta2[site_igbp[i]], tau2)
            beta3[i] ~ dnorm(eta3[site_igbp[i]], tau3)
        }
        for (i in 1:nigbp) {
            eta0[i] ~ dnorm(mu[1], lambda[1])
            eta1[i] ~ dnorm(mu[2], lambda[2])
            eta2[i] ~ dnorm(mu[3], lambda[3])
            eta3[i] ~ dnorm(mu[4], lambda[4])
    
            # tau0[i] ~ dgamma(0.1, 0.1)
            # tau1[i] ~ dgamma(0.1, 0.1)
            # tau2[i] ~ dgamma(0.1, 0.1)
            # tau3[i] ~ dgamma(0.1, 0.1)
        }
    
        tau0 ~ dgamma(0.1, 0.1)
        tau1 ~ dgamma(0.1, 0.1)
        tau2 ~ dgamma(0.1, 0.1)
        tau3 ~ dgamma(0.1, 0.1)
    
        for (i in 1:4) { 
            mu[i] ~ dnorm(0, gamma[i])
            gamma[i] ~ dgamma(0.1, 0.1)
        }
        for (i in 1:4) { lambda[i] ~ dgamma(0.1, 0.1) }
        for (i in 1:nigbp) { tau_y[i] ~ dgamma(0.1, 0.1) }
        # tau_y ~ dgamma(0.1, 0.1)
    
        # predict
        for (i in 1:n_test) {
            test[i] <- beta0[site_test[i]] + beta1[site_test[i]] * emax_test[i] + beta2[site_test[i]] * emin_test[i] + beta3[site_test[i]] * gsl_test[i]
        }
    }"

    data <- list(
        igbp = igbp, ns = ns, nigbp = nigbp, site_igbp = site_igbp,
        Y = Y_train, site = site_train, n = n_train, emin = emin_train, emax = emax_train, gsl = gsl_train,
        site_test = site_test, n_test = n_test, emin_test = emin_test, emax_test = emax_test, gsl_test = gsl_test
    )
    inits <- list(beta0 = rep(-500, ns), beta1 = rep(1300, ns), beta2 = rep(2500, ns), beta3 = rep(5, ns))
    model <- jags.model(textConnection(model_string), data = data, inits = inits, n.chains = 2, quiet = TRUE)
    update(model, 20000, progress.bar = "none")
    test_samp <- coda.samples(model,
        variable.names = c("test"),
        n.iter = 50000,
        thin = 10,
        progress.bar = "none"
    )

    # ~ Plot test result
    # ~~~~~~~~~~~~~~~~~
    test <- summary(test_samp)$quantiles[1:n_test, 3]
    test_df <- data.frame(test = test, agpp = Y_test, igbp = north_sites_dt[IGBP != "CSH"]$IGBP[ind])

    cv_figs <- append(cv_figs, list(test_df))
}

figs <- list()
years <- 2001:2014
for (i in 1:length(cv_figs)) {
    fig <- ggplot(data.frame(test = cv_figs[[i]]$test, agpp = cv_figs[[i]]$agpp, igbp = cv_figs[[i]]$igbp)) +
        geom_point(aes(x = test, y = agpp, color = igbp)) +
        xlim(0, 3500) +
        ylim(0, 3500) +
        xlab("Predicted") +
        ylab("True") +
        scale_color_brewer(palette = "Set3") +
        geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
        ggtitle(years[i])
    figs[[i]] <- fig
}
wrap_plots(figs, nrow = 4)



#endregion [Cross validation leave one year out]




#region [Cross validation leave one site out] ####
#************************************************************
cv_res <- list()
progress <- 1

no_csh_dt <- north_sites_dt[IGBP != "CSH", ]
unique_sites <- unique(no_csh_dt[, site])
site_igbp <- as.factor(unique(no_csh_dt[, .(site, IGBP)])$IGBP)
nigbp <- length(unique(site_igbp))
for (st in unique_sites) {
    # left one site out
    # st <- unique_sites[2]

    test_ind <- grep(st, no_csh_dt$site)

    # ~ Train data
    # ~~~~~~~~~~~~~~~~~
    site_train <- as.numeric(as.factor(no_csh_dt[site != st, site]))
    igbp_train <- as.numeric(as.factor(no_csh_dt[site != st, IGBP]))
    site_igbp_train <- as.numeric(factor(unique(no_csh_dt[site != st, .(site, IGBP)])$IGBP, levels = levels(site_igbp)))

    Y_train <- no_csh_dt$annual_gpp[-test_ind]
    n_train <- length(Y_train) # number of samples
    ns_train <- length(site_igbp_train)

    # EVI2-based
    emin_train <- no_csh_dt$m_EVImin_1[-test_ind] * 100
    emax_train <- no_csh_dt$m_EVImax_1[-test_ind] * 100
    gsl_train <- no_csh_dt$m_gsl[-test_ind]

    # GPP-based
    emin_train <- no_csh_dt$f_gppmin_1[-test_ind]
    emax_train <- no_csh_dt$f_gppmax_1[-test_ind]
    gsl_train <- no_csh_dt$f_gsl[-test_ind]


    # ~ Test data
    # ~~~~~~~~~~~~~~~~~
    Y_test <- no_csh_dt$annual_gpp[test_ind]
    igbp_test <- as.numeric(as.factor(no_csh_dt$IGBP))[test_ind]

    # EVI2-based
    emin_test <- no_csh_dt$m_EVImin_1[test_ind]
    emax_test <- no_csh_dt$m_EVImax_1[test_ind]
    gsl_test <- no_csh_dt$m_gsl[test_ind]

    # GPP-based
    emin_test <- no_csh_dt$f_gppmin_1[test_ind]
    emax_test <- no_csh_dt$f_gppmax_1[test_ind]
    gsl_test <- no_csh_dt$f_gsl[test_ind]

    # ~ Fit Bayesian heriarchical model ####
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    model_string <- "model {
        # likelihood
        for (i in 1:n) {
            Y[i] ~ dnorm(f[i], tau_y[igbp[i]])
            f[i] <- beta0[site[i]] + beta2[site[i]] * emin[i] + beta3[site[i]] * gsl[i]
        }
    
        # priors
        for (i in 1:ns) {
            # beta0[i] ~ dnorm(eta0[site_igbp[i]], tau0[site_igbp[i]])
            beta0[i] ~ dnorm(eta0[site_igbp[i]], tau0)
            # beta1[i] ~ dnorm(eta1[site_igbp[i]], tau1)
            beta2[i] ~ dnorm(eta2[site_igbp[i]], tau2)
            beta3[i] ~ dnorm(eta3[site_igbp[i]], tau3)
        }
        for (i in 1:nigbp) {
            eta0[i] ~ dnorm(mu[1], lambda[1])
            # eta1[i] ~ dnorm(mu[2], lambda[2])
            eta2[i] ~ dnorm(mu[3], lambda[3])
            eta3[i] ~ dnorm(mu[4], lambda[4])
    
            # tau0[i] ~ dgamma(0.1, 0.1)
            # tau1[i] ~ dgamma(0.1, 0.1)
            # tau2[i] ~ dgamma(0.1, 0.1)
            # tau3[i] ~ dgamma(0.1, 0.1)
        }
    
        tau0 ~ dgamma(0.1, 0.1)
        # tau1 ~ dgamma(0.1, 0.1)
        tau2 ~ dgamma(0.1, 0.1)
        tau3 ~ dgamma(0.1, 0.1)
    
        for (i in 1:4) { 
            mu[i] ~ dnorm(0, gamma[i])
            gamma[i] ~ dgamma(0.1, 0.1)
        }
        for (i in 1:4) { lambda[i] ~ dgamma(0.1, 0.1) }
        for (i in 1:nigbp) { tau_y[i] ~ dgamma(0.1, 0.1) }
        # tau_y ~ dgamma(0.1, 0.1)
    }"

    data <- list(
        igbp = igbp_train, ns = ns_train, nigbp = nigbp, site_igbp = site_igbp_train,
        Y = Y_train, site = site_train, n = n_train, emin = emin_train, emax = emax_train, gsl = gsl_train
    )
    inits <- list(beta0 = rep(-500, ns_train), beta2 = rep(2500, ns_train), beta3 = rep(5, ns_train))
    model <- jags.model(textConnection(model_string), data = data, inits = inits, n.chains = 2, quiet = TRUE)
    update(model, 20000, progress.bar = "none")
    eta_samp <- coda.samples(model,
        variable.names = c("eta0", "eta1", "eta2", "eta3"),
        n.iter = 500000,
        # thin = 10,
        progress.bar = "none"
    )
    beta_samp <- coda.samples(model,
        variable.names = c("beta0", "beta1", "beta2", "beta3"),
        n.iter = 50000,
        # thin = 10,
        progress.bar = "none"
    )
    eta_samp <- rbind(eta_samp[[1]], eta_samp[[2]])
    beta_samp <- rbind(beta_samp[[1]], beta_samp[[2]])

eta0_i <- summary(eta_samp)[, 6]
eta1_i <- summary(eta_samp)[, 16]
eta2_i <- summary(eta_samp)[, 26]
eta3_i <- summary(eta_samp)[, 36]

eta_sum <- cbind(eta0_i, eta1_i, eta2_i, eta3_i)

this_st_igbp <- grep(igbp_test[1], site_igbp_train)
beta0_i <- summary(beta_samp)[, colnames(summary(beta_samp)) %like% "beta0"][, this_st_igbp]
beta1_i <- summary(beta_samp)[, colnames(summary(beta_samp)) %like% "beta1"][, this_st_igbp]
beta2_i <- summary(beta_samp)[, colnames(summary(beta_samp)) %like% "beta2"][, this_st_igbp]
beta3_i <- summary(beta_samp)[, colnames(summary(beta_samp)) %like% "beta3"][, this_st_igbp]
beta_sum <- cbind(beta0_i, beta1_i, beta2_i, beta3_i)

effectiveSize(eta_samp)
gelman.diag(eta_samp)

eta0_i <- eta_samp[, 6]
eta1_i <- eta_samp[, 16]
eta2_i <- eta_samp[, 26]
eta3_i <- eta_samp[, 36]

cor(cbind(eta0_i, eta1_i))

plot(emin_train)
plot(emax_train)
plot(gsl_train)


    # ~ Make prediction
    # ~~~~~~~~~~~~~~~~~
    pred <- data.frame(matrix(NA, nrow = length(Y_test), ncol = 3))
    colnames(pred) <- c("lower", "median", "upper")
    for (i in 1:length(Y_test)) {
        pred_new <- apply(eta_samp, 1, function(x, igbp_type = igbp_test[i], emin_new = emin_test[i], emax_new = emax_test[i], gsl_new = gsl_test[i]) {
            eta0 <- x[[paste0("eta0", "[", igbp_type, "]")]]
            eta1 <- x[[paste0("eta1", "[", igbp_type, "]")]]
            eta2 <- x[[paste0("eta2", "[", igbp_type, "]")]]
            eta3 <- x[[paste0("eta3", "[", igbp_type, "]")]]
            eta0 + eta1 * emax_new + eta2 * emin_new + eta3 * gsl_new
        })
        pred[i, ] <- quantile(pred_new, c(0.025, 0.6, 0.975))
    }
    pred$true <- Y_test

    cv_res <<- append(cv_res, list(pred))

    progress <- progress + 1
    print(progress)
}
names(cv_res) <- unique_sites

saveRDS(cv_res, file = "cv_res_gpp.Rds")



figs <- list()
point_color <- brewer.pal(10, "Set3")
for (i in unique(no_csh_dt[, IGBP])) {
    cur_igbp_sites <- unique(no_csh_dt[IGBP == i, site])
    igbp_test <- grep(i, levels(site_igbp))
    for (j in cur_igbp_sites) {
        fig <- ggplot(data.frame(pred = cv_res[[j]]$median, agpp = cv_res[[j]]$true, lower = cv_res[[j]]$lower, upper = cv_res[[j]]$upper)) +
            geom_point(aes(x = pred, y = agpp), color = point_color[igbp_test]) +
            geom_errorbar(aes(x = pred, y = agpp, xmin = lower, xmax = upper), color = point_color[igbp_test]) +
            xlim(0, 3500) +
            ylim(0, 3500) +
            xlab("Predicted") +
            ylab("True") +
            scale_color_brewer(palette = "Set3") +
            geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
            ggtitle(paste(j, "(", i, ")"))
        figs <- append(figs, list(fig))
    }
}
pdf("leave_one_site_out_cv_gpp.pdf", width = 40, height = 35)
wrap_plots(figs, nrow = 10)
dev.off()


#endregion [Cross validation leave one site out]




#region [Isolate midgup and midgdown effects] ####
#************************************************************
# ~ Fit Bayesian heriarchical model ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
model_string <- "model {
    # likelihood
    for (i in 1:n) {
        Y[i] ~ dnorm(f[i], tau_y[igbp[i]])
        f[i] <- beta0[site[i]] + beta1[site[i]] * emax[i] + beta2[site[i]] * emin[i] + beta3[site[i]] * midgup[i] + beta4[site[i]] * midgdown[i]
    }

    # priors
    for (i in 1:ns) {
        beta0[i] ~ dnorm(eta0[site_igbp[i]], tau0)
        beta1[i] ~ dnorm(eta1[site_igbp[i]], tau1)
        beta2[i] ~ dnorm(eta2[site_igbp[i]], tau2)
        beta3[i] ~ dnorm(eta3[site_igbp[i]], tau3)
        beta4[i] ~ dnorm(eta4[site_igbp[i]], tau4)
    }
    for (i in 1:nigbp) {
        eta0[i] ~ dnorm(mu[1], lambda[1])
        eta1[i] ~ dnorm(mu[2], lambda[2])
        eta2[i] ~ dnorm(mu[3], lambda[3])
        eta3[i] ~ dnorm(mu[4], lambda[4])
        eta4[i] ~ dnorm(mu[5], lambda[5])

        # tau0[i] ~ dgamma(0.1, 0.1)
        # tau1[i] ~ dgamma(0.1, 0.1)
        # tau2[i] ~ dgamma(0.1, 0.1)
        # tau3[i] ~ dgamma(0.1, 0.1)
    }

    tau0 ~ dgamma(0.1, 0.1)
    tau1 ~ dgamma(0.1, 0.1)
    tau2 ~ dgamma(0.1, 0.1)
    tau3 ~ dgamma(0.1, 0.1)
    tau4 ~ dgamma(0.1, 0.1)

    for (i in 1:5) { 
        mu[i] ~ dnorm(0, gamma[i])
        gamma[i] ~ dgamma(0.1, 0.1)
    }
    for (i in 1:5) { lambda[i] ~ dgamma(0.1, 0.1) }
    for (i in 1:nigbp) { tau_y[i] ~ dgamma(0.1, 0.1) }
    # tau_y ~ dgamma(0.1, 0.1)

    # predict
    for (i in 1:n) {
        fit[i] <- beta0[site[i]] + beta1[site[i]] * emax[i] + beta2[site[i]] * emin[i] + beta3[site[i]] * midgup[i] + beta4[site[i]] * midgdown[i]
    }
}"

# ~ data ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Y <- north_sites_dt[IGBP != "CSH"]$annual_gpp

site <- as.numeric(as.factor(north_sites_dt[IGBP != "CSH"]$site))
igbp <- as.numeric(as.factor(north_sites_dt[IGBP != "CSH"]$IGBP))
n <- length(Y) # number of samples
ns <- length(unique(site)) # number of sites
nigbp <- length(unique(igbp)) # number of igbp types

site_igbp <- as.numeric(as.factor(unique(north_sites_dt[IGBP != "CSH", .(site, IGBP)])$IGBP))


# ~ GPP-based ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
emin <- as.vector(scale(north_sites_dt[IGBP != "CSH"]$f_gppmin_1, scale = FALSE))
emax <- as.vector(scale(north_sites_dt[IGBP != "CSH"]$f_gppmax_1, scale = FALSE))
midgup <- as.vector(scale(north_sites_dt[IGBP != "CSH"]$f_midgup_1, scale = FALSE))
midgdown <- as.vector(scale(north_sites_dt[IGBP != "CSH"]$f_midgdown_1, scale = FALSE))

data <- list(
    Y = Y, site = site, igbp = igbp, n = n, ns = ns, nigbp = nigbp, site_igbp = site_igbp,
    emin = emin, emax = emax, midgup = midgup, midgdown = midgdown
)
inits <- list(beta0 = rep(1300, ns), beta1 = rep(100, ns), beta2 = rep(100, ns), beta3 = rep(-5, ns), beta4 = rep(5, ns))
model <- jags.model(textConnection(model_string), data = data, inits = inits, n.chains = 2, quiet = TRUE)
update(model, 20000, progress.bar = "none")
beta_samp <- coda.samples(model,
    variable.names = c("beta0", "beta1", "beta2", "beta3", "beta4"),
    n.iter = 100000,
    thin = 10,
    progress.bar = "none"
)
eta_samp <- coda.samples(model,
    variable.names = c("eta0", "eta1", "eta2", "eta3", "eta4"),
    n.iter = 100000,
    thin = 10,
    progress.bar = "none"
)
mu_samp <- coda.samples(model,
    variable.names = c("mu"),
    n.iter = 100000,
    thin = 10,
    progress.bar = "none"
)
fit_samp <- coda.samples(model,
    variable.names = c("fit"),
    n.iter = 100000,
    thin = 10,
    progress.bar = "none"
)


# ~ Check convergence
# ~~~~~~~~~~~~~~~~~
# convergence of beta
gelman.diag(beta_samp)
geweke.diag(beta_samp)
effectiveSize(beta_samp)
plot(beta_samp)

# convergence of eta
gelman.diag(eta_samp)
geweke.diag(eta_samp)
effectiveSize(eta_samp)
plot(eta_samp)

# convergence of mu
gelman.diag(mu_samp)
geweke.diag(mu_samp)
effectiveSize(mu_samp)
plot(mu_samp)

# ~ Draw figures
# ~~~~~~~~~~~~~~~~~
beta_quan <- list()
beta_quan$beta0 <- summary(beta_samp)$quantiles[1:ns, c(1, 3, 5)]
beta_quan$beta1 <- summary(beta_samp)$quantiles[1 * ns + (1:ns), c(1, 3, 5)]
beta_quan$beta2 <- summary(beta_samp)$quantiles[2 * ns + (1:ns), c(1, 3, 5)]
beta_quan$beta3 <- summary(beta_samp)$quantiles[3 * ns + (1:ns), c(1, 3, 5)]
beta_quan$beta4 <- summary(beta_samp)$quantiles[4 * ns + (1:ns), c(1, 3, 5)]

eta_quan <- list()
eta_quan$eta0 <- summary(eta_samp)$quantiles[1:nigbp, c(1, 3, 5)]
eta_quan$eta1 <- summary(eta_samp)$quantiles[1 * nigbp + (1:nigbp), c(1, 3, 5)]
eta_quan$eta2 <- summary(eta_samp)$quantiles[2 * nigbp + (1:nigbp), c(1, 3, 5)]
eta_quan$eta3 <- summary(eta_samp)$quantiles[3 * nigbp + (1:nigbp), c(1, 3, 5)]
eta_quan$eta4 <- summary(eta_samp)$quantiles[4 * nigbp + (1:nigbp), c(1, 3, 5)]

mu_quan <- summary(mu_samp)$quantiles[, c(1, 3, 5)]
rownames(mu_quan) <- c("(Intercept)", "Max", "Min", "Midgup", "Midgdown")


# fig: figures
g_beta0_df <- data.table(y = beta_quan$beta0, unique(north_sites_dt[IGBP != "CSH", .(site, IGBP)]))
colnames(g_beta0_df) <- c("gsl_min", "gsl_median", "gsl_max", "Site", "IGBP")
setorder(g_beta0_df, cols = IGBP)
setkey(g_beta0_df, "IGBP")
g_beta0 <- ggplot(g_beta0_df) +
    geom_pointrange(aes(x = 1:ns, y = gsl_median, ymax = gsl_max, ymin = gsl_min, color = IGBP)) +
    xlab("Site") +
    ylab("Intercept") +
    ggtitle("GPP-based model")

g_beta1_df <- data.table(y = beta_quan$beta1, unique(north_sites_dt[IGBP != "CSH", .(site, IGBP)]))
colnames(g_beta1_df) <- c("gsl_min", "gsl_median", "gsl_max", "Site", "IGBP")
setorder(g_beta1_df, cols = IGBP)
setkey(g_beta1_df, "IGBP")
g_beta1 <- ggplot(g_beta1_df) +
    geom_pointrange(aes(x = 1:ns, y = gsl_median, ymax = gsl_max, ymin = gsl_min, color = IGBP)) +
    xlab("Site") +
    ylab("GPPmax coefficient") +
    ggtitle("GPP-based model")

g_beta2_df <- data.table(y = beta_quan$beta2, unique(north_sites_dt[IGBP != "CSH", .(site, IGBP)]))
colnames(g_beta2_df) <- c("gsl_min", "gsl_median", "gsl_max", "Site", "IGBP")
setorder(g_beta2_df, cols = IGBP)
setkey(g_beta2_df, "IGBP")
g_beta2 <- ggplot(g_beta2_df) +
    geom_pointrange(aes(x = 1:ns, y = gsl_median, ymax = gsl_max, ymin = gsl_min, color = IGBP)) +
    xlab("Site") +
    ylab("GPPmin coefficient") +
    ggtitle("GPP-based model")

g_beta3_df <- data.table(y = beta_quan$beta3, unique(north_sites_dt[IGBP != "CSH", .(site, IGBP)]))
colnames(g_beta3_df) <- c("gsl_min", "gsl_median", "gsl_max", "Site", "IGBP")
setorder(g_beta3_df, cols = IGBP)
setkey(g_beta3_df, "IGBP")
g_beta3 <- ggplot(g_beta3_df) +
    geom_pointrange(aes(x = 1:ns, y = gsl_median, ymax = gsl_max, ymin = gsl_min, color = IGBP)) +
    xlab("Site") +
    ylab("Midgup coefficient") +
    ggtitle("GPP-based model")

g_beta4_df <- data.table(y = beta_quan$beta4, unique(north_sites_dt[IGBP != "CSH", .(site, IGBP)]))
colnames(g_beta4_df) <- c("gsl_min", "gsl_median", "gsl_max", "Site", "IGBP")
setorder(g_beta4_df, cols = IGBP)
setkey(g_beta4_df, "IGBP")
g_beta4 <- ggplot(g_beta4_df) +
    geom_pointrange(aes(x = 1:ns, y = gsl_median, ymax = gsl_max, ymin = gsl_min, color = IGBP)) +
    xlab("Site") +
    ylab("Midgdown coefficient") +
    ggtitle("GPP-based model")

g_beta0 / g_beta1 / g_beta2 / g_beta3 / g_beta4


# fig: eta figures
g_eta0 <- ggplot(data.frame(x = levels(as.factor(north_sites_dt[IGBP != "CSH"]$IGBP)), y = eta_quan$eta0)) +
    geom_pointrange(aes(x = x, y = y.50., ymax = y.97.5., ymin = y.2.5.)) +
    xlab("IGBP") +
    ylab("Intercept") +
    ggtitle("GPP-based model")
g_eta1 <- ggplot(data.frame(x = levels(as.factor(north_sites_dt[IGBP != "CSH"]$IGBP)), y = eta_quan$eta1)) +
    geom_pointrange(aes(x = x, y = y.50., ymax = y.97.5., ymin = y.2.5.)) +
    xlab("IGBP") +
    ylab("GPPmax coefficient") +
    ggtitle("GPP-based model")
g_eta2 <- ggplot(data.frame(x = levels(as.factor(north_sites_dt[IGBP != "CSH"]$IGBP)), y = eta_quan$eta2)) +
    geom_pointrange(aes(x = x, y = y.50., ymax = y.97.5., ymin = y.2.5.)) +
    xlab("IGBP") +
    ylab("GPPmin coefficient") +
    ggtitle("GPP-based model")
g_eta3 <- ggplot(data.frame(x = levels(as.factor(north_sites_dt[IGBP != "CSH"]$IGBP)), y = eta_quan$eta3)) +
    geom_pointrange(aes(x = x, y = y.50., ymax = y.97.5., ymin = y.2.5.)) +
    xlab("IGBP") +
    ylab("Midgup coefficient") +
    ggtitle("GPP-based model")
g_eta4 <- ggplot(data.frame(x = levels(as.factor(north_sites_dt[IGBP != "CSH"]$IGBP)), y = eta_quan$eta4)) +
    geom_pointrange(aes(x = x, y = y.50., ymax = y.97.5., ymin = y.2.5.)) +
    xlab("IGBP") +
    ylab("Midgdown coefficient") +
    ggtitle("GPP-based model")

g_eta0 / g_eta1 / g_eta2 / g_eta3 / g_eta4

# ~ Plot fitted
# ~~~~~~~~~~~~~~~~~
fit <- summary(fit_samp)$quantiles[1:n, 3]
ggplot(data.frame(fit = fit, agpp = Y, igbp = north_sites_dt[IGBP != "CSH"]$IGBP)) +
    geom_point(aes(x = fit, y = agpp, color = igbp)) +
    xlim(0, 3500) +
    ylim(0, 3500) +
    xlab("Fitted") +
    ylab("True") +
    scale_color_brewer(palette = "Set3") +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    ggtitle("EVI2-based model")


plot(1)





# ~ EVI2 based ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
emin <- as.vector(scale(north_sites_dt[IGBP != "CSH"]$m_EVImin_1, scale = FALSE)) * 100
emax <- as.vector(scale(north_sites_dt[IGBP != "CSH"]$m_EVImax_1, scale = FALSE)) * 100
midgup <- as.vector(scale(north_sites_dt[IGBP != "CSH"]$m_midgup_1, scale = FALSE))
midgdown <- as.vector(scale(north_sites_dt[IGBP != "CSH"]$m_midgdown_1, scale = FALSE))

data <- list(
    Y = Y, site = site, igbp = igbp, n = n, ns = ns, nigbp = nigbp, site_igbp = site_igbp,
    emin = emin, emax = emax, midgup = midgup, midgdown = midgdown
)
inits <- list(beta0 = rep(-500, ns), beta1 = rep(1300, ns), beta2 = rep(2500, ns), beta3 = rep(5, ns), beta4 = rep(5, ns))
model <- jags.model(textConnection(model_string), data = data, inits = inits, n.chains = 2, quiet = TRUE)
update(model, 20000, progress.bar = "none")
beta_samp <- coda.samples(model,
    variable.names = c("beta0", "beta1", "beta2", "beta3", "beta4"),
    n.iter = 100000,
    thin = 10,
    progress.bar = "none"
)
eta_samp <- coda.samples(model,
    variable.names = c("eta0", "eta1", "eta2", "eta3", "eta4"),
    n.iter = 100000,
    thin = 10,
    progress.bar = "none"
)
mu_samp <- coda.samples(model,
    variable.names = c("mu"),
    n.iter = 100000,
    thin = 10,
    progress.bar = "none"
)
fit_samp <- coda.samples(model,
    variable.names = c("fit"),
    n.iter = 100000,
    thin = 10,
    progress.bar = "none"
)



# ~ Check convergence
# ~~~~~~~~~~~~~~~~~
# convergence of beta
gelman.diag(beta_samp)
geweke.diag(beta_samp)
effectiveSize(beta_samp)
plot(beta_samp)

# convergence of eta
gelman.diag(eta_samp)
geweke.diag(eta_samp)
effectiveSize(eta_samp)
plot(eta_samp)

# convergence of mu
gelman.diag(mu_samp)
geweke.diag(mu_samp)
effectiveSize(mu_samp)
plot(mu_samp)


# ~ Draw figures
# ~~~~~~~~~~~~~~~~~
beta_quan <- list()
beta_quan$beta0 <- summary(beta_samp)$quantiles[1:ns, c(1, 3, 5)]
beta_quan$beta1 <- summary(beta_samp)$quantiles[1 * ns + (1:ns), c(1, 3, 5)]
beta_quan$beta2 <- summary(beta_samp)$quantiles[2 * ns + (1:ns), c(1, 3, 5)]
beta_quan$beta3 <- summary(beta_samp)$quantiles[3 * ns + (1:ns), c(1, 3, 5)]
beta_quan$beta4 <- summary(beta_samp)$quantiles[4 * ns + (1:ns), c(1, 3, 5)]

eta_quan <- list()
eta_quan$eta0 <- summary(eta_samp)$quantiles[1:nigbp, c(1, 3, 5)]
eta_quan$eta1 <- summary(eta_samp)$quantiles[1 * nigbp + (1:nigbp), c(1, 3, 5)]
eta_quan$eta2 <- summary(eta_samp)$quantiles[2 * nigbp + (1:nigbp), c(1, 3, 5)]
eta_quan$eta3 <- summary(eta_samp)$quantiles[3 * nigbp + (1:nigbp), c(1, 3, 5)]
eta_quan$eta4 <- summary(eta_samp)$quantiles[4 * nigbp + (1:nigbp), c(1, 3, 5)]

mu_quan <- summary(mu_samp)$quantiles[, c(1, 3, 5)]
rownames(mu_quan) <- c("(Intercept)", "Max", "Min", "Midgup", "Midgdown")


# fig: figures
e_beta0_df <- data.table(y = beta_quan$beta0, unique(north_sites_dt[IGBP != "CSH", .(site, IGBP)]))
colnames(e_beta0_df) <- c("gsl_min", "gsl_median", "gsl_max", "Site", "IGBP")
setorder(e_beta0_df, cols = IGBP)
setkey(e_beta0_df, "IGBP")
e_beta0 <- ggplot(e_beta0_df) +
    geom_pointrange(aes(x = 1:ns, y = gsl_median, ymax = gsl_max, ymin = gsl_min, color = IGBP)) +
    xlab("Site") +
    ylab("Intercept") +
    ggtitle("EVI2-based model")

e_beta1_df <- data.table(y = beta_quan$beta1, unique(north_sites_dt[IGBP != "CSH", .(site, IGBP)]))
colnames(e_beta1_df) <- c("gsl_min", "gsl_median", "gsl_max", "Site", "IGBP")
setorder(e_beta1_df, cols = IGBP)
setkey(e_beta1_df, "IGBP")
e_beta1 <- ggplot(e_beta1_df) +
    geom_pointrange(aes(x = 1:ns, y = gsl_median, ymax = gsl_max, ymin = gsl_min, color = IGBP)) +
    xlab("Site") +
    ylab("EVI2max coefficient") +
    ggtitle("EVI2-based model")

e_beta2_df <- data.table(y = beta_quan$beta2, unique(north_sites_dt[IGBP != "CSH", .(site, IGBP)]))
colnames(e_beta2_df) <- c("gsl_min", "gsl_median", "gsl_max", "Site", "IGBP")
setorder(e_beta2_df, cols = IGBP)
setkey(e_beta2_df, "IGBP")
e_beta2 <- ggplot(e_beta2_df) +
    geom_pointrange(aes(x = 1:ns, y = gsl_median, ymax = gsl_max, ymin = gsl_min, color = IGBP)) +
    xlab("Site") +
    ylab("EVI2min coefficient") +
    ggtitle("EVI2-based model")

e_beta3_df <- data.table(y = beta_quan$beta3, unique(north_sites_dt[IGBP != "CSH", .(site, IGBP)]))
colnames(e_beta3_df) <- c("gsl_min", "gsl_median", "gsl_max", "Site", "IGBP")
setorder(e_beta3_df, cols = IGBP)
setkey(e_beta3_df, "IGBP")
e_beta3 <- ggplot(e_beta3_df) +
    geom_pointrange(aes(x = 1:ns, y = gsl_median, ymax = gsl_max, ymin = gsl_min, color = IGBP)) +
    xlab("Site") +
    ylab("Midgup coefficient") +
    ggtitle("EVI2-based model")

e_beta4_df <- data.table(y = beta_quan$beta4, unique(north_sites_dt[IGBP != "CSH", .(site, IGBP)]))
colnames(e_beta4_df) <- c("gsl_min", "gsl_median", "gsl_max", "Site", "IGBP")
setorder(e_beta4_df, cols = IGBP)
setkey(e_beta4_df, "IGBP")
e_beta4 <- ggplot(e_beta4_df) +
    geom_pointrange(aes(x = 1:ns, y = gsl_median, ymax = gsl_max, ymin = gsl_min, color = IGBP)) +
    xlab("Site") +
    ylab("Midgdown coefficient") +
    ggtitle("EVI2-based model")

e_beta0 / e_beta1 / e_beta2 / e_beta3 / e_beta4


# fig: eta figures
e_eta0 <- ggplot(data.frame(x = levels(as.factor(north_sites_dt[IGBP != "CSH"]$IGBP)), y = eta_quan$eta0)) +
    geom_pointrange(aes(x = x, y = y.50., ymax = y.97.5., ymin = y.2.5.)) +
    xlab("IGBP") +
    ylab("Intercept") +
    ggtitle("EVI2-based model")
e_eta1 <- ggplot(data.frame(x = levels(as.factor(north_sites_dt[IGBP != "CSH"]$IGBP)), y = eta_quan$eta1)) +
    geom_pointrange(aes(x = x, y = y.50., ymax = y.97.5., ymin = y.2.5.)) +
    xlab("IGBP") +
    ylab("EVI2max coefficient") +
    ggtitle("EVI2-based model")
e_eta2 <- ggplot(data.frame(x = levels(as.factor(north_sites_dt[IGBP != "CSH"]$IGBP)), y = eta_quan$eta2)) +
    geom_pointrange(aes(x = x, y = y.50., ymax = y.97.5., ymin = y.2.5.)) +
    xlab("IGBP") +
    ylab("EVI2min coefficient") +
    ggtitle("EVI2-based model")
e_eta3 <- ggplot(data.frame(x = levels(as.factor(north_sites_dt[IGBP != "CSH"]$IGBP)), y = eta_quan$eta3)) +
    geom_pointrange(aes(x = x, y = y.50., ymax = y.97.5., ymin = y.2.5.)) +
    xlab("IGBP") +
    ylab("Midgup coefficient") +
    ggtitle("EVI2-based model")
e_eta4 <- ggplot(data.frame(x = levels(as.factor(north_sites_dt[IGBP != "CSH"]$IGBP)), y = eta_quan$eta4)) +
    geom_pointrange(aes(x = x, y = y.50., ymax = y.97.5., ymin = y.2.5.)) +
    xlab("IGBP") +
    ylab("Midgdown coefficient") +
    ggtitle("EVI2-based model")

e_eta0 / e_eta1 / e_eta2 / e_eta3 / e_eta4


# ~ Plot fitted
# ~~~~~~~~~~~~~~~~~
fit <- summary(fit_samp)$quantiles[1:n, 3]
ggplot(data.frame(fit = fit, agpp = Y, igbp = north_sites_dt[IGBP != "CSH"]$IGBP)) +
    geom_point(aes(x = fit, y = agpp, color = igbp)) +
    xlim(0, 3500) +
    ylim(0, 3500) +
    xlab("Fitted") +
    ylab("True") +
    scale_color_brewer(palette = "Set3") +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    ggtitle("EVI2-based model")




#endregion [Isolate midgup and midgdown effects]



#region [Only GSL Bayesian] ####
#************************************************************
onlyGSLBayesian <- function(data, inits, site_igbp_names, igbp_levels) {
    model_string <- "model {
        # likelihood
        for (i in 1:n) {
            Y[i] ~ dnorm(f[i], tau_y)
            f[i] <- eta1[igbp[i]] + eta2[igbp[i]] * gsl[i]
        }

        # priors
        for (i in 1:nigbp) {
            eta1[i] ~ dnorm(mu[1], lambda[1])
            eta2[i] ~ dnorm(mu[2], lambda[2])
        }

        tau_y ~ dgamma(0.1, 0.1)

        for (i in 1:2) {
            mu[i] ~ dnorm(0, 0.0001)
            lambda[i] ~ dgamma(0.1, 0.1)
        }

        sigma_y <- 1 / sqrt(tau_y)
        sigma_beta <- 1 / sqrt(lambda)
    }"

    model <- jags.model(textConnection(model_string), data = data, inits = inits, n.chains = 2, quiet = TRUE)
    update(model, 20000, progress.bar = "none")
    iter <- 1
    repeat({
        samp <- coda.samples(model,
            variable.names = c(
                "eta1", "eta2",
                "mu"
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

    samp_df <- as.data.frame(rbind(samp[[1]], samp[[2]]))
    samp_names <- colnames(samp_df)

    eta1 <- samp_df[, grep("^eta1", samp_names)] * 100
    eta2 <- samp_df[, grep("^eta2", samp_names)] * 100

    mu <- samp_df[, grep("^mu", samp_names)] * 100

    eta1_dt <- SummarySample(eta1)
    eta2_dt <- SummarySample(eta2)

    eta1_dt$IGBP <- igbp_levels
    eta2_dt$IGBP <- igbp_levels

    mu_dt <- SummarySample(mu)

    fit_samp <- coda.samples(model,
        variable.names = c("f"),
        n.iter = 50000,
        thin = 10,
        progress.bar = "none"
    )
    fit <- summary(fit_samp)$quantiles[, c(1, 3, 5)]
    colnames(fit) <- c("fit_lwr", "fit_median", "fit_upr")

    return(list(
        eta1 = eta1_dt, eta2 = eta2_dt,
        mu = mu_dt, fit = fit
    ))
}

# ~ models only consider growing season length
# ~~~~~~~~~~~~~~~~~
f_model_gsl_only <- onlyGSLBayesian(
    data = list(
        Y = Y, n = n, nigbp = f_nigbp, igbp = f_igbp,
        gsl = f_gsl
    ), inits = list(
        eta1 = runif(f_nigbp, 10, 30), eta2 = runif(f_nigbp, 0, 20)
    ), site_igbp_names = unique(north_sites_dt[, .(site, IGBP)]),
    igbp_levels = levels(f_igbp_factor)
)

m_model_gsl_only <- onlyGSLBayesian(
    data = list(
        Y = Y, n = n, nigbp = f_nigbp, igbp = f_igbp,
        gsl = m_gsl
    ), inits = list(
        eta1 = runif(f_nigbp, 10, 30), eta2 = runif(f_nigbp, 0, 20)
    ), site_igbp_names = unique(north_sites_dt[, .(site, IGBP)]),
    igbp_levels = levels(f_igbp_factor)
)

f_model_gsl_only$mu
m_model_gsl_only$mu

load("model_3_fit.RData")
f_model_gsl_only$mu$med[2] - mean(unlist(f_model_3$mu_beta[2])) * 100
m_model_gsl_only$mu$med[2] - mean(unlist(m_model_3$mu_beta[2])) * 100

# out: GSL only model
save(f_model_gsl_only, m_model_gsl_only, file = "gsl_only_model.RData")

lm_mod <- lm(Y ~ f_igbp_factor * f_gsl)
summary(lm_mod)
confint(lm_mod)
mean((coef(lm_mod)[grep("f_gsl", names(coef(lm_mod)))][1] + coef(lm_mod)[grep("f_gsl", names(coef(lm_mod)))]) * 100, na.rm = TRUE)

#endregion [Only GSL Bayesian]




