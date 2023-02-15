#************************************************************************
# Description: I'm going to do cross-validation for all the fitted models 
#   so to understand how their test error changes. The models in this file are
#   selected models for the manuscript.
# Author: Xiaojie(J) Gao
# Date: 2022-01-06
#************************************************************************

source("Code/base.R")
source("Code/mod_base.R")


# ~ Leave one year out ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
uni_yrs <- sort(unique(north_sites_dt$year))

# ~ Model 1
# ~~~~~~~~~~~~~~~~~
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

f_cv_mse_1 <- rep(0, length(uni_yrs))
m_cv_mse_1 <- rep(0, length(uni_yrs))

for (i in 1:length(uni_yrs)) {
    yr <- uni_yrs[i]

    # ~ GPP based
    f_model_train_1 <- FitBayesian(
        model_string = model_1,
        params = c("beta1", "beta2", "beta3", "beta4"),
        data = list(
            Y = north_sites_dt[year != yr]$annual_gpp / 100, n = nrow(north_sites_dt[year != yr]),
            gsl = as.vector(scale(north_sites_dt[year != yr]$f_gsl, scale = FALSE)),
            emax = as.vector(scale(north_sites_dt[year != yr]$f_gppmax_1, scale = FALSE)),
            emin = as.vector(scale(north_sites_dt[year != yr]$f_gppmin_1, scale = FALSE))
        ), inits = list(
            beta1 = runif(1, 10, 30), beta2 = runif(1, 0, 20),
            beta3 = runif(1, 10, 30), beta4 = runif(1, 0, 30)
        )
    )
    mu_f_gsl <- attr(scale(north_sites_dt[year != yr]$f_gsl, scale = FALSE), "scaled:center")
    mu_f_max <- attr(scale(north_sites_dt[year != yr]$f_gppmax_1, scale = FALSE), "scaled:center")
    mu_f_min <- attr(scale(north_sites_dt[year != yr]$f_gppmin_1, scale = FALSE), "scaled:center")

    f_beta2 <- median(f_model_train_1$beta2)
    f_beta3 <- median(f_model_train_1$beta3)
    f_beta4 <- median(f_model_train_1$beta4)
    f_beta1 <- median(f_model_train_1$beta1) - (f_beta2 * mu_f_gsl + f_beta3 * mu_f_max + f_beta4 * mu_f_min) # scale intercept to the original scale

    f_y_pred <- f_beta1 + f_beta2 * north_sites_dt[year == yr, f_gsl] +
        f_beta3 * north_sites_dt[year == yr, f_gppmax_1] + f_beta4 * north_sites_dt[year == yr, f_gppmin_1]
    # MSE
    f_cv_mse_1[i] <- mean((f_y_pred - north_sites_dt[year == yr, annual_gpp])^2)


    # ~ EVI2 based
    m_model_train_1 <- FitBayesian(
        model_string = model_1,
        params = c("beta1", "beta2", "beta3", "beta4"),
        data = list(
            Y = north_sites_dt[year != yr]$annual_gpp / 100, n = nrow(north_sites_dt[year != yr]),
            gsl = as.vector(scale(north_sites_dt[year != yr]$m_gsl, scale = FALSE)),
            emax = as.vector(scale(north_sites_dt[year != yr]$m_EVImax_1, scale = FALSE)) * 10,
            emin = as.vector(scale(north_sites_dt[year != yr]$m_EVImin_1, scale = FALSE)) * 10
        ), inits = list(
            beta1 = runif(1, 10, 30), beta2 = runif(1, 0, 20),
            beta3 = runif(1, 10, 30), beta4 = runif(1, 0, 30)
        )
    )
    mu_m_gsl <- attr(scale(north_sites_dt[year != yr]$m_gsl, scale = FALSE), "scaled:center")
    mu_m_max <- attr(scale(north_sites_dt[year != yr]$m_EVImax_1, scale = FALSE), "scaled:center")
    mu_m_min <- attr(scale(north_sites_dt[year != yr]$m_EVImin_1, scale = FALSE), "scaled:center")

    # lm(annual_gpp ~ m_gsl + m_EVImax_1 + m_EVImin_1, data = north_sites_dt[year != yr])
    m_beta2 <- median(m_model_train_1$beta2)
    m_beta3 <- median(m_model_train_1$beta3) * 10
    m_beta4 <- median(m_model_train_1$beta4) * 10
    m_beta1 <- median(m_model_train_1$beta1) - (m_beta2 * mu_m_gsl + m_beta3 * mu_m_max + m_beta4 * mu_m_min) # scale intercept to the original scale

    m_y_pred <- m_beta1 + m_beta2 * north_sites_dt[year == yr, m_gsl] +
        m_beta3 * north_sites_dt[year == yr, m_EVImax_1] + m_beta4 * north_sites_dt[year == yr, m_EVImin_1]
    # MSE
    m_cv_mse_1[i] <- mean((m_y_pred - north_sites_dt[year == yr, annual_gpp])^2)
}


plot(uni_yrs, sqrt(f_cv_mse_1), type = "b", ylim = c(100, 500))
lines(uni_yrs, sqrt(m_cv_mse_1), type = "b", col = "blue")



# ~ Model 2
# ~~~~~~~~~~~~~~~~~
model_2 <- "model {
    # likelihood
    for (i in 1:n) {
        Y[i] ~ dnorm(f[i], tau_y)
        f[i] <- beta1[igbp[i]] + beta2 * gsl[i] + beta3 * emax[i] + beta4 * emin[i]
    }

    # priors
    for (i in 1:nigbp) {
        beta1[i] ~ dnorm(mu_beta[1], tau_beta[1])
    }
    beta2 ~ dnorm(0, 0.001) # GSL effect
    beta3 ~ dnorm(0, 0.001) # Zmax effect
    beta4 ~ dnorm(0, 0.001) # Zmin effect

    mu_beta[1] ~ dunif(0, 30)
    for (i in 1:1) { tau_beta[1] ~ dgamma(0.1, 0.1) }
    tau_y ~ dgamma(0.1, 0.1)

    # prediction
    for (i in 1:ntest) {
        y_pred[i] <- beta1[igbp_test[i]] + beta2 * gsl_test[i] + beta3 * emax_test[i] + beta4 * emin_test[i]
    }
}"

f_cv_mse_2 <- rep(0, length(uni_yrs))
m_cv_mse_2 <- rep(0, length(uni_yrs))

for (i in 1:length(uni_yrs)) {
    yr <- uni_yrs[i]

    y_train <- north_sites_dt[year != yr]$annual_gpp / 100
    n_train <- nrow(north_sites_dt[year != yr])
    igbp_train <- f_igbp[north_sites_dt$year != yr]
    nigbp <- length(unique(f_igbp))

    y_test <- north_sites_dt[year == yr, annual_gpp / 100]
    n_test <- length(y_test)
    igbp_test <- f_igbp[north_sites_dt$year == yr]


    # ~ GPP based
    f_gsl_train <- as.vector(scale(north_sites_dt[year != yr]$f_gsl, scale = FALSE))
    f_emax_train <- as.vector(scale(north_sites_dt[year != yr]$f_gppmax_1, scale = FALSE))
    f_emin_train <- as.vector(scale(north_sites_dt[year != yr]$f_gppmin_1, scale = FALSE))

    mu_f_gsl <- attr(scale(north_sites_dt[year != yr]$f_gsl, scale = FALSE), "scaled:center")
    mu_f_max <- attr(scale(north_sites_dt[year != yr]$f_gppmax_1, scale = FALSE), "scaled:center")
    mu_f_min <- attr(scale(north_sites_dt[year != yr]$f_gppmin_1, scale = FALSE), "scaled:center")

    f_gsl_test <- north_sites_dt[year == yr, f_gsl] - mu_f_gsl
    f_max_test <- north_sites_dt[year == yr, f_gppmax_1] - mu_f_max
    f_min_test <- north_sites_dt[year == yr, f_gppmin_1] - mu_f_min


    f_model <- jags.model(textConnection(model_2),
        data = list(
            Y = y_train, n = n_train,
            gsl = f_gsl_train, emax = f_emax_train, emin = f_emin_train, igbp = igbp_train, nigbp = nigbp,
            ntest = n_test, gsl_test = f_gsl_test, emax_test = f_max_test, emin_test = f_min_test, igbp_test = igbp_test
        ),
        inits = list(
            beta1 = runif(nigbp, 10, 30), beta2 = runif(1, 0, 20),
            beta3 = runif(1, 10, 30), beta4 = runif(1, 0, 30)
        ),
        n.chains = 2, quiet = TRUE
    )
    update(f_model, 2000, progress.bar = "none")
    repeat({
        samp <- coda.samples(f_model,
            variable.names = c("beta1", "beta2", "beta3", "beta4", "tau_y", "tau_beta", "mu_beta"),
            n.iter = 10000,
            thin = 10,
            progress.bar = "none"
        )
        mpsrf <- gelman.diag(samp)$mpsrf

        if (abs(mpsrf) < 1.1 & max(abs(gelman.diag(samp)$psrf)) < 1.1) break
    })
    y_pred_samp <- coda.samples(f_model,
        variable.names = c("y_pred"),
        n.iter = 5000,
        thin = 10,
        progress.bar = "none"
    )
    f_pred <- summary(y_pred_samp)$quantiles[, 3]
    f_cv_mse_2[i] <- mean(((f_pred - y_test) * 100)^2)


    # ~ EVI2 based
    m_gsl_train <- as.vector(scale(north_sites_dt[year != yr]$m_gsl, scale = FALSE))
    m_emax_train <- as.vector(scale(north_sites_dt[year != yr]$m_EVImax_1, scale = FALSE)) * 10
    m_emin_train <- as.vector(scale(north_sites_dt[year != yr]$m_EVImin_1, scale = FALSE)) * 10

    mu_m_gsl <- attr(scale(north_sites_dt[year != yr]$m_gsl, scale = FALSE), "scaled:center")
    mu_m_max <- attr(scale(north_sites_dt[year != yr]$m_EVImax_1, scale = FALSE), "scaled:center")
    mu_m_min <- attr(scale(north_sites_dt[year != yr]$m_EVImin_1, scale = FALSE), "scaled:center")

    m_gsl_test <- north_sites_dt[year == yr, m_gsl] - mu_m_gsl
    m_max_test <- (north_sites_dt[year == yr, m_EVImax_1] - mu_m_max) * 10
    m_min_test <- (north_sites_dt[year == yr, m_EVImin_1] - mu_m_min) * 10


    m_model <- jags.model(textConnection(model_2),
        data = list(
            Y = y_train, n = n_train,
            gsl = m_gsl_train, emax = m_emax_train, emin = m_emin_train, igbp = igbp_train, nigbp = nigbp,
            ntest = n_test, gsl_test = m_gsl_test, emax_test = m_max_test, emin_test = m_min_test, igbp_test = igbp_test
        ),
        inits = list(
            beta1 = runif(nigbp, 10, 30), beta2 = runif(1, 0, 20),
            beta3 = runif(1, 10, 30), beta4 = runif(1, 0, 30)
        ),
        n.chains = 2, quiet = TRUE
    )
    update(m_model, 2000, progress.bar = "none")
    repeat({
        samp <- coda.samples(m_model,
            variable.names = c("beta1", "beta2", "beta3", "beta4", "tau_y", "tau_beta", "mu_beta"),
            n.iter = 10000,
            thin = 10,
            progress.bar = "none"
        )
        mpsrf <- gelman.diag(samp)$mpsrf

        if (abs(mpsrf) < 1.1 & max(abs(gelman.diag(samp)$psrf)) < 1.1) break
    })
    y_pred_samp <- coda.samples(m_model,
        variable.names = c("y_pred"),
        n.iter = 5000,
        thin = 10,
        progress.bar = "none"
    )
    m_pred <- summary(y_pred_samp)$quantiles[, 3]
    m_cv_mse_2[i] <- mean(((m_pred - y_test) * 100)^2)
}

lines(uni_yrs, sqrt(f_cv_mse_2), type = "b", ylim = c(100, 500))
lines(uni_yrs, sqrt(m_cv_mse_2), type = "b", col = "blue")


# ~ Model 3
# ~~~~~~~~~~~~~~~~~
model_3 <- "model {
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

    # prediction
    for (i in 1:ntest) {
        y_pred[i] <- beta1[igbp_test[i]] + beta2[igbp_test[i]] * gsl_test[i] + beta3[igbp_test[i]] * emax_test[i] + beta4[igbp_test[i]] * emin_test[i]
    }
}"

f_cv_mse_3 <- rep(0, length(uni_yrs))
m_cv_mse_3 <- rep(0, length(uni_yrs))

for (i in 1:length(uni_yrs)) {
    yr <- uni_yrs[i]

    y_train <- north_sites_dt[year != yr]$annual_gpp / 100
    n_train <- nrow(north_sites_dt[year != yr])
    igbp_train <- f_igbp[north_sites_dt$year != yr]
    nigbp <- length(unique(f_igbp))

    y_test <- north_sites_dt[year == yr, annual_gpp / 100]
    n_test <- length(y_test)
    igbp_test <- f_igbp[north_sites_dt$year == yr]


    # ~ GPP based
    f_gsl_train <- as.vector(scale(north_sites_dt[year != yr]$f_gsl, scale = FALSE))
    f_emax_train <- as.vector(scale(north_sites_dt[year != yr]$f_gppmax_1, scale = FALSE))
    f_emin_train <- as.vector(scale(north_sites_dt[year != yr]$f_gppmin_1, scale = FALSE))

    mu_f_gsl <- attr(scale(north_sites_dt[year != yr]$f_gsl, scale = FALSE), "scaled:center")
    mu_f_max <- attr(scale(north_sites_dt[year != yr]$f_gppmax_1, scale = FALSE), "scaled:center")
    mu_f_min <- attr(scale(north_sites_dt[year != yr]$f_gppmin_1, scale = FALSE), "scaled:center")

    f_gsl_test <- north_sites_dt[year == yr, f_gsl] - mu_f_gsl
    f_max_test <- north_sites_dt[year == yr, f_gppmax_1] - mu_f_max
    f_min_test <- north_sites_dt[year == yr, f_gppmin_1] - mu_f_min


    f_model <- jags.model(textConnection(model_3),
        data = list(
            Y = y_train, n = n_train,
            gsl = f_gsl_train, emax = f_emax_train, emin = f_emin_train, igbp = igbp_train, nigbp = nigbp,
            ntest = n_test, gsl_test = f_gsl_test, emax_test = f_max_test, emin_test = f_min_test, igbp_test = igbp_test
        ),
        inits = list(
            beta1 = runif(nigbp, 10, 30), beta2 = runif(nigbp, 0, 20),
            beta3 = runif(nigbp, 10, 30), beta4 = runif(nigbp, 0, 30)
        ),
        n.chains = 2, quiet = TRUE
    )
    update(f_model, 2000, progress.bar = "none")
    repeat({
        samp <- coda.samples(f_model,
            variable.names = c("beta1", "beta2", "beta3", "beta4", "tau_y", "tau_beta", "mu_beta"),
            n.iter = 10000,
            thin = 10,
            progress.bar = "none"
        )
        mpsrf <- gelman.diag(samp)$mpsrf

        if (abs(mpsrf) < 1.1 & max(abs(gelman.diag(samp)$psrf)) < 1.1) break
    })
    y_pred_samp <- coda.samples(f_model,
        variable.names = c("y_pred"),
        n.iter = 5000,
        thin = 10,
        progress.bar = "none"
    )
    f_pred <- summary(y_pred_samp)$quantiles[, 3]
    f_cv_mse_3[i] <- mean(((f_pred - y_test) * 100)^2)


    # ~ EVI2 based
    m_gsl_train <- as.vector(scale(north_sites_dt[year != yr]$m_gsl, scale = FALSE))
    m_emax_train <- as.vector(scale(north_sites_dt[year != yr]$m_EVImax_1, scale = FALSE)) * 10
    m_emin_train <- as.vector(scale(north_sites_dt[year != yr]$m_EVImin_1, scale = FALSE)) * 10

    mu_m_gsl <- attr(scale(north_sites_dt[year != yr]$m_gsl, scale = FALSE), "scaled:center")
    mu_m_max <- attr(scale(north_sites_dt[year != yr]$m_EVImax_1, scale = FALSE), "scaled:center")
    mu_m_min <- attr(scale(north_sites_dt[year != yr]$m_EVImin_1, scale = FALSE), "scaled:center")

    m_gsl_test <- north_sites_dt[year == yr, m_gsl] - mu_m_gsl
    m_max_test <- (north_sites_dt[year == yr, m_EVImax_1] - mu_m_max) * 10
    m_min_test <- (north_sites_dt[year == yr, m_EVImin_1] - mu_m_min) * 10


    m_model <- jags.model(textConnection(model_3),
        data = list(
            Y = y_train, n = n_train,
            gsl = m_gsl_train, emax = m_emax_train, emin = m_emin_train, igbp = igbp_train, nigbp = nigbp,
            ntest = n_test, gsl_test = m_gsl_test, emax_test = m_max_test, emin_test = m_min_test, igbp_test = igbp_test
        ),
        inits = list(
            beta1 = runif(nigbp, 10, 30), beta2 = runif(nigbp, 0, 20),
            beta3 = runif(nigbp, 10, 30), beta4 = runif(nigbp, 0, 30)
        ),
        n.chains = 2, quiet = TRUE
    )
    update(m_model, 2000, progress.bar = "none")
    repeat({
        samp <- coda.samples(m_model,
            variable.names = c("beta1", "beta2", "beta3", "beta4", "tau_y", "tau_beta", "mu_beta"),
            n.iter = 10000,
            thin = 10,
            progress.bar = "none"
        )
        mpsrf <- gelman.diag(samp)$mpsrf

        if (abs(mpsrf) < 1.1 & max(abs(gelman.diag(samp)$psrf)) < 1.1) break
    })
    y_pred_samp <- coda.samples(m_model,
        variable.names = c("y_pred"),
        n.iter = 5000,
        thin = 10,
        progress.bar = "none"
    )
    m_pred <- summary(y_pred_samp)$quantiles[, 3]
    m_cv_mse_3[i] <- mean(((m_pred - y_test) * 100)^2)
}
lines(uni_yrs, sqrt(f_cv_mse_3), type = "b", ylim = c(100, 500))
lines(uni_yrs, sqrt(m_cv_mse_3), type = "b", col = "blue")


# ~ Model 4
# ~~~~~~~~~~~~~~~~~
model_4 <- "model {
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

    # prediction
    for (i in 1:ntest) {
        y_pred[i] <- beta1[site_test[i]] + beta2[igbp_test[i]] * gsl_test[i] + beta3[igbp_test[i]] * emax_test[i] + beta4[igbp_test[i]] * emin_test[i]
    }
}"

f_cv_mse_4 <- rep(0, length(uni_yrs))
m_cv_mse_4 <- rep(0, length(uni_yrs))

for (i in 1:length(uni_yrs)) {
    yr <- uni_yrs[i]

    y_train <- north_sites_dt[year != yr]$annual_gpp / 100
    n_train <- nrow(north_sites_dt[year != yr])
    igbp_train <- f_igbp[north_sites_dt$year != yr]
    nigbp <- length(unique(f_igbp))
    site_train <- site[north_sites_dt$year != yr]
    
    y_test <- north_sites_dt[year == yr, annual_gpp / 100]
    n_test <- length(y_test)
    igbp_test <- f_igbp[north_sites_dt$year == yr]
    site_test <- site[north_sites_dt$year == yr]

    # ~ GPP based
    f_gsl_train <- as.vector(scale(north_sites_dt[year != yr]$f_gsl, scale = FALSE))
    f_emax_train <- as.vector(scale(north_sites_dt[year != yr]$f_gppmax_1, scale = FALSE))
    f_emin_train <- as.vector(scale(north_sites_dt[year != yr]$f_gppmin_1, scale = FALSE))

    mu_f_gsl <- attr(scale(north_sites_dt[year != yr]$f_gsl, scale = FALSE), "scaled:center")
    mu_f_max <- attr(scale(north_sites_dt[year != yr]$f_gppmax_1, scale = FALSE), "scaled:center")
    mu_f_min <- attr(scale(north_sites_dt[year != yr]$f_gppmin_1, scale = FALSE), "scaled:center")

    f_gsl_test <- north_sites_dt[year == yr, f_gsl] - mu_f_gsl
    f_max_test <- north_sites_dt[year == yr, f_gppmax_1] - mu_f_max
    f_min_test <- north_sites_dt[year == yr, f_gppmin_1] - mu_f_min


    f_model <- jags.model(textConnection(model_4),
        data = list(
            Y = y_train, n = n_train,
            gsl = f_gsl_train, emax = f_emax_train, emin = f_emin_train, igbp = igbp_train, nigbp = nigbp, site = site_train, ns = ns,
            ntest = n_test, gsl_test = f_gsl_test, emax_test = f_max_test, emin_test = f_min_test, igbp_test = igbp_test, site_test = site_test
        ),
        inits = list(
            beta1 = runif(ns, 10, 30), beta2 = runif(nigbp, 0, 20),
            beta3 = runif(nigbp, 10, 30), beta4 = runif(nigbp, 0, 30)
        ),
        n.chains = 2, quiet = TRUE
    )
    update(f_model, 2000, progress.bar = "none")
    repeat({
        samp <- coda.samples(f_model,
            variable.names = c("beta1", "beta2", "beta3", "beta4", "tau_y", "tau_beta", "mu_beta"),
            n.iter = 10000,
            thin = 10,
            progress.bar = "none"
        )
        mpsrf <- gelman.diag(samp)$mpsrf

        if (abs(mpsrf) < 1.1 & max(abs(gelman.diag(samp)$psrf)) < 1.1) break
    })
    y_pred_samp <- coda.samples(f_model,
        variable.names = c("y_pred"),
        n.iter = 5000,
        thin = 10,
        progress.bar = "none"
    )
    f_pred <- summary(y_pred_samp)$quantiles[, 3]
    f_cv_mse_4[i] <- mean(((f_pred - y_test) * 100)^2)


    # ~ EVI2 based
    m_gsl_train <- as.vector(scale(north_sites_dt[year != yr]$m_gsl, scale = FALSE))
    m_emax_train <- as.vector(scale(north_sites_dt[year != yr]$m_EVImax_1, scale = FALSE)) * 10
    m_emin_train <- as.vector(scale(north_sites_dt[year != yr]$m_EVImin_1, scale = FALSE)) * 10

    mu_m_gsl <- attr(scale(north_sites_dt[year != yr]$m_gsl, scale = FALSE), "scaled:center")
    mu_m_max <- attr(scale(north_sites_dt[year != yr]$m_EVImax_1, scale = FALSE), "scaled:center")
    mu_m_min <- attr(scale(north_sites_dt[year != yr]$m_EVImin_1, scale = FALSE), "scaled:center")

    m_gsl_test <- north_sites_dt[year == yr, m_gsl] - mu_m_gsl
    m_max_test <- (north_sites_dt[year == yr, m_EVImax_1] - mu_m_max) * 10
    m_min_test <- (north_sites_dt[year == yr, m_EVImin_1] - mu_m_min) * 10


    m_model <- jags.model(textConnection(model_4),
        data = list(
            Y = y_train, n = n_train,
            gsl = m_gsl_train, emax = m_emax_train, emin = m_emin_train, igbp = igbp_train, nigbp = nigbp, site = site_train, ns = ns,
            ntest = n_test, gsl_test = m_gsl_test, emax_test = m_max_test, emin_test = m_min_test, igbp_test = igbp_test, site_test = site_test
        ),
        inits = list(
            beta1 = runif(ns, 10, 30), beta2 = runif(nigbp, 0, 20),
            beta3 = runif(nigbp, 10, 30), beta4 = runif(nigbp, 0, 30)
        ),
        n.chains = 2, quiet = TRUE
    )
    update(m_model, 2000, progress.bar = "none")
    repeat({
        samp <- coda.samples(m_model,
            variable.names = c("beta1", "beta2", "beta3", "beta4", "tau_y", "tau_beta", "mu_beta"),
            n.iter = 10000,
            thin = 10,
            progress.bar = "none"
        )
        mpsrf <- gelman.diag(samp)$mpsrf

        if (abs(mpsrf) < 1.1 & max(abs(gelman.diag(samp)$psrf)) < 1.1) break
    })
    y_pred_samp <- coda.samples(m_model,
        variable.names = c("y_pred"),
        n.iter = 5000,
        thin = 10,
        progress.bar = "none"
    )
    m_pred <- summary(y_pred_samp)$quantiles[, 3]
    m_cv_mse_4[i] <- mean(((m_pred - y_test) * 100)^2)
}

lines(uni_yrs, sqrt(f_cv_mse_4), type = "b", ylim = c(100, 500))
lines(uni_yrs, sqrt(m_cv_mse_4), type = "b", col = "blue")



# ~ Model 5
# ~~~~~~~~~~~~~~~~~
model_5 <- "model {
    # likelihood
    for (i in 1:n) {
        Y[i] ~ dnorm(f[i], tau_y)
        f[i] <- beta1[site[i]] + beta2[site[i]] * gsl[i] + beta3[site[i]] * emax[i] + beta4[site[i]] * emin[i]
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

    # prediction
    for (i in 1:ntest) {
        y_pred[i] <- beta1[site_test[i]] + beta2[site_test[i]] * gsl_test[i] + beta3[site_test[i]] * emax_test[i] + beta4[site_test[i]] * emin_test[i]
    }
}"

f_cv_mse_5 <- rep(0, length(uni_yrs))
m_cv_mse_5 <- rep(0, length(uni_yrs))

for (i in 1:length(uni_yrs)) {
    yr <- uni_yrs[i]

    y_train <- north_sites_dt[year != yr]$annual_gpp / 100
    n_train <- nrow(north_sites_dt[year != yr])
    igbp_train <- m_igbp[north_sites_dt$year != yr]
    nigbp <- length(unique(m_igbp))
    site_train <- site[north_sites_dt$year != yr]
    site_igbp <- as.numeric(as.factor(unique(north_sites_dt[, .(site, IGBP)])$IGBP))

    y_test <- north_sites_dt[year == yr, annual_gpp / 100]
    n_test <- length(y_test)
    igbp_test <- m_igbp[north_sites_dt$year == yr]
    site_test <- site[north_sites_dt$year == yr]
    site_igbp_test <- as.numeric(as.factor(unique(north_sites_dt[year == yr, .(site, IGBP)])$IGBP))

    # ~ GPP based
    f_gsl_train <- as.vector(scale(north_sites_dt[year != yr]$f_gsl, scale = FALSE))
    f_emax_train <- as.vector(scale(north_sites_dt[year != yr]$f_gppmax_1, scale = FALSE))
    f_emin_train <- as.vector(scale(north_sites_dt[year != yr]$f_gppmin_1, scale = FALSE))

    mu_f_gsl <- attr(scale(north_sites_dt[year != yr]$f_gsl, scale = FALSE), "scaled:center")
    mu_f_max <- attr(scale(north_sites_dt[year != yr]$f_gppmax_1, scale = FALSE), "scaled:center")
    mu_f_min <- attr(scale(north_sites_dt[year != yr]$f_gppmin_1, scale = FALSE), "scaled:center")

    f_gsl_test <- north_sites_dt[year == yr, f_gsl] - mu_f_gsl
    f_max_test <- north_sites_dt[year == yr, f_gppmax_1] - mu_f_max
    f_min_test <- north_sites_dt[year == yr, f_gppmin_1] - mu_f_min

    f_model <- jags.model(textConnection(model_5),
        data = list(
            Y = y_train, n = n_train,
            gsl = f_gsl_train, emax = f_emax_train, emin = f_emin_train, nigbp = nigbp, site = site_train, ns = ns,
            ntest = n_test, gsl_test = f_gsl_test, emax_test = f_max_test, emin_test = f_min_test, site_test = site_test,
            site_igbp = site_igbp
        ),
        inits = list(
            beta1 = runif(ns, 10, 30), beta2 = runif(ns, 0, 20),
            beta3 = runif(ns, 10, 30), beta4 = runif(ns, 0, 30)
        ),
        n.chains = 2, quiet = TRUE
    )
    update(f_model, 2000, progress.bar = "none")
    iter <- 1
    repeat({
        samp <- coda.samples(f_model,
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
    y_pred_samp <- coda.samples(f_model,
        variable.names = c("y_pred"),
        n.iter = 5000,
        thin = 10,
        progress.bar = "none"
    )
    f_pred <- summary(y_pred_samp)$quantiles[, 3]
    f_cv_mse_5[i] <- mean(((f_pred - y_test) * 100)^2)
    

    # ~ EVI2 based
    m_gsl_train <- as.vector(scale(north_sites_dt[year != yr]$m_gsl, scale = FALSE))
    m_emax_train <- as.vector(scale(north_sites_dt[year != yr]$m_EVImax_1, scale = FALSE)) * 10
    m_emin_train <- as.vector(scale(north_sites_dt[year != yr]$m_EVImin_1, scale = FALSE)) * 10

    mu_m_gsl <- attr(scale(north_sites_dt[year != yr]$m_gsl, scale = FALSE), "scaled:center")
    mu_m_max <- attr(scale(north_sites_dt[year != yr]$m_EVImax_1, scale = FALSE), "scaled:center")
    mu_m_min <- attr(scale(north_sites_dt[year != yr]$m_EVImin_1, scale = FALSE), "scaled:center")

    m_gsl_test <- north_sites_dt[year == yr, m_gsl] - mu_m_gsl
    m_max_test <- (north_sites_dt[year == yr, m_EVImax_1] - mu_m_max) * 10
    m_min_test <- (north_sites_dt[year == yr, m_EVImin_1] - mu_m_min) * 10

    m_model <- jags.model(textConnection(model_5),
        data = list(
            Y = y_train, n = n_train,
            gsl = m_gsl_train, emax = m_emax_train, emin = m_emin_train, nigbp = nigbp, site = site_train, ns = ns,
            ntest = n_test, gsl_test = m_gsl_test, emax_test = m_max_test, emin_test = m_min_test, site_test = site_test,
            site_igbp = site_igbp
        ),
        inits = list(
            beta1 = runif(ns, 10, 30), beta2 = runif(ns, 0, 20),
            beta3 = runif(ns, 10, 30), beta4 = runif(ns, 0, 30)
        ),
        n.chains = 2, quiet = TRUE
    )
    update(m_model, 2000, progress.bar = "none")
    iter <- 1
    repeat({
        samp <- coda.samples(m_model,
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
    y_pred_samp <- coda.samples(m_model,
        variable.names = c("y_pred"),
        n.iter = 5000,
        thin = 10,
        progress.bar = "none"
    )
    m_pred <- summary(y_pred_samp)$quantiles[, 3]
    m_cv_mse_5[i] <- mean(((m_pred - y_test) * 100)^2)
}

lines(uni_yrs, sqrt(m_cv_mse_5), type = "b", col = "red")


# out: save cross validation results
save(f_cv_mse_1, m_cv_mse_1, f_cv_mse_2, m_cv_mse_2, f_cv_mse_3, m_cv_mse_3, 
    f_cv_mse_4, m_cv_mse_4, f_cv_mse_5, m_cv_mse_5, file = "Pipeline/cv.RData")



#fig: Maybe show the in boxplot
{
    png("Output/leave_year_cv.png", width = 1600, height = 1300, res = 300)
    par(mgp = c(1.5, 0.5, 0), mar = c(3, 5, 1, 1))
    boxplot(sqrt(f_cv_mse_1), sqrt(m_cv_mse_1), sqrt(f_cv_mse_2), sqrt(m_cv_mse_2), sqrt(f_cv_mse_3), sqrt(m_cv_mse_3), sqrt(
        f_cv_mse_4), sqrt(m_cv_mse_4), sqrt(m_f_cv_mse_4), col = c(rep(c("blue", "red"), 4),"red"), horizontal = TRUE, 
        yaxt = "n", xlab = "RMSE")

    sapply(seq(1, 7, by = 2), function(i) {
        x <- rep(c(0, 500), 4)
        y <- rep(seq(2.5, 8.5, by = 2), each = 2)
        df <- cbind(x, y)
        lines(df[c(i, i+1), 1], df[c(i, i+1), 2], lty = 2, col = "grey70")
    })
    axis(2, at = c(seq(1.5, 7.5, by = 2), 9), labels = c("Model 1", "Model 2", "Model 3", "Model 4", "Model 4*"), las = 2)
    legend("bottomright", legend = c("MODIS EVI2 based", "Flux GPP based"), fill = c("red", "blue"), bty = "n")
    dev.off()
}


