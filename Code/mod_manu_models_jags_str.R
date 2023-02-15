#*******************************************************************************
# Description: This file contains all the model structures selected for 
# the manuscript.
#
# Author: Xiaojie Gao
# Date: 2022-01-06
#*******************************************************************************

# The simplest model
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
model_1_params <- c("beta1", "beta2", "beta3", "beta4", "tau_y")



# Model with biome-level intercept
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
}"
model_2_params <- c("beta1", "beta2", "beta3", "beta4", "tau_y", "tau_beta", "mu_beta")



# Model with biome-level intercepts and slopes
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
}"
model_3_params <- c("beta1", "beta2", "beta3", "beta4", "tau_y", "tau_beta", "mu_beta")



# Model with biome-level slopes and site-level intercepts
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
}"
model_4_params <- c("beta1", "beta2", "beta3", "beta4", "tau_y", "tau_beta", "mu_beta")


# Model with site-level intercepts and slopes
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
}"
model_5_params <- c("beta1", "beta2", "beta3", "beta4",
    "eta1", "eta2", "eta3", "eta4",
    "sigma_y", "sigma_beta", "sigma_eta", "mu"
)