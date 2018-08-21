library(fields) # To generate covariance matrix for simulated GPs
library(mvtnorm) # To simulate GPs w/ a given covariance matrix

# Simulation of spatial locations
# n: Number of spatial locations to simulate
# lim: Spatial locations are simulated uniformly on the cube lim x lim
# seed: Random seed to be set before generating spatial locations
sim_locations <- function(n=100, lim=c(0, 1), seed=72){
    set.seed(seed)
    return(matrix(runif(2*n, min=lim[1], max=lim[2]), ncol=2))
}

#### Simulate various types of glm data w/ and w/o spatial random effects
#### Common Inputs
# n: Number of data points to simulate
# loc: Spatial locations of data points to be simulated
# p: Number of covariates to simulate (not including an intercept)
# num_zero: Number of covariates that are set to 0 (<=p)
# sd: Standard deviation for likelihood of normal models
# n_trial: Number of trials for binomial models
# covariate.sigma: Margincal variance of GP for simulated spatial covariates
# covariate.rho: Spatial range of GP for simulated spatial covariates
# sre.sigma: Margincal variance of GP for simulated spatial random effects
# sre.rho: Spatial range of GP for simulated spatial random effects
# cor_matrix: Matrix that describes the correlation between covariates
# seed: Random seed to be set before generating data

# Simulation for normal linear regression data
sim_linear_regression <- function(n=100, p=10, num_zero=3, sd=3, seed=72){
    set.seed(seed)
    # Design matrix is standard normal
    X <- cbind(rep(1, n), matrix(rnorm(n*p), nrow=n, ncol=p))
    # Regression coefficients are standard t w/ df=1, i.e. Cauchy
    beta <- rt(p+1, df=1)
    # Zero out the last num_zero coefficients
    beta[(p+2-num_zero):(p+1)] <- 0
    # Generate responses
    noise <- rnorm(n, sd=sd)
    y <- X %*% beta + noise
    return(list(n=n, p=p, y=y, X=data.frame(X[, 2:(p+1)]), beta=beta, sd=sd))
}

# Simulation for binomial regression data
sim_binomial_regression <- function(n=100, p=10, num_zero=3, n_trial=1, 
                                    seed=72){
    set.seed(seed)
    # Design matrix is standard normal
    X <- cbind(rep(1, n), matrix(rnorm(n*p), nrow=n, ncol=p))
    # Regression coefficients are standard t w/ df=1, i.e. Cauchy
    beta <- rt(p+1, df=1)
    # Zero out the last num_zero coefficients
    beta[(p+2-num_zero):(p+1)] <- 0
    # Generate responses
    probs <- 1/(1+exp(-(X %*% beta)))
    y <- rbinom(n, size=n_trial, prob=probs)
    return(list(n=n, p=p, y=y, X=data.frame(X[, 2:(p+1)]), n_trial=n_trial,
                beta=beta))
}

# Simulation for normal linear regression data w/ spatial covariates and 
# spatial random effects
sim_linear_regression_spatial <- function(loc, p=10, num_zero=3, sd=3,
                                          covariate.sigma=1, covariate.rho=1,
                                          sre.sigma=1, sre.rho=1, 
                                          cor_matrix=diag(p), seed=72){
    set.seed(seed)
    n <- nrow(loc)
    # Generate covariates from a GP w/ Matern covariance 
    # nu=1, sigma^2=covariate.sigma, rho=covariate.rho
    X <- matrix(rep(0, n*p), nrow=n, ncol=p)
    for(i in 1:p){
        temp_cov <- stationary.cov(loc, loc, theta=1, Distance="rdist",
                                   Covariance="Matern", nu=1, 
                                   phi=covariate.sigma,
                                   alpha=sqrt(8)/covariate.rho)
        X[, i] <- c(rmvnorm(1, sigma=temp_cov))
    }
    # Use cor_matrix to make covariates correlated
    X <- t(chol(cor_matrix) %*% t(X))
    X <- cbind(rep(1, n), X)
    # Regression coefficients are standard t w/ df=1, i.e. Cauchy
    beta <- rt(p+1, df=1)
    # Zero out the last num_zero coefficients
    beta[(p+2-num_zero):(p+1)] <- 0
    # Generate spatial random effects from a GP w/ Matern covariance
    # nu=1, sigma^2=sre.sigma, rho=sre.rho
    sre_cov <- stationary.cov(loc, loc, theta=1, Distance= "rdist",
                              Covariance="Matern", nu=1, phi=sre.sigma,
                              alpha=sqrt(8)/sre.rho)
    sre <- c(rmvnorm(1, sigma=sre_cov))
    # Generate responses
    nugget <- rnorm(n, sd=sd)
    y <- X %*% beta + sre + nugget
    return(list(loc=loc, n=n, p=p, y=y, X=data.frame(X[, 2:(p+1)]), beta=beta, 
                covariate.sigma=covariate.sigma, covariate.rho=covariate.rho,
                sd=sd, sre.sigma=sre.sigma, sre.rho=sre.rho,
                # Include the generated sre and nugget for debugging purposes
                sre=sre, nugget=nugget))
}

# Simulation for binomial regression data w/ spatial covariates and 
# spatial random effects
sim_binomial_regression_spatial <- function(loc, p=10, num_zero=3, n_trial=1, 
                                            covariate.sigma=1, covariate.rho=1,
                                            sre.sigma=1, sre.rho=1, 
                                            cor_matrix=diag(p), 
                                            seed=72){
    set.seed(seed)
    n <- nrow(loc)
    # Generate covariates from a GP w/ Matern covariance 
    # nu=1, sigma^2=covariate.sigma, rho=covariate.rho
    X <- matrix(rep(0, n*p), nrow=n, ncol=p)
    for(i in 1:p){
        temp_cov <- stationary.cov(loc, loc, theta=1, Distance="rdist",
                                   Covariance="Matern", nu=1, 
                                   phi=covariate.sigma,
                                   alpha=sqrt(8)/covariate.rho)
        X[, i] <- c(rmvnorm(1, sigma=temp_cov))
    }
    # Use cor_matrix to make covariates correlated
    X <- t(chol(cor_matrix) %*% t(X))
    X <- cbind(rep(1, n), X)
    # Regression coefficients are standard t w/ df=1, i.e. Cauchy
    beta <- rt(p+1, df=1)
    # Zero out the last num_zero coefficients
    beta[(p+2-num_zero):(p+1)] <- 0
    # Generate spatial random effects from a GP w/ Matern covariance
    # nu=1, sigma^2=sre.sigma, rho=sre.rho
    sre_cov <- stationary.cov(loc, loc, theta=1, Distance= "rdist",
                              Covariance="Matern", nu=1, phi=sre.sigma,
                              alpha=sqrt(8)/sre.rho)
    sre <- c(rmvnorm(1, sigma=sre_cov))
    # Generate responses
    probs <- 1/(1+exp(-(X %*% beta + sre)))
    y <- rbinom(n, size=n_trial, prob=probs)
    return(list(loc=loc, n=n, p=p, y=y, X=data.frame(X[, 2:(p+1)]), beta=beta,
                covariate.sigma=covariate.sigma, covariate.rho=covariate.rho,
                sre.sigma=sre.sigma, sre.rho=sre.rho, n_trial=n_trial,
                # Include the generated sre for debugging purposes
                sre=sre))
}












