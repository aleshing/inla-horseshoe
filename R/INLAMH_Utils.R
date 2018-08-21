#### Helper functions for MCMC to be passed to INLAMH

# Log density of a log transformed half-Cauchy random variable
# val: Value to evaluate the prior density at
log_half_cauchy_density <- function(val){
    return(val+log(2/pi)-log(1+exp(2*val)))
}
log_half_cauchy_density <- Vectorize(log_half_cauchy_density)

# Log density of multiple log transformed half-Cauchy random variables
# samps: Samples to evaluate the prior density at
prior_density <- function(samps){
    return(sum(log_half_cauchy_density(samps)))
}

# Naive Metropolis-Hastings proposal from a normal distribution
# samps: Samples which parameterize the mean of the M-H proposal
# s: Standard deviation of the M-H proposal
proposal_naive <- function(samps, s=scale){
    return(rnorm(length(samps), mean=samps, sd=s))
}

# Log density of a Metropolis-Hastings proposal from a normal distribution
# curr_samps: Samples which parameterize the mean of the M-H proposal
# prop_samps: Proposed M-H samples
# s: Standard deviation of the M-H proposal
proposal_density_naive <- function(curr_samps, prop_samps, s=scale){
    return(sum(dnorm(prop_samps, mean=curr_samps, sd=s, log=TRUE)))
}

#### Helper functions for INLA to be passed to INLAMH
#### Given samples for parameters of Horseshoe not fittable in INLA, these
#### functions fit the given models condtional on these samples in INLA
#### Currently use simplified.laplace, need to compare to Gaussian and Laplace
#### Common Inputs
# dat: Data needed to fit models in INLA
# samps: Proposed M-H samples

# Normal linear regression 
# Use default INLA prior for precision of normal likelihood
hs_linear_regression_inla <- function(dat, samps){
    # Global shrinkage parameter tau is the first sample
    tau=exp(samps[1])
    lambda=exp(samps[2:length(samps)])
    inla.dat <- cbind(dat$X, data.frame(y=dat$y))
    # Build up formula using MCMC samples
    formula <- "y~1"
    for(i in 1:dat$p){
        formula <- paste(formula, "+", "f(", colnames(dat$X)[i], 
                         ", model='linear', mean.linear=", 0, ", prec.linear=",
                         1/(tau*lambda[i])^2, ")", sep="")
    }
    # Fit model in INLA
    result <- inla(formula(formula), data=inla.dat, family="Gaussian",
                   control.fixed=list(mean.intercept=0, prec.intercept=1),
                   control.compute=list(mlik=TRUE),
                   control.inla=list(strategy="simplified.laplace"))
    return(list(model.sim=result, mlik=result$mlik))
}

# Binomial regression 
hs_binomial_regression_inla <- function(dat, samps){
    # Global shrinkage parameter tau is the first sample
    tau=exp(samps[1])
    lambda=exp(samps[2:length(samps)])
    inla.dat <- cbind(dat$X, data.frame(y=dat$y))
    # Build up formula using MCMC samples
    formula <- "y~1"
    for(i in 1:dat$p){
        formula <- paste(formula, "+", "f(", colnames(dat$X)[i], 
                         ", model='linear', mean.linear=", 0, ", prec.linear=", 
                         1/(tau*lambda[i])^2, ")", sep="")
    }
    # Fit model in INLA
    result <- inla(formula(formula), data=inla.dat, family="Binomial",
                   Ntrials=dat$n_trial,
                   control.fixed=list(mean.intercept=0, prec.intercept=1),
                   control.compute=list(mlik=TRUE),
                   control.inla=list(strategy="simplified.laplace"))
    return(list(model.sim=result, mlik=result$mlik))
}

# Normal linear regression w/ spatial random effects
# Use default INLA prior for precision of normal likelihood
hs_linear_regression_spatial_inla <- function(dat, samps){
    # Global shrinkage parameter tau is the first sample
    tau=exp(samps[1])
    lambda=exp(samps[2:length(samps)])
    # Build up formula using MCMC samples
    # Assume that the SPDE approximation for the spatial random effects
    # has already been built up, and is passed in through dat w/ stk and spde
    formula <- "y~0+f(int, model='linear', mean.linear=0, prec.linear= 1)"
    for(i in 1:dat$p){
        formula <- paste(formula, "+", "f(", colnames(dat$X)[i], 
                         ", model='linear', mean.linear=", 0, ", prec.linear=",
                         1/(tau*lambda[i])^2, ")", sep="")
    }
    formula <- paste(formula, "+ f(spatial, model=dat$spde)", sep="")
    # Fit model in INLA
    result <- inla(formula(formula), family="Gaussian", 
                   data=inla.stack.data(dat$stk),
                   control.compute=list(config=TRUE, mlik=TRUE),
                   control.inla=list(strategy="simplified.laplace"),
                   control.predictor=list(A=inla.stack.A(dat$stk)))
    return(list(model.sim=result, mlik=result$mlik))
}

# Binomial regression w/ spatial random effects
hs_binomial_regression_spatial_inla <- function(dat, samps){
    # Global shrinkage parameter tau is the first sample
    tau=exp(samps[1])
    lambda=exp(samps[2:length(samps)])
    # Build up formula using MCMC samples
    # Assume that the SPDE approximation for the spatial random effects
    # has already been built up, and is passed in through dat w/ stk and spde
    formula <- "y~0+f(int, model='linear', mean.linear=0, prec.linear= 1)"
    for(i in 1:dat$p){
        formula <- paste(formula, "+", "f(", colnames(dat$X)[i], 
                         ", model='linear', mean.linear=", 0, ", prec.linear=", 
                         1/(tau*lambda[i])^2, ")", sep="")
    }
    formula <- paste(formula, "+ f(spatial, model=dat$spde)", sep="")
    # Fit model in INLA
    result <- inla(formula(formula), family="Binomial",
                   Ntrials=dat$n_trial,
                   data=inla.stack.data(dat$stk),
                   control.compute=list(mlik=TRUE),
                   control.inla=list(strategy="simplified.laplace"),
                   control.predictor=list(A=inla.stack.A(dat$stk)))
    return(list(model.sim=result, mlik=result$mlik))
}

# Binomial regression w/ normal random effects (non-spatial)
hs_binomial_regression_random_effects_inla <- function(dat, samps){
    # Global shrinkage parameter tau is the first sample
    tau=exp(samps[1])
    lambda=exp(samps[2:length(samps)])
    inla.dat <- cbind(dat$X, data.frame(y=dat$y), data.frame(z=1:dat$n))
    # Build up formula using MCMC samples
    formula <- "y~1"
    for(i in 1:dat$p){
        formula <- paste(formula, "+", "f(", colnames(dat$X)[i], 
                         ", model='linear', mean.linear=", 0, ", prec.linear=", 
                         1/(tau*lambda[i])^2, ")", sep="")
    }
    formula <- paste(formula, "+ f(z, model=\"iid\")", sep="")
    # Fit model in INLA
    result <- inla(formula(formula), family="Binomial",
                   Ntrials=dat$n_trial,
                   data=inla.dat,
                   control.fixed=list(mean.intercept=0, prec.intercept=1),
                   control.compute=list(mlik=TRUE),
                   control.inla=list(strategy="simplified.laplace"))
    return(list(model.sim=result, mlik=result$mlik))
}







