library(INLABMA) # To average INLA posteriors

#### Functions to post-process the INLA and MCMC results
#### Common Inputs
# dat: data used to fit models
# iter: number of MCMC iterations

# Transform list of MCMC samples of log parameters given by INLAMH output into 
# matrix of samples of parameters
# samps: List of MCMC samples from INLAMH
process_MCMC <- function(dat, iter, samps){
    samps_temp <- matrix(rep(0, iter*(dat$p+1)), ncol=(dat$p+1), nrow=iter)
    for(i in 1:iter){ samps_temp[i, ] <- exp(samps[[i]]) }
    return(samps_temp)
}

# Transform list of INLA results given by INLAMH output into posterior marginal 
# densities by averaging over the posterior marginals at each iteration
# inla_results: List of fit INLA models from INLAMH
# hyper: If TRUE, retrieve posterior marginals for hyperparameters
process_INLA <- function(dat, iter, inla_results, hyper=TRUE){
    clean <- vector("list", iter)
    for(i in 1:iter){ clean[[i]] <- inla_results[[i]][[1]] }
    bma <- INLABMA:::fitmargBMA2(clean, ws=rep(1/iter,iter), 
                                 item=c("marginals.fixed"))
    if(hyper){
        bma_hyper <- INLABMA:::fitmargBMA2(clean, ws=rep(1/iter,iter),
                                           item=c("marginals.hyperpar"))
        bma <- append(bma, bma_hyper)
    }
    return(bma)
}

#### Functions to plot INLA and MCMC results
#### Common Inputs
# dat: data used to fit models
# iter: number of MCMC iterations
# samps: Processed MCMC samples
# marginals: Processed posterior marginals
# layout: Size of grid to plot posteriors on
# stan: If TRUE, also plot posteriors obtained from Stan
# stan_samps: Posterior samples obtained from Stan

# Plot traceplots of log tau and log lambda MCMC samples
plot_MCMC_trace <- function(dat, iter, samps, layout=c(3, 4)){
    par(mfrow=layout)
    plot(1:iter, log(samps[1:iter, 1]), type="l", main="Log Tau", ylab="val", 
         xlab="Iteration")
    for(i in 1:dat$p){
        plot(1:iter, log(samps[1:iter, 1+i]), type="l", 
             main=paste("Log Lambda", i), ylab="val", xlab="Iteration")
    }
    par(mfrow=c(1,1))
}

# Plot densities of tau and lambda MCMC samples 
plot_MCMC_hyperpar <- function(dat, samps, layout=c(3, 4), 
                               stan=FALSE, stan_samps=NULL){
    par(mfrow=layout)
    plot(density(samps[, 1]), col="red", xlab="val", main="Tau", ylab="Density")
    if(stan){ lines(density(stan_samps$tau)) }
    for(i in 1:dat$p){
        plot(density(samps[, i+1]), col="red", xlab="val", ylab="Density",
             main=paste("Lambda", i))
        if(stan){ lines(density(stan_samps$lambda[, i])) }
    }
    par(mfrow=c(1,1))
}

# Plot posterior marginals for regression coefficients from INLA
plot_INLA_beta <- function(dat, marginals, layout=c(3, 4), 
                           stan=FALSE, stan_samps=NULL){
    par(mfrow=layout)
    plot(marginals[[1]], col="red", xlab="val", main="Beta 0", ylab="Density",
         type="l")
    if(stan){ lines(density(stan_samps$beta_0)) }
    if(!is.na(dat$beta)){ abline(v=dat$beta[1]) }
    for(i in 1:dat$p){
        plot(marginals[[i+1]], col="red", xlab="val", main=paste("Beta", i),
             type="l", ylab="Density")
        if(stan){ lines(density(stan_samps$beta[, i])) }
        if(!is.na(dat$beta)){ abline(v=dat$beta[i+1]) }
        abline(v=0, col="gray47", lty="dotted", lwd=1.5)
    }
    par(mfrow=c(1,1))
}

# Plot posterior marginals for hyperparameters from INLA
# normal_precision: If TRUE, plot posterior marginal for precision of normal
# likelihood
# spatial: If TRUE, plot posterior marginals for spatial hyperparameters
# normal_random_effect: If TRUE, plot posterior marginals for precision of 
# normal random effect
plot_INLA_hyperpar <- function(dat, marginals, layout=c(1, 1), 
                               normal_precision=FALSE, spatial=FALSE,
                               normal_random_effect=FALSE,
                               stan=FALSE, stan_samps=NULL){
    par(mfrow=layout)
    if(normal_precision){
        plot(marginals$`Precision for the Gaussian observations`, col="red", 
             xlab="val", main="Precision", ylab="Density",
             type="l")
        if(stan){lines(density(stan_samps$precision))}
        if(!is.na(dat$sd)){ abline(v=1/(dat$sd^2)) }
    }
    if(spatial){
        plot(marginals$`Range for spatial`, col="red", xlab="val", 
             main="Spatial Range", type="l", ylab="Density")
        if(!is.na(dat$sre.rho)){ abline(v=dat$sre.rho)}
        plot(marginals$`Stdev for spatial`, col="red", xlab="val", 
             main="Spatial Standard Deviation", type="l", ylab="Density")
        if(!is.na(dat$sre.sigma)){abline(v=sqrt(dat$sre.sigma))}
    }
    if(normal_random_effect){
        plot(marginals$`Precision for z`, col="red", xlab="val", 
             main="Random Effect Precision", type="l", ylab="Density")
    }
    par(mfrow=c(1,1))
}



