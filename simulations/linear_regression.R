library(INLA)
library(INLABMA)

# Load helper functions
source("../R/Data_Generation.R")
source("../R/INLAMH_Utils.R")
source("../R/Post_Processing.R")

# Simulate data
dat <- sim_linear_regression(n=1000, p=10, num_zero=3, sd=3, seed=72)

# Fit model using INLA w/in MCMC
burn <- 100 # Number of burn in MCMC iterations
iter <- 2000 # Number of MCMC iterations after burn in  
scale <- 0.25 # Standard deviation for M-H proposals
set.seed(80) # Set random seed before fitting model
time <- system.time(
    results <- INLAMH(d=dat, fit.inla=hs_linear_regression_inla,
                      b.init=rep(0, dat$p+1), rq=proposal_naive,
                      dq=proposal_density_naive, prior=prior_density, 
                      n.burnin=burn, n.sim=iter, verbose=TRUE)
)

# Process results
accepts <- results$acc.sim
samps <- process_MCMC(dat, iter, results$b.sim)
marginals <- process_INLA(dat, iter, results$model.sim, hyper=TRUE)

# Plot results
plot_MCMC_trace(dat, iter=iter, samps)
plot_MCMC_hyperpar(dat, samps)
plot_INLA_beta(dat, marginals)
plot_INLA_hyperpar(dat, marginals, normal_precision=TRUE)

# Save full results
save(dat, iter, scale, time, results, accepts, samps, marginals,
     file="results/linear_regression_results_full.RData")

# Save only processed results
save(dat, iter, scale, time, accepts, samps, marginals, 
     file="results/linear_regression_results_processed.RData")

