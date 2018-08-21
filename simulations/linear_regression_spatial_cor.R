library(INLA)
library(INLABMA)

# Load helper functions
source("../R/Data_Generation.R")
source("../R/INLAMH_Utils.R")
source("../R/Post_Processing.R")

# Simulate data
loc <- sim_locations(n=1000, lim=c(0,1), seed=70)
dat <- sim_linear_regression_spatial(loc, p=10, num_zero=3, sd=3,
                                     sre.sigma=1, sre.rho=0.5, 
                                     cor_matrix=diag(0.9, 10)+
                                         matrix(rep(0.1, 100), nrow=10),
                                     seed=72)
plot(dat$loc)

# Create mesh for SPDE and specify PC priors
bound <- inla.nonconvex.hull(dat$loc, convex=0.25)
mesh <- inla.mesh.2d(dat$loc, boundary=bound, max.edge=c(0.1, 0.25), 
                     cutoff=0.05)
plot(mesh, asp=1)
points(dat$loc, pch=19, col=2)
A <- inla.spde.make.A(mesh, dat$loc)
spde <- inla.spde2.pcmatern(mesh, alpha=2, 
                            prior.range=c(0.1, 0.05), # P(rho<0.1)=0.05
                            prior.sigma=c(5, 0.05)) # P(sigma>5)=0.05
stk <- inla.stack(data=data.frame(y=dat$y), A=list(A, 1),
                  effect=list(spatial=1:spde$n.spde,
                              cbind(dat$X, data.frame(int=rep(1, dat$n)))))
dat$stk <- stk
dat$spde <- spde

# Fit model using INLA w/in MCMC
burn <- 100 # Number of burn in MCMC iterations
iter <- 2000 # Number of MCMC iterations after burn in  
scale <- 0.25 # Standard deviation for M-H proposals
set.seed(80) # Set random seed before fitting model
time <- system.time(
    results <- INLAMH(d=dat, fit.inla=hs_linear_regression_spatial_inla,
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
plot_INLA_hyperpar(dat, marginals, normal_precision=TRUE, spatial=TRUE,
                   layout=c(2,2))

# Save full results
save(dat, bound, mesh, A, iter, scale, time, results, accepts, samps, marginals,
     file="results/linear_regression_spatial_cor_results_full.RData")

# Save only processed results
save(dat, bound, mesh, A, iter, scale, time, accepts, samps, marginals,
     file="results/linear_regression_spatial_cor_results_processed.RData")


