library("tvvarss")
library("broom")
library("rstan")

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

n_species <- 4
n_year <- 60
n_simulations = 100
## min log-density threshold
dens_min <- -3
## max log-density threshold
dens_max <- 3

## topo matrix for linear food chain
B0_lfc <- matrix(list(0),n_species,n_species)
diag(B0_lfc) <- "dd"
for(i in 1:(n_species-1)) {
  B0_lfc[i,i+1] <- "td"
  B0_lfc[i+1,i] <- "bu"
}

sim_config = expand.grid("B_diag" = c("low", "med", "high"),
  "sim" = 1:n_simulations, "site" = c(1,2,3),
  "obs_error" = c("low", "med"))

## empty list for results
saved_output = vector("list", n_simulations)

set.seed(123)

for(ns in 1:nrow(sim_config)) {

  B0_init <- matrix(0,n_species,n_species)
  diag(B0_init) = runif(n_species, c(0.4, 0.6, 0.8)[as.numeric(sim_config$B_diag[ns])],
    c(0.6, 0.8, 1)[as.numeric(sim_config$B_diag[ns])])
  B0_init[1,2] = runif(1, -0.05, 0)
  B0_init[2,1] = runif(1, 0, 0.05)
  B0_init[2,3] <- runif(1, -0.05, 0)
  B0_init[3,2] <- runif(1, 0, 0.05)
  B0_init[3,4] <- runif(1, -0.05, 0)
  B0_init[4,3] <- runif(1, 0, 0.05)

  ## simulate process. var_QX is process error on states.
  ## var_QB is process var on B -- ignored for static B models.
  lfc <- simTVVAR(Bt = B0_init,
    topo = B0_lfc,
    TT = n_year,
    var_QX = 0.02^2,
    cov_QX = 0,
    var_QB = 0,
    cov_QB = 0)

  while(max(lfc$states) > dens_max | min(lfc$states) < dens_min) {
    ## simulate process. var_QX is process error on states.
    ## var_QB is process var on B -- ignored for static B models.
    lfc <- simTVVAR(Bt = B0_init,
      topo = B0_lfc,
      TT = n_year,
      var_QX = 0.02^2,
      cov_QX = 0,
      var_QB = 0,
      cov_QB = 0)
  }

  ## add obs error
  Y <- sim2fit(lfc, n_sims = sim_config$site[ns],
    sd=c(0.01,0.1)[as.numeric(sim_config$obs_error[ns])])

  # Burn-in the data, fitting model to 2nd half
  y_thin = array(0, dim=c(dim(Y)[1], dim(Y)[2]/2, dim(Y)[3]))
  for(ii in 1:dim(y_thin)[1]) {y_thin[ii,,] = Y[ii,-c(1:dim(Y)[2]/2),]}

  fitted_model <- tvvarss(y=y_thin, topo=B0_lfc,
    shared_r=matrix(1, n_species, dim(y_thin)[1]),
    shared_q=matrix(1, n_species, 1),
    mcmc_chain=1, mcmc_iter=400,
    mcmc_warmup=100, mcmc_thin = 1,
    dynamicB = FALSE, process = rep(1, dim(y_thin)[1]))

  coefs = tidy(fitted_model, conf.int = TRUE, conf.level=0.9,
    rhat = TRUE, ess = TRUE)

  coef = coefs[which(startsWith(coefs$term, "B[1,")),]
  coef = rbind(coef, coefs[which(startsWith(coefs$term, "resid_process_sd")),])
  coef = rbind(coef, coefs[which(startsWith(coefs$term, "obs_sd")),])
  ## save data (Y), simulation output (lfc), model coefficients
  ## add check for convergence
  saved_output[[ns]] <- list('data' = y_thin, 'sim_output' = lfc, 'estimate' = coef,
    'converged' = ifelse(max(coef$rhat, na.rm=T) < 1.1, 1, 0))
  print(ns)
  save(saved_output, file = "linear_chain_staticB_multisite.Rdata")
}

## save results


