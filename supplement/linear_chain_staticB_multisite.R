library("tvvarss")
library("broom")
library("rstan")

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

n_species <- 4
n_year <- 60
n_site <- 5
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
  "sim" = 1:n_simulations, "site" = c(1,2,4))

## empty list for results
saved_output = vector("list", n_simulations)

## function for fitting/extracting
get_coefs <- function(y, topo, shared_r, shared_q, mcmc_chain, mcmc_iter, mcmc_warmup, process) {
  fitted_model <- tvvarss(y=y, topo=topo, shared_r=shared_r, shared_q=shared_q,
    mcmc_chain=mcmc_chain, mcmc_iter=mcmc_iter, mcmc_warmup=mcmc_warmup,
    dynamicB = FALSE, process = process)
  #coef <- tidy(fitted_model, intervals=TRUE, prob=0.9)
  return(fitted_model)
}

set.seed(123)
for(ns in 1:n_simulations) {

  B0_init <- matrix(0,n_species,n_species)
  diag(B0_init) = runif(n_species, c(0.2, 0.4, 0.6)[as.numeric(sim_config$B_diag[ns])],
    c(0.4, 0.7, 0.9)[as.numeric(sim_config$B_diag[ns])])
  B0_init[1,2] = runif(1, -0.1, 0)
  B0_init[2,1] = runif(1, 0, 0.1)
  B0_init[2,3] <- runif(1, -0.1, 0)
  B0_init[3,2] <- runif(1, 0, 0.1)
  B0_init[3,4] <- runif(1, -0.1, 0)
  B0_init[4,3] <- runif(1, 0, 0.1)

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
  Y <- sim2fit(lfc, n_sims = sim_config$site[ns], sd=0.01)

  # Burn-in the data, fitting model to 2nd half
  y_thin = array(0, dim=c(dim(Y)[1], dim(Y)[2]/2, dim(Y)[3]))
  for(ii in 1:dim(y_thin)[1]) {y_thin[ii,,] = Y[ii,-c(1:dim(Y)[2]/2),]}

  fitted <- get_coefs(y = y_thin, topo = B0_lfc,
    shared_r = matrix(1, n_species, dim(y_thin)[1]),
    shared_q = matrix(1, n_species, 1),
    process = rep(1, dim(y_thin)[1]),
    mcmc_chain = 1,
    mcmc_iter = 2000,
    mcmc_warmup = 1000)

  coef = tidy(fitted, conf.int = TRUE, conf.level=0.9,
    rhat = TRUE, ess = TRUE)
  coef = coef[which(startsWith(coef$term, "B[1,")),]
  ## save data (Y), simulation output (lfc), model coefficients
  saved_output[[ns]] <- list('data' = y_thin, 'sim_output' = lfc, 'estimate' = coef)
  print(ns)
  save(saved_output, file = "output_linear_chain_staticB.Rdata")
}
## save results


