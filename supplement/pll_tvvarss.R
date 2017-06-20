if(!require("tvvarss")) {
  devtools::install_github("nwfsc-timeseries/tvvarss")
  library("tvvarss")
}
if(!require("broom")) {
  install.packages("broom")
  library("broom")
}
if(!require("rstan")) {
  install.packages("rstan")
  library("rstan")
}
if(!require("foreach")) {
  install.packages("foreach")
  library("foreach")
}

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

## number of species/guilds
n_species <- 4
## number of years to simulate
nY <- 60
## number of sites
nS <- 1
## number of MC simulations
n_sims <- 4
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

## initial B
B0_init <- matrix(0,n_species,n_species)
diag(B0_init) <- 0.6
B0_init[1,2] <- -0.03
B0_init[2,1] <- 0.1
B0_init[2,3] <- -0.1
B0_init[3,2] <- 0.1
B0_init[3,4] <- -0.3
B0_init[4,3] <- 0.1

## empty list for results
saved_output <- vector("list", n_sims)

## define list of inputs for simulation
sim_list <- list(B0_init = B0_init,
                 B0_lfc = B0_lfc,
                 n_year = nY,
                 n_site = nS,
                 var_QX = rev(seq(1, 4) / 40),
                 cov_QX = 0,
                 var_QB = 0.01,
                 cov_QB = 0,
                 dens_max = dens_max,
                 dens_min = dens_min)

## define list of inputs for simulations
fit_list <- list(topo = B0_lfc,
                 shared_r = matrix(1, n_species, n_site),
                 mcmc_chain = 3,
                 mcmc_iter = 3000,
                 mcmc_warmup = 2000,
                 intervals = TRUE,
                 prob = 0.9)

## define function for simulating/fitting
simfit <- function(sim, fit) {
  ## extract list elements
  for(i in 1:length(sim)) {
    assign(names(sim[i]), sim[[i]], inherits = TRUE)
  }
  print(var_QX)
  for(i in 1:length(fit)) {
    assign(names(fit[i]), fit[[i]], inherits = TRUE)
  }
  ## simulate process
  lfc <- simTVVAR(Bt = B0_init, topo = B0_lfc, TT = n_year,
                  var_QX = var_QX, cov_QX = cov_QX,
                  var_QB = var_QB, cov_QB = cov_QB)
  while(max(lfc$states) > dens_max | min(lfc$states) < dens_min) {
    lfc <- simTVVAR(Bt = B0_init, topo = B0_lfc, TT = n_year,
                    var_QX = var_QX, cov_QX = cov_QX,
                    var_QB = var_QB, cov_QB = cov_QB)
  }
  ## add obs error
  Y <- sim2fit(lfc, n_site, sd=0.1)
  ## fit the model to first half of data
  fitted_model <- tvvarss(y=Y[,-c(1:(n_year/2)),], topo=topo, shared_r=shared_r,
                          mcmc_chain=mcmc_chain, mcmc_iter=mcmc_iter, mcmc_warmup=mcmc_warmup)
  coef <- tidy(fitted_model, intervals, prob)
  ## save data (Y), simulation output (lfc), model coefficients
  return(list('data' = Y, 'sim_output' = lfc, 'estimate' = coef))
}

saved_output <- times(n_sims) %dopar% simfit(sim_list, fit_list)

