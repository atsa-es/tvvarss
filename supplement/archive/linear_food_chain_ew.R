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

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

n_species <- 2
n_year <- 50
n_site <- 2
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

## initial B
B0_init <- matrix(0,n_species,n_species)
diag(B0_init) <- 0.6
B0_init[1,2] <- -0.03
B0_init[2,1] <- 0.1

set.seed(1)
## simulate process
  lfc <- simTVVAR(Bt = B0_init,
                  topo = B0_lfc,
                  TT = n_year,
                  var_QX = 0.01,
                  cov_QX = 0,
                  var_QB = 0,
                  cov_QB = 0)

  ## add obs error
  Y <- sim2fit(lfc, n_site, sd=0.000001)
  lfc$states = lfc$states[,-1]
  Y[1,,] = lfc$states + matrix(rnorm(dim(lfc$states)[1]*dim(lfc$states)[2],0,0.001), dim(lfc$states)[1], dim(lfc$states)[2])
  Y[2,,] = lfc$states + matrix(rnorm(dim(lfc$states)[1]*dim(lfc$states)[2],0,0.001), dim(lfc$states)[1], dim(lfc$states)[2])

  # drop burn-in
  yy = Y[,-c(1:10),]

  ## fit model to only 2nd half of data & save param summaries
  #marss_constrained <- tvvarss(y=yy, topo=B0_lfc, shared_r = matrix(1, n_species, n_site),
  #  mcmc_chain=1, mcmc_iter=2000, mcmc_warmup=1000, dynamicB = FALSE, site_indices = c(1,1))
  #apply(extract(marss_constrained, "B")$B, c(2,3,4), mean)[1,,]

  marss_constrained <- tvvarss(y=yy, topo=B0_lfc, shared_r = matrix(1, n_species, n_site),
    mcmc_chain=3, mcmc_iter=1000, mcmc_warmup=500, dynamicB = TRUE)
  apply(extract(marss_constrained, "B")$B, c(2,3,4), mean)[1,,]

# Change up constraints on B matrix
  B0_lfc[1,2] = "cf"
  B0_lfc[2,1] = "cf"
  ## fit model to only 2nd half of data & save param summaries
  marss_unconstrained <- tvvarss(y=yy, topo=B0_lfc, shared_r = matrix(1, n_species, n_site),
    mcmc_chain=3, mcmc_iter=3000, mcmc_warmup=2000, dynamicB = FALSE)
  apply(extract(marss_unconstrained, "B")$B, c(2,3,4), mean)[1,,]

# Try MARSS
  library(MARSS)
  MARSS(t(yy), model = list("B"="unconstrained","U"="zero"))
