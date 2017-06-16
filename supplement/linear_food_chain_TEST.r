if(!require("tvvarss")) {
  devtools::install_github("nwfsc-timeseries/tvvarss")
  library("tvvarss")
}
if(!require("broom")) {
  install.packages("broom")
  library("broom")
}
if(!require("stan")) {
  install.packages("stan")
  library("stan")
}

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

n_species <- 4
n_year <- 60
n_site <- 1
n_simulations = 10
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
saved_output = vector("list", n_simulations)

## function for fitting/extracting
get_coefs <- function(y, topo, shared_r, mcmc_chain, mcmc_iter, mcmc_warmup,
                      intervals = TRUE, prob = 0.9) {
  fitted_model <- tvvarss(y, topo, shared_r, mcmc_chain, mcmc_iter, mcmc_warmup)
  coef <- tidy(fitted_model, intervals, prob)
}

for(ns in 1:n_simulations) {
  ## simulate process
  lfc <- simTVVAR(Bt = B0_init,
                  topo = B0_lfc,
                  TT = n_year,
                  var_QX = rev(seq(1, 4) / 40),
                  cov_QX = 0,
                  var_QB = 0.01,
                  cov_QB = 0)
  while(max(lfc$states) > dens_max | min(lfc$states) < dens_min) {
    lfc <- simTVVAR(Bt = B0_init,
                    topo = B0_lfc,
                    TT = n_year,
                    var_QX = rev(seq(1, 4) / 40),
                    cov_QX = 0,
                    var_QB = 0.01,
                    cov_QB = 0)
  }
  ## add obs error
  Y <- sim2fit(lfc, n_site, sd=0.1)
  ## fit model to only 2nd half of data & save param summaries
  coef <- get_coefs(y = Y[,-c(1:(n_year/2)),], topo = B0_lfc,
                    shared_r = matrix(1, n_species, n_site),
                    mcmc_chain = 3, mcmc_iter = 3000, mcmc_warmup = 2000)
  ## save data (Y), simulation output (lfc), model coefficients
  saved_output[[ns]] <- list('data' = Y, 'sim_output' = lfc, 'estimate' = coef)
}
## save results
save(saved_output, file = "output_linear_chain.Rdata")


## NOT RUN
## plot states
## lvls <- c("TC","SC","PC","PP")
## clr <- c("purple","darkred","blue","darkgreen")
## par(mai=c(0.9,0.9,0,0), omi=c(0.1,0.1,0.1,1.5))
## matplot(Y[,,], type="l", lty="solid", lwd=2, xpd=NA,
##         col=clr, ylab="Log density")
## legend("right", legend=lvls, lty="solid", lwd=2, bty="n",
##        col=clr, inset=-0.2, xpd=NA, cex=0.9)
